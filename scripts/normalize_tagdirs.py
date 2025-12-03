#!/usr/bin/env python3
import os
import shutil
import subprocess
import argparse
import pandas as pd
from pathlib import Path

def normalize_tagdirs(metadata_file: Path, genome_dirs: dict):
    # genome_dirs keys = target + spike species
    # assume the first key is target species, rest are spike species
    species_list = list(genome_dirs.keys())
    target_species = species_list[0]
    spike_species = species_list[1:]
    """
    Normalize HOMER tag directories of the target species based on spike-in factors.

    Parameters
    ----------
    metadata_file : Path
        Path to sample_metadata.norm.tsv
    genome_dirs : dict
        Dictionary of genome -> Path to genome_data directory
    target_species : str
        Target genome whose tagdirs should be normalized
    spike_species : list
        List of spike genomes (e.g., ["dm6", "sac3"]) to use for normalization factors
    """

    seqstats_file = Path(metadata_file).resolve()
    seqstats = pd.read_csv(seqstats_file, sep="\t")

    # Clean column names
    seqstats.columns = (
        seqstats.columns
        .str.replace(r"[\u200b\u00a0\r\n]", "", regex=True)
        .str.strip()
        .str.replace(" ", ".", regex=False)
    )

    print("\nAvailable columns:", seqstats.columns.tolist())

    # Directories
    target_dir = genome_dirs[target_species]
    aligned_dir = target_dir / f"{target_species}_aligned"
    base_dir = target_dir / f"{target_species}_tagdirs"
    output_dir = target_dir / f"{target_species}_normalized_tagdirs"

    os.makedirs(base_dir, exist_ok=True)
    os.makedirs(output_dir, exist_ok=True)

    suffix_to_remove = f".{target_species}-tagdir"

    print(f"\nProcessing target species: {target_species}")
    print("Tagdirs will be written to:", base_dir)
    print("Normalized tagdirs will be written to:", output_dir)
    print("Looking for SAM files in:", aligned_dir)

    # ================= Step 1: Make tag directories =================
    for sam_file in aligned_dir.glob(f"*.{target_species}.nosuffx2.sam"):
        sample_id = sam_file.name.replace(f".{target_species}.nosuffx2.sam", "")
        tagdir_name = f"{sample_id}.{target_species}-tagdir"
        tagdir_path = base_dir / tagdir_name

        if tagdir_path.exists():
            print(f" Tagdir exists for {sample_id}, skipping makeTagDirectory.")
            continue

        print(f" Creating tag directory for {sample_id}...")
        cmd = [
            "makeTagDirectory",
            str(tagdir_path),
            str(sam_file),
            "-genome", target_species,
            "-checkGC",
            "-fragLength", "150",
        ]
        try:
            subprocess.run(cmd, check=True)
            print(" Created:", tagdir_name)
        except subprocess.CalledProcessError as e:
            print(f" Failed for {sample_id}: {e}")

    # ================= Step 2: Normalization =================
    for tagdir in base_dir.iterdir():
        if not (tagdir.is_dir() and tagdir.name.endswith(suffix_to_remove)):
            continue

        sample_id = tagdir.name.replace(suffix_to_remove, "")
        if sample_id not in seqstats["library.ID"].values:
            print(f" Sample {sample_id} not found in metadata, skipping.")
            continue

        row = seqstats.loc[seqstats["library.ID"] == sample_id].iloc[0]

        # Normalize once per spike genome
        for spike in spike_species:
            spike_col = f"{spike}.normfactor.ipeff.adj"
            if spike_col not in seqstats.columns:
                print(f"Skipping {spike}, missing column {spike_col}")
                continue

            factor = row[spike_col]
            new_tagdir_name = f"{sample_id}.{spike}.normalized-tagdir"
            new_tagdir_path = output_dir / new_tagdir_name

            if new_tagdir_path.exists():
                print(f" Normalized tagdir exists: {new_tagdir_name}, skipping.")
                continue

            shutil.copytree(tagdir, new_tagdir_path)

            # Modify tagInfo.txt
            taginfo_path = new_tagdir_path / "tagInfo.txt"
            if not taginfo_path.exists():
                print(f" Missing tagInfo.txt in {new_tagdir_name}, skipping.")
                continue

            with open(taginfo_path, "r") as f:
                lines = f.readlines()

            try:
                cols = lines[1].strip().split("\t")
                original = float(cols[2])
                cols[2] = str(original * factor)
                lines[1] = "\t".join(cols) + "\n"
            except Exception as e:
                print(f"  Error modifying {taginfo_path}: {e}")
                continue

            with open(taginfo_path, "w") as f:
                f.writelines(lines)

            print(f" {sample_id} ({spike}) normalized: factor={factor}, {original} → {cols[2]}")


# ====================== CLI Interface ======================
def main():
    parser = argparse.ArgumentParser(
        description="Create HOMER tag directories for the target species and apply spike-in normalization."
    )
    parser.add_argument("--metadata_file", type=Path, help="Path to sample_metadata.norm.tsv")
    parser.add_argument("--user_dir", type=Path, default=Path.cwd())
    parser.add_argument("--target_species", required=True)
    parser.add_argument("--spike_species", nargs="+", required=True,
                        help="List of spike genomes for normalization (e.g. dm6 sac3)")

    args = parser.parse_args()

    user_dir = args.user_dir.resolve()
    target = args.target_species
    spikes = args.spike_species

    # Construct genome_dirs dictionary
    genome_dirs = {target: user_dir / f"{target}_data"}
    for spike in spikes:
        genome_dirs[spike] = user_dir / f"{spike}_data"

    metadata_file = args.metadata_file or user_dir / "sample_metadata.norm.tsv"

    normalize_tagdirs(metadata_file=metadata_file, genome_dirs=genome_dirs,
                      target_species=target, spike_species=spikes)


if __name__ == "__main__":
    main()
