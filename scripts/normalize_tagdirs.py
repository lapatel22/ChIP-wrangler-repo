#!/usr/bin/env python3
import os
import shutil
import subprocess
import argparse
import pandas as pd
from pathlib import Path

def normalize_tagdirs(
    metadata_file: Path,
    output_dir: Path,
    target_species: str,
    spike1_species: str = None,
    spike2_species: str = None,
    force_overwrite: bool = False
):
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
    target_dir = output_dir / f"{target_species}_data"
    aligned_dir = target_dir / f"{target_species}_aligned"
    base_dir = target_dir / f"{target_species}_tagdirs"
    output_norm_dir = target_dir / f"{target_species}_normalized_tagdirs"

    spike_species = [s for s in [spike1_species, spike2_species] if s is not None]

    os.makedirs(base_dir, exist_ok=True)
    os.makedirs(output_norm_dir, exist_ok=True)

    suffix_to_remove = f".{target_species}-tagdir"

    print(f"\nProcessing target species: {target_species}")
    print("Tagdirs will be written to:", base_dir)
    print("Normalized tagdirs will be written to:", output_norm_dir)
    print("Looking for SAM files in:", aligned_dir)

    # ================= Step 1: Make tag directories =================
    for sam_file in aligned_dir.glob(f"*.{target_species}.nosuffx2.sam"):
        sample_id = sam_file.name.replace(f".{target_species}.nosuffx2.sam", "")
        tagdir_name = f"{sample_id}.{target_species}-tagdir"
        tagdir_path = base_dir / tagdir_name

        if tagdir_path.exists() and not force_overwrite:
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
            new_tagdir_path = output_norm_dir / new_tagdir_name

    # NOTE: previously here i had an if statement, and skipped normalizing a tagdir if it already existed. given normalization is fast, we now remake every time

            # remove old one
            if new_tagdir_path.exists():
                print(f" Overwriting existing normalized tagdir: {new_tagdir_name}")
                shutil.rmtree(new_tagdir_path)

            # copy original tagdir
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
    parser.add_argument("--output_dir", type=Path, default=Path.cwd())
    parser.add_argument("force_overwrite", action="store_true", 
                       help = "Global option to force overwriting of every step")

    args = parser.parse_args()

    output_dir = args.output_dir.resolve()
    force_overwrite=args.force_overwrite
    metadata_file = args.metadata_file or output_dir / "sample_metadata.norm.tsv"

    # Construct genome_dirs dictionary
    genome_dirs = {target: output_dir / f"{target_species}_data"}
    genome_dirs[spike] = output_dir / f"{spike}_data"

    normalize_tagdirs(metadata_file=metadata_file,
                      target_species=target, spike_species=spikes,
                     force_overwrite=force_overwrite)


if __name__ == "__main__":
    main()


### NOTES TO ADD: add the commands run as an output/log file
### add a line in the tagInfo file that records the normalization factor used (and a record that ChIP-wrangler modified it!)
### add a note that the genoems must be configured with HOMER!! esp spike-ins
