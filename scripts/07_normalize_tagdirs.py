#!/usr/bin/env python3
import os
import shutil
import subprocess
import argparse
import pandas as pd
from pathlib import Path


# ====================== Parse Arguments ======================
def get_args():
    parser = argparse.ArgumentParser(
        description="Create HOMER tag directories and apply spike-in normalization."
    )

    parser.add_argument("--user_dir", required=True, help="Base project directory")
    parser.add_argument("--target_species", required=True, help="Target species (e.g. hg38)")
    parser.add_argument("--spike1_species", required=True, help="Spike-in species 1 (e.g. dm6)")
    parser.add_argument("--spike2_species", required=True, help="Spike-in species 2 (e.g. sac3)")

    return parser.parse_args()


# ====================== Main Script ======================
def main():

    args = get_args()

    user_dir = Path(args.user_dir).resolve()
    target = args.target_species
    spike1 = args.spike1_species
    spike2 = args.spike2_species

    # ================ Construct Directory Layout ================
    aligned_dir = user_dir / f"{target}_data" / f"{target}_aligned"
    base_dir = user_dir / f"{target}_data" / f"{target}_tagdirs"
    output_dir = user_dir / f"{target}_data" / f"{target}_normalized_tagdirs"

    os.makedirs(base_dir, exist_ok=True)
    os.makedirs(output_dir, exist_ok=True)

    suffix_to_remove = f".{target}-tagdir"

    print("Tagdirs will be written to:", base_dir)
    print("Normalized tagdirs will be written to:", output_dir)
    print("Looking for SAM files in:", aligned_dir)

    # ================= Step 1: Make tag directories =================
    for sam_file in aligned_dir.glob(f"*.{target}.nosuffx2.sam"):
        sample_id = sam_file.name.replace(f".{target}.nosuffx2.sam", "")
        tagdir_name = f"{sample_id}.{target}-tagdir"
        tagdir_path = base_dir / tagdir_name

        if tagdir_path.exists():
            print(f" Tagdir exists for {sample_id}, skipping makeTagDirectory.")
            continue

        print(f" Creating tag directory for {sample_id}...")

        cmd = [
            "makeTagDirectory",
            str(tagdir_path),
            str(sam_file),
            "-genome", target,
            "-checkGC",
            "-fragLength", "150",
        ]

        try:
            subprocess.run(cmd, check=True)
            print(" Created:", tagdir_name)
        except subprocess.CalledProcessError as e:
            print(f" Failed for {sample_id}: {e}")

    # ================= Step 2: Load metadata =================
    seqstats_file = user_dir / "sample_metadata.norm.tsv"
    seqstats = pd.read_csv(seqstats_file, sep="\t")

    seqstats.columns = (
        seqstats.columns
        .str.replace(r"[\u200b\u00a0\r\n]", "", regex=True)
        .str.strip()
        .str.replace(" ", ".", regex=False)
    )

    print("\nAvailable columns:", seqstats.columns.tolist())

    # Required columns
    spike1_col = f"{spike1}.normfactor.ipeff.adj"
    spike2_col = f"{spike2}.normfactor.ipeff.adj"

    required = ["library.ID", spike1_col, spike2_col]
    missing = [c for c in required if c not in seqstats.columns]
    if missing:
        raise ValueError(f"ERROR: Missing required columns: {missing}")

    # ================= Step 3: Normalization =================
    for tagdir in base_dir.iterdir():
        if not (tagdir.is_dir() and tagdir.name.endswith(suffix_to_remove)):
            continue

        sample_id = tagdir.name.replace(suffix_to_remove, "")

        if sample_id not in seqstats["library.ID"].values:
            print(f" Sample {sample_id} not found in metadata, skipping.")
            continue

        row = seqstats.loc[seqstats["library.ID"] == sample_id].iloc[0]

        spike1_factor = row[spike1_col]
        spike2_factor = row[spike2_col]

        normalization_modes = [
            (spike1, spike1_factor),
            (spike2, spike2_factor),
        ]

        for species, factor in normalization_modes:
            new_tagdir_name = f"{sample_id}.{species}.normalized-tagdir"
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

            print(
                f" {sample_id} ({species}) normalized: "
                f"factor={factor}, {original} → {cols[2]}"
            )


if __name__ == "__main__":
    main()
