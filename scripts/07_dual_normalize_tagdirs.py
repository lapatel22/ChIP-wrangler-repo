#!/usr/bin/env python3
import os
import shutil
import subprocess
import argparse
import pandas as pd
from pathlib import Path


# ======================== Parse Arguments ========================
def get_args():
    parser = argparse.ArgumentParser(description="Generate and normalize HOMER tag directories.")
    parser.add_argument("--user_dir", required=True, help="Base project directory")
    parser.add_argument("--target_species", required=True, help="Target genome species, e.g. hg38")
    parser.add_argument("--spike1_species", required=True, help="Spike-in species 1")
    parser.add_argument("--spike2_species", required=True, help="Spike-in species 2")

    return parser.parse_args()


# ======================== Main Script ========================
def main():
    args = get_args()

    user_dir = Path(args.user_dir).resolve()
    target_species = args.target_species

    # ======================== Config Paths ========================
    aligned_dir = user_dir / f"{target_species}_data" / f"{target_species}_aligned"
    base_tagdir_dir = user_dir / f"{target_species}_data" / f"{target_species}_tagdirs"
    normalized_tagdir_dir = user_dir / f"{target_species}_data" / f"{target_species}_normalized_tagdirs"

    # Suffixes
    sam_suffix = ".nosuffx2.sam"
    tagdir_suffix = f".{target_species}-tagdir"
    normalized_suffix = ".normalized-tagdir"

    # Seqstats file
    seqstats_file = user_dir / "sample_metadata.norm.tsv"

    # Column with normalization factor
    norm_col = "dual.normfactor.ipeff.adj"

    # Create output dir
    normalized_tagdir_dir.mkdir(parents=True, exist_ok=True)

    # ======================== Load seqstats ========================
    seqstats = pd.read_csv(seqstats_file, sep="\t")

    seqstats.columns = (
        seqstats.columns
        .str.replace(r"[\u200b\u00a0\r\n]", "", regex=True)
        .str.strip()
        .str.replace(" ", ".", regex=False)
    )

    print("Available columns in seqstats:", seqstats.columns.tolist())

    # ======================== Step 1: Create tagdirs ========================
    if base_tagdir_dir.exists() and any(base_tagdir_dir.iterdir()):
        print(f"{base_tagdir_dir} exists and is not empty — skipping tagdir creation.")
    else:
        base_tagdir_dir.mkdir(parents=True, exist_ok=True)
        for sam_file in aligned_dir.glob(f"*{sam_suffix}"):

            sample_id = sam_file.name.replace(f".{target_species}.nosuffx2.sam", "")
            tagdir_name = f"{sample_id}{tagdir_suffix}"
            tagdir_path = base_tagdir_dir / tagdir_name

            if tagdir_path.exists():
                print(f"Tagdir already exists for {sample_id}, skipping.")
                continue

            print(f" Creating tagdir for {sample_id}...")

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
                print(f"Created {tagdir_name}")
            except subprocess.CalledProcessError as e:
                print(f" Failed to create tagdir for {sample_id}: {e}")

    # ======================== Step 2: Normalize tagdirs ========================
    for tagdir in base_tagdir_dir.glob(f"*{tagdir_suffix}"):
        sample_id = tagdir.name.replace(tagdir_suffix, "")

        if sample_id not in seqstats["library.ID"].values:
            print(f" Sample {sample_id} missing in seqstats metadata — skipping.")
            continue

        norm_factor = seqstats.loc[
            seqstats["library.ID"] == sample_id, norm_col
        ].values[0]

        if pd.isna(norm_factor):
            print(f" Norm factor for {sample_id} is NaN — skipping.")
            continue

        normalized_tagdir_name = f"{sample_id}{normalized_suffix}"
        normalized_tagdir_path = normalized_tagdir_dir / normalized_tagdir_name

        if normalized_tagdir_path.exists():
            print(f"✓ Normalized tagdir already exists for {sample_id}, skipping.")
            continue

        shutil.copytree(tagdir, normalized_tagdir_path)

        # Modify tagInfo.txt
        taginfo_path = normalized_tagdir_path / "tagInfo.txt"
        if taginfo_path.exists():
            with open(taginfo_path, "r") as f:
                lines = f.readlines()

            cols = lines[1].strip().split("\t")
            original_val = float(cols[2])
            new_val = original_val * norm_factor
            cols[2] = str(new_val)
            lines[1] = "\t".join(cols) + "\n"

            with open(taginfo_path, "w") as f:
                f.writelines(lines)

            print(
                f"🔧 Modified {sample_id}: factor={norm_factor}, "
                f"original={original_val} → new={new_val}"
            )
        else:
            print(f" tagInfo.txt missing in {normalized_tagdir_path}, skipping.")


if __name__ == "__main__":
    main()
