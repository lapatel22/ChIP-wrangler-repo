#!/usr/bin/env python3
import os
import shutil
import subprocess
import pandas as pd
from pathlib import Path

# ======================== Config ========================
# Paths
aligned_dir = Path("../hg38_data/hg38_aligned")        # SAM files
base_tagdir_dir = Path("../hg38_data/hg38_tagdirs")   # Original tagdirs
normalized_tagdir_dir = Path("../hg38_data/hg38_normalized_tagdirs")  # Output normalized tagdirs

# Genome
target_genome = "hg38"

# Suffixes
sam_suffix = ".nosuffx2.sam"
tagdir_suffix = ".hg38-tagdir"
normalized_suffix = ".normalized-tagdir"

# Seqstats file
seqstats_file = Path("../sample_metadata.norm.tsv")

# Column with normalization factor
norm_col = "dual.normfactor.ipeff.adj"

# Create output directory if it does not exist
normalized_tagdir_dir.mkdir(parents=True, exist_ok=True)

# ======================== Load seqstats ========================
seqstats = pd.read_csv(seqstats_file, sep="\t")

# Clean column names of hidden characters and whitespace
seqstats.columns = (
    seqstats.columns
    .str.replace(r"[\u200b\u00a0\r\n]", "", regex=True)
    .str.strip()
    .str.replace(" ", ".", regex=False)
)

print("Available columns in seqstats:", seqstats.columns.tolist())

# ======================== Step 1: Create tag directories from SAM files ========================
if base_tagdir_dir.exists() and any(base_tagdir_dir.iterdir()):
    print(f"📂 {base_tagdir_dir} already exists and is not empty. Skipping tag directory creation.")
else:
    base_tagdir_dir.mkdir(parents=True, exist_ok=True)
    for sam_file in aligned_dir.glob(f"*{sam_suffix}"):
        # Remove both .hg38 and .nosuffx2 to get the sample ID
        sample_id = sam_file.name.replace(".hg38.nosuffx2.sam", "")
        tagdir_name = f"{sample_id}{tagdir_suffix}"
        tagdir_path = base_tagdir_dir / tagdir_name

        if tagdir_path.exists():
            print(f"Tag directory already exists for {sample_id}, skipping makeTagDirectory.")
            continue

        print(f"🧬 Creating tag directory for {sample_id}...")
        cmd = [
            "makeTagDirectory",
            str(tagdir_path),
            str(sam_file),
            "-genome", target_genome,
            "-checkGC",
            "-fragLength", "150",
        ]
        try:
            subprocess.run(cmd, check=True)
            print(f"Created {tagdir_name}")
        except subprocess.CalledProcessError as e:
            print(f"Failed to create tag directory for {sample_id}: {e}")

# ======================== Step 2: Normalize tag directories ========================
for tagdir in base_tagdir_dir.glob(f"*{tagdir_suffix}"):
    # Remove .hg38-tagdir to match seqstats library.ID
    sample_id = tagdir.name.replace(tagdir_suffix, "")

    if sample_id not in seqstats["library.ID"].values:
        print(f"Warning: Sample ID {sample_id} not found in seqstats, skipping.")
        continue

    norm_factor = seqstats.loc[seqstats["library.ID"] == sample_id, norm_col].values[0]

    if pd.isna(norm_factor):
        print(f"Warning: Normalization factor for {sample_id} is NaN, skipping.")
        continue

    # Copy to normalized output directory
    normalized_tagdir_name = f"{sample_id}{normalized_suffix}"
    normalized_tagdir_path = normalized_tagdir_dir / normalized_tagdir_name
    if normalized_tagdir_path.exists():
        print(f"Success: Normalized tag directory already exists for {sample_id}, skipping copy.")
        continue

    shutil.copytree(tagdir, normalized_tagdir_path)

    # Modify tagInfo.txt
    taginfo_path = normalized_tagdir_path / "tagInfo.txt"
    if taginfo_path.exists():
        with open(taginfo_path, "r") as f:
            lines = f.readlines()

        # 2nd row (index 1), 3rd column (index 2)
        cols = lines[1].strip().split("\t")
        original_val = float(cols[2])
        cols[2] = str(original_val * norm_factor)
        lines[1] = "\t".join(cols) + "\n"

        with open(taginfo_path, "w") as f:
            f.writelines(lines)

        print(f"Modified: {sample_id}, factor {norm_factor}, original {original_val} → new {cols[2]}")
    else:
        print(f"Warning: tagInfo.txt not found in {normalized_tagdir_path}, skipping modification.")
