#!/usr/bin/env python3
import os
import shutil
import subprocess
import pandas as pd

# ================== Required Inputs ====================

# metadata.tsv containing columns: library.ID, dm6.tagdir.normfactor.adj, sac3.tagdir.normfactor.adj
# target bam file
# suffix to remove

# === Base and output directories ===
base_dir = "../hg38_data/hg38_tagdirs"  # where tagdirs will be made
output_dir = "../hg38_data/hg38_normalized_tagdirs"
aligned_dir = "../hg38_data/hg38_aligned"  # where SAMs are located
target_genome = "hg38"
histogram_dir = "../hg38_data/hg38_histograms"

os.makedirs(base_dir, exist_ok=True)
os.makedirs(output_dir, exist_ok=True)

suffix_to_remove = ".hg38-tagdir"  # used later for normalization

print("Will output normalized tag directories in:", output_dir)
print("Will create tag directories from:", aligned_dir)

# === Step 1: Make tag directories from SAM files ===
for sam_file in os.listdir(aligned_dir):
    if not sam_file.endswith(".hg38.nosuffx2.sam"):
        continue

    sample_id = sam_file.replace(".hg38.nosuffx2.sam", "")
    tagdir_name = f"{sample_id}.hg38-tagdir"
    tagdir_path = os.path.join(base_dir, tagdir_name)
    sam_path = os.path.join(aligned_dir, sam_file)

    if os.path.exists(tagdir_path):
        print(f" Tag directory already exists for {sample_id}, skipping makeTagDirectory.")
        continue

    print(f" Creating tag directory for {sample_id}...")
    cmd = [
        "makeTagDirectory",
        tagdir_path,
        sam_path,
        "-genome", target_genome, 
        "-checkGC",
        "-fragLength", "150",
    ]

    try:
        subprocess.run(cmd, check=True)
        print(f" Created {tagdir_name}")
    except subprocess.CalledProcessError as e:
        print(f"Warning: Failed to create tag directory for {sample_id}: {e}")

# === Step 2: Load seqstats DataFrame ===
seqstats = pd.read_csv("../sample_metadata.norm.tsv", sep="\t")   

# === Clean column names ===
seqstats.columns = (
    seqstats.columns
    .str.replace(r"[\u200b\u00a0\r\n]", "", regex=True)
    .str.strip()
    .str.replace(" ", ".", regex=False)
)

print("Available columns in metadata:", seqstats.columns.tolist())

# === Verify that required columns exist ===
required_cols = ["library.ID", "dm6.normfactor.ipeff.adj", "sac3.normfactor.ipeff.adj"]
missing = [c for c in required_cols if c not in seqstats.columns]
if missing:
    raise ValueError(f"Error: Missing required columns in metadata: {missing}")

# === Step 3: Normalization loop ===
for tagdir in os.listdir(base_dir):
    full_tagdir_path = os.path.join(base_dir, tagdir)
    if not (os.path.isdir(full_tagdir_path) and tagdir.endswith(suffix_to_remove)):
        continue

    sample_id = tagdir.replace(suffix_to_remove, "")

    if sample_id not in seqstats["library.ID"].values:
        print(f"Warning: Sample ID {sample_id} not found in metadata, skipping.")
        continue

    # Extract normalization factors
    row = seqstats.loc[seqstats["library.ID"] == sample_id].iloc[0]
    fly_factor = row["dm6.normfactor.ipeff.adj"]
    yeast_factor = row["sac3.normfactor.ipeff.adj"]

    # Define normalization modes
    normalization_modes = [
        ("fly", fly_factor),
        ("yeast", yeast_factor),
    ]

    for mode, norm_factor in normalization_modes:
        new_tagdir_name = f"{sample_id}.{mode}.normalized-tagdir"
        new_tagdir_path = os.path.join(output_dir, new_tagdir_name)

        # Copy directory
        if os.path.exists(new_tagdir_path):
            print(f" Normalized tagdir already exists: {new_tagdir_name}, skipping.")
            continue

        shutil.copytree(full_tagdir_path, new_tagdir_path)

        # Modify tagInfo.txt
        taginfo_path = os.path.join(new_tagdir_path, "tagInfo.txt")
        if not os.path.exists(taginfo_path):
            print(f"Warning: tagInfo.txt missing in {new_tagdir_name}, skipping normalization.")
            continue

        with open(taginfo_path, "r") as f:
            lines = f.readlines()

        # Modify 2nd row (index 1), 3rd column (index 2)
        try:
            cols = lines[1].strip().split("\t")
            original_val = float(cols[2])
            cols[2] = str(original_val * norm_factor)
            lines[1] = "\t".join(cols) + "\n"
        except Exception as e:
            print(f" Error modifying {taginfo_path}: {e}")
            continue

        with open(taginfo_path, "w") as f:
            f.writelines(lines)

        print(f" {sample_id} ({mode}) normalized: factor={norm_factor}, {original_val} → {cols[2]}")
