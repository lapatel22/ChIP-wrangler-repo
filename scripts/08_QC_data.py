#!/usr/bin/env python3
import os
import pandas as pd
from pathlib import Path
import numpy as np

# ===================== CONFIG ===================== #
species_info = {
    "hg38": {"tagdir_root": "../hg38_data/hg38_tagdirs", "suffix_pattern": ".hg38-tagdir"},
    "dm6":  {"tagdir_root": "../dm6_data/dm6_tagdirs", "suffix_pattern": ".dm6-tagdir"},
    "sac3": {"tagdir_root": "../sac3_data/sac3_tagdirs", "suffix_pattern": ".sac3-tagdir"}
}

output_dir = Path("basedist_and_gc_outputs")
output_dir.mkdir(exist_ok=True)

# ===================== PART 1: BASE COMPOSITION ===================== #
def generate_basedist(tagdir_root, suffix_pattern):
    tagdirs = [d for d in os.listdir(tagdir_root) if d.endswith("-tagdir")]
    if not tagdirs:
        print(f"Error: No tagdirs found in {tagdir_root}")
        return None

    print(f"Processing base composition for {tagdir_root} ({len(tagdirs)} tagdirs)")

    all_dfs = []
    for tagdir in tagdirs:
        freq_path = Path(tagdir_root) / tagdir / "tagFreqUniq.txt"
        if not freq_path.exists():
            print(f"Warning: Missing {freq_path}, skipping. Make sure to run makeTagDirectories with -checkGC option.")
            continue

        # Read tag frequency table
        df = pd.read_csv(freq_path, sep="\t", comment="#")
        if "Offset" not in df.columns:
            df.rename(columns={df.columns[0]: "Offset"}, inplace=True)
        df_long = df.melt(id_vars=["Offset"], value_vars=["A", "C", "G", "T"],
                          var_name="Base", value_name="Frequency")

        # Clean sample name
        sample = tagdir.replace(suffix_pattern, "")
        sample = sample.lstrip("0123456789_").replace("-", "")
        df_long["Sample"] = sample
        all_dfs.append(df_long)

    if not all_dfs:
        print(f"No usable tagFreqUniq.txt found in {tagdir_root}. Make sure to run makeTagDirectories with -checkGC option.")
        return None

    combined = pd.concat(all_dfs, ignore_index=True)
    return combined


# ===================== PART 2: GC CONTENT STATS ===================== #
def weighted_median(gc, pdf):
    df = pd.DataFrame({"gc": gc, "pdf": pdf}).sort_values("gc")
    df["cumPDF"] = df["pdf"].cumsum() / df["pdf"].sum()
    return np.interp(0.5, df["cumPDF"], df["gc"])

def extract_gc_stats(tagdir_root, suffix_pattern):
    tagdirs = [d for d in os.listdir(tagdir_root) if d.endswith("-tagdir")]
    if not tagdirs:
        print(f"No tagdirs found in {tagdir_root}")
        return None

    rows = []
    for tagdir in tagdirs:
        folder = Path(tagdir_root) / tagdir
        sample_gc_path = folder / "tagGCcontent.txt"
        genome_gc_path = folder / "genomeGCcontent.txt"

        if not sample_gc_path.exists() or not genome_gc_path.exists():
            print(f"Warning: Missing GC files in {tagdir}, skipping. Make sure to run makeTagDirectories with -checkGC option.")
            continue

        try:
            sample_gc = pd.read_csv(sample_gc_path, sep="\t")
            genome_gc = pd.read_csv(genome_gc_path, sep="\t")

            # Convert to numeric and clean
            for df in [sample_gc, genome_gc]:
                df["GC%"] = pd.to_numeric(df["GC%"], errors="coerce")
                df["Normalized Fraction(PDF)"] = pd.to_numeric(df["Normalized Fraction(PDF)"], errors="coerce")
                df.dropna(subset=["GC%", "Normalized Fraction(PDF)"], inplace=True)

            mean_sample = (sample_gc["GC%"] * sample_gc["Normalized Fraction(PDF)"]).sum() / sample_gc["Normalized Fraction(PDF)"].sum()
            mean_genome = (genome_gc["GC%"] * genome_gc["Normalized Fraction(PDF)"]).sum() / genome_gc["Normalized Fraction(PDF)"].sum()

            median_sample = weighted_median(sample_gc["GC%"], sample_gc["Normalized Fraction(PDF)"])
            median_genome = weighted_median(genome_gc["GC%"], genome_gc["Normalized Fraction(PDF)"])

            sample_name = tagdir.replace(suffix_pattern, "")
            sample_name = sample_name.lstrip("0123456789_").replace("-", "")

            rows.append({
                "Sample": sample_name,
                "Sample_MeanGC": mean_sample,
                "Sample_MedianGC": median_sample,
                "Genome_MeanGC": mean_genome,
                "Genome_MedianGC": median_genome
            })
        except Exception as e:
            print(f"Error reading GC data for {tagdir}: {e} ")

    if not rows:
        print(f"No valid GC files found in {tagdir_root}. Make sure to run makeTagDirectories with -checkGC option.")
        return None

    return pd.DataFrame(rows)


# ===================== MAIN LOOP ===================== #
basedist_results = {}
gc_stats_results = {}

for sp, info in species_info.items():
    root = info["tagdir_root"]
    suffix = info["suffix_pattern"]

    # --- Base composition ---
    bd = generate_basedist(root, suffix)
    if bd is not None:
        out_path = output_dir / f"{sp}_basedist.tsv"
        bd.to_csv(out_path, sep="\t", index=False)
        basedist_results[sp] = out_path
        print(f"Success: Saved {sp} base distribution → {out_path}")

    # --- GC stats ---
    gc_df = extract_gc_stats(root, suffix)
    if gc_df is not None:
        gc_df["Genome"] = sp
        out_path = output_dir / f"{sp}_gc_stats.tsv"
        gc_df.to_csv(out_path, sep="\t", index=False)
        gc_stats_results[sp] = out_path
        print(f"Success: Saved {sp} GC stats → {out_path}")

# Combine GC stats for all species into one file
if gc_stats_results:
    all_gc = pd.concat([pd.read_csv(p, sep="\t") for p in gc_stats_results.values()], ignore_index=True)
    all_gc["GC_Ratio"] = all_gc["Sample_MeanGC"] / all_gc["Genome_MeanGC"]
    combined_path = output_dir / "all_gc_stats_combined.tsv"
    all_gc.to_csv(combined_path, sep="\t", index=False)
    print(f"Success: Combined GC stats → {combined_path}")

print("Success: Completed all species base composition and GC bias extraction.")
