#!/usr/bin/env python3
import os
import pandas as pd
from pathlib import Path
import numpy as np
import argparse

# ============================================================
# ARGUMENTS
# ============================================================
def parse_args():
    parser = argparse.ArgumentParser(description="Compute GC bias, base composition, and QC report.")
    parser.add_argument("--user_dir", required=True, help="Base directory where *_data/ directories exist.")
    parser.add_argument("--target_species", required=True, help="Primary genome (e.g. hg38)")
    parser.add_argument("--spike1_species", required=True, help="Spike-in species 1 (e.g. dm6)")
    parser.add_argument("--spike2_species", required=True, help="Spike-in species 2 (e.g. sac3)")
    return parser.parse_args()


# ============================================================
# QC REPORT FUNCTION (added)
# ============================================================
def generate_sample_qc_report(
    gc_stats_path: str,
    sample_metadata_path: str,
    spike1_species: str,
    spike2_species: str,
    gc_hi: float = 1.2,
    gc_lo: float = 0.8,
    norm_hi: float = 1.1,
    norm_lo: float = 0.9,
    norm_diff_fraction: float = 0.20
):
    gc = pd.read_csv(gc_stats_path, sep="\t")
    meta = pd.read_csv(sample_metadata_path, sep="\t")

    meta_lookup = meta.set_index("library.ID")
    sample_ids = gc["Sample"].unique()

    print("\n\n==========================")
    print(" QC REPORT")
    print("==========================")

    for sid in sample_ids:
        print(f"\n----------------------------------")
        print(f"Sample: {sid}")
        print("----------------------------------")

        sub = gc[gc["Sample"] == sid]
        is_input = ("input" in sid.lower())

        # ---- GC RATIOS ----
        for _, row in sub.iterrows():
            genome = row["Genome"]
            ratio = row["GC_Ratio"]
            msg = f"  GC ratio ({genome}): {ratio:.3f}"

            if is_input:
                if ratio > gc_hi:
                    msg += "  **RED FLAG: GC high**"
                elif ratio < gc_lo:
                    msg += "  **RED FLAG: GC low**"
                else:
                    msg += " ---> GC content acceptable"

            print(msg)

        # ---- QC FLAG ----
        if sid in meta_lookup.index:
            qc = meta_lookup.loc[sid, "QC_flag"]
            qc_display = qc if pd.notna(qc) and qc != "" else "within acceptable range"
            print(f"  Spike-in/target ratio: {qc_display}")
        else:
            print("  QC_flag: NOT FOUND")

        # ---- Normalization factors ----
        if sid in meta_lookup.index:
            dm6_nf = meta_lookup.loc[sid, "dm6.normfactor.ipeff.adj"]
            sac3_nf = meta_lookup.loc[sid, "sac3.normfactor.ipeff.adj"]
            dual_nf = meta_lookup.loc[sid, "dual.normfactor.ipeff.adj"]

            print(f"  {spike1_species} Normalization factor:  {dm6_nf:.3f}")
            print(f"  {spike2_species} Normalization factor:  {sac3_nf:.3f}")
            print(f"  dual Normalization factor: {dual_nf:.3f}")

            # Out-of-range red flags
            dm6_low  = dm6_nf < norm_lo
            dm6_high = dm6_nf > norm_hi
            sac3_low  = sac3_nf < norm_lo
            sac3_high = sac3_nf > norm_hi

            # Trigger red flag only if they disagree in opposite directions
            if (dm6_low and sac3_high) or (dm6_high and sac3_low):
                print("    **RED FLAG: dm6 and sac3 normalization disagree (opposite directions)**")

            # Difference between dm6 and sac3
            if dm6_nf > 0 and sac3_nf > 0:
                diff = abs(dm6_nf - sac3_nf) / ((dm6_nf + sac3_nf) / 2)
                if diff > norm_diff_fraction:
                    print(f"    **WARNING: dm6 vs sac3 NF differ by >{norm_diff_fraction*100:.0f}%**")
                elif diff < norm_diff_fraction:
                    print(f"  Spike-in normalization factors in similar range")
        else:
            print("  Normalization factors: NOT FOUND")



# ============================================================
# BASE COMPOSITION + GC CONTENT FUNCTIONS
# (your existing code unchanged)
# ============================================================
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

        df = pd.read_csv(freq_path, sep="\t", comment="#")
        if "Offset" not in df.columns:
            df.rename(columns={df.columns[0]: "Offset"}, inplace=True)
        df_long = df.melt(id_vars=["Offset"], value_vars=["A", "C", "G", "T"],
                          var_name="Base", value_name="Frequency")

        sample = tagdir.replace(suffix_pattern, "")
        sample = sample.lstrip("0123456789_").replace("-", "")
        df_long["Sample"] = sample
        all_dfs.append(df_long)

    if not all_dfs:
        print(f"No usable tagFreqUniq.txt found in {tagdir_root}. Make sure to run makeTagDirectories with -checkGC option.")
        return None

    combined = pd.concat(all_dfs, ignore_index=True)
    return combined


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
            print(f"Error reading GC data for {tagdir}: {e}")

    if not rows:
        print(f"No valid GC files found in {tagdir_root}. Make sure to run makeTagDirectories with -checkGC option.")
        return None

    return pd.DataFrame(rows)



# ============================================================
# MAIN
# ============================================================
def main():
    args = parse_args()
    user_dir = Path(args.user_dir)

    # species_info now dynamically constructed
    species_info = {
        args.target_species: {
            "tagdir_root": user_dir / f"{args.target_species}_data" / f"{args.target_species}_tagdirs",
            "suffix_pattern": f".{args.target_species}-tagdir"
        },
        args.spike1_species: {
            "tagdir_root": user_dir / f"{args.spike1_species}_data" / f"{args.spike1_species}_tagdirs",
            "suffix_pattern": f".{args.spike1_species}-tagdir"
        },
        args.spike2_species: {
            "tagdir_root": user_dir / f"{args.spike2_species}_data" / f"{args.spike2_species}_tagdirs",
            "suffix_pattern": f".{args.spike2_species}-tagdir"
        }
    }

    output_dir = Path("basedist_and_gc_outputs")
    output_dir.mkdir(exist_ok=True)

    gc_stats_results = {}

    # ---- MAIN LOOP ----
    for sp, info in species_info.items():
        root = info["tagdir_root"]
        suffix = info["suffix_pattern"]

        # GC stats
        gc_df = extract_gc_stats(root, suffix)
        if gc_df is not None:
            gc_df["Genome"] = sp
            out_path = output_dir / f"{sp}_gc_stats.tsv"
            gc_df.to_csv(out_path, sep="\t", index=False)
            gc_stats_results[sp] = out_path
            print(f"Success: Saved {sp} GC stats → {out_path}")

    # ---- COMBINE ----
    if gc_stats_results:
        all_gc = pd.concat([pd.read_csv(p, sep="\t") for p in gc_stats_results.values()], ignore_index=True)
        all_gc["GC_Ratio"] = all_gc["Sample_MeanGC"] / all_gc["Genome_MeanGC"]
        combined_path = output_dir / "all_gc_stats_combined.tsv"
        all_gc.to_csv(combined_path, sep="\t", index=False)
        print(f"Success: Combined GC stats → {combined_path}")

        # ---- NEW QC REPORT SECTION ----
        metadata_norm_path = user_dir / "sample_metadata.norm.tsv"
        if metadata_norm_path.exists():
            generate_sample_qc_report(
                gc_stats_path=combined_path,
                sample_metadata_path=metadata_norm_path,
                spike1_species=args.spike1_species,
                spike2_species=args.spike2_species
            )
        else:
            print("WARNING: sample_metadata.norm.tsv not found → QC report skipped.")


if __name__ == "__main__":
    main()
