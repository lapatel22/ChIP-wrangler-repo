#!/usr/bin/env python3
import os
import pandas as pd
from pathlib import Path
import numpy as np
import argparse
import sys
from contextlib import contextmanager
import warnings

# ============================================================
# Logging / Tee
# ============================================================
class Tee(object):
    def __init__(self, *streams):
        self.streams = streams
    
    def write(self, data):
        for s in self.streams:
            s.write(data)
    
    def flush(self):
        for s in self.streams:
            s.flush()

@contextmanager
def tee_stdout(log_path):
    with open(log_path, "w") as logfile:
        tee = Tee(sys.stdout, logfile)
        old_stdout = sys.stdout
        sys.stdout = tee
        try:
            yield
        finally:
            sys.stdout = old_stdout

# ============================================================
# QC REPORT FUNCTION
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

        if sid in meta_lookup.index:
            qc = meta_lookup.loc[sid, "QC_flag"]
            qc_display = qc if pd.notna(qc) and qc != "" else "within acceptable range"
            print(f"  Spike-in/target ratio: {qc_display}")
        else:
            print("  QC_flag: NOT FOUND")

        if sid in meta_lookup.index:
            dm6_nf = meta_lookup.loc[sid, f"{spike1_species}.normfactor.ipeff.adj"]
            sac3_nf = meta_lookup.loc[sid, f"{spike2_species}.normfactor.ipeff.adj"]
            dual_nf = meta_lookup.loc[sid, "dual.normfactor.ipeff.adj"]

            print(f"  {spike1_species} Normalization factor: {dm6_nf:.3f}")
            print(f"  {spike2_species} Normalization factor: {sac3_nf:.3f}")
            print(f"  dual Normalization factor: {dual_nf:.3f}")

            dm6_low, dm6_high = dm6_nf < norm_lo, dm6_nf > norm_hi
            sac3_low, sac3_high = sac3_nf < norm_lo, sac3_nf > norm_hi

            if (dm6_low and sac3_high) or (dm6_high and sac3_low):
                print("    **RED FLAG: dm6 and sac3 normalization disagree (opposite directions)**")

            if dm6_nf > 0 and sac3_nf > 0:
                diff = abs(dm6_nf - sac3_nf) / ((dm6_nf + sac3_nf) / 2)
                if diff > norm_diff_fraction:
                    print(f"    **WARNING: dm6 vs sac3 NF differ by >{norm_diff_fraction*100:.0f}%**")
                else:
                    print(f"  Spike-in normalization factors in similar range")
        else:
            print("  Normalization factors: NOT FOUND")

# ============================================================
# GC / Base Composition Functions
# ============================================================
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
            print(f"Warning: Missing GC files in {tagdir}, skipping.")
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
            sample_name = tagdir.replace(suffix_pattern, "").lstrip("0123456789_").replace("-", "")
            rows.append({
                "Sample": sample_name,
                "Sample_MeanGC": mean_sample,
                "Sample_MedianGC": median_sample,
                "Genome_MeanGC": mean_genome,
                "Genome_MedianGC": median_genome
            })
        except Exception as e:
            warnings.warn(f"Error reading GC data for {tagdir}: {e}")

    if not rows:
        return None
    return pd.DataFrame(rows)

# ============================================================
# WRAPPER-COMPATIBLE QC FUNCTION
# ============================================================

def generate_qc(metadata_file: Path, output_dir: Path,
                target_species: str, spike1_species: str, spike2_species: str):

    species_info = {
        target_species: {
            "tagdir_root": output_dir / f"{target_species}_data" / f"{target_species}_tagdirs",
            "suffix_pattern": f".{target_species}-tagdir"
        },
        spike1_species: {
            "tagdir_root": output_dir / f"{spike1_species}_data" / f"{spike1_species}_tagdirs",
            "suffix_pattern": f".{spike1_species}-tagdir"
        },
        spike2_species: {
            "tagdir_root": output_dir / f"{spike2_species}_data" / f"{spike2_species}_tagdirs",
            "suffix_pattern": f".{spike2_species}-tagdir"
        }
    }

    gc_stats_results = {}
    
    qc_output_dir = output_dir / "basedist_and_gc_outputs"
    qc_output_dir.mkdir(exist_ok=True, parents=True)
    
    for sp, info in species_info.items():
        gc_df = extract_gc_stats(info["tagdir_root"], info["suffix_pattern"])
        if gc_df is not None:
            gc_df["Genome"] = sp
            out_path = qc_output_dir / f"{sp}_gc_stats.tsv"   # <-- use subfolder
            gc_df.to_csv(out_path, sep="\t", index=False)
            gc_stats_results[sp] = out_path
            print(f"Saved {sp} GC stats → {out_path}")


    if gc_stats_results:
        all_gc = pd.concat([pd.read_csv(p, sep="\t") for p in gc_stats_results.values()], ignore_index=True)
        all_gc["GC_Ratio"] = all_gc["Sample_MeanGC"] / all_gc["Genome_MeanGC"]
        combined_path = qc_output_dir / "all_gc_stats_combined.tsv"
        all_gc.to_csv(combined_path, sep="\t", index=False)
        print(f"Combined GC stats → {combined_path}")

        if metadata_file.exists():
            generate_sample_qc_report(
                gc_stats_path=combined_path,
                sample_metadata_path=metadata_file,
                spike1_species=spike1_species,
                spike2_species=spike2_species
            )
        else:
            print("WARNING: sample_metadata.norm.tsv not found → QC report skipped.")

# ============================================================
# CLI ENTRY POINT
# ============================================================

def parse_args():
    parser = argparse.ArgumentParser(description="Compute GC bias, base composition, and QC report.")
    parser.add_argument("--target_species", required=True)
    parser.add_argument("--spike1_species", required=True)
    parser.add_argument("--spike2_species", required=True)
    parser.add_argument("--output_dir", type=Path, default=".")
    return parser.parse_args()

def main():
    args = parse_args()

    # Resolve base output directory
    base_output_dir = args.output_dir.resolve()
    metadata_file = base_output_dir / "sample_metadata.norm.tsv"

    # Create QC output folder
    output_dir = base_output_dir / "basedist_and_gc_outputs"
    output_dir.mkdir(exist_ok=True, parents=True)

    # Log file inside output_dir
    log_file = output_dir / "GC_bias_QC_report.log"

    # Tee stdout (both console + log)
    with tee_stdout(log_file):
        # Capture warnings as well
        with warnings.catch_warnings():
            warnings.simplefilter("always")
            print(f"Starting QC report. Log file: {log_file}\n")
            
            generate_qc(
                metadata_file=metadata_file,
                output_dir=output_dir,
                target_species=args.target_species,
                spike1_species=args.spike1_species,
                spike2_species=args.spike2_species
            )
            print("\nQC workflow complete!")



if __name__ == "__main__":
    main()
