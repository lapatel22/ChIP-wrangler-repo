#!/usr/bin/env python3

import argparse
import subprocess
from pathlib import Path
import pandas as pd
from typing import Optional


# ============================================================
# Core function (callable from wrapper)
# ============================================================

def create_counts_matrix(
    target: str,
    metadata_tsv: Path,
    conditions: str,
    output_dir: Path,
    annot_size: str = "given",
    out_tsv: Optional[Path] = None,
) -> Path:
    
    """
    Create HOMER counts matrix for DESeq2 using mergePeaks + annotatePeaks.

    Returns
    -------
    Path
        Path to output counts matrix TSV
    """

    output_dir = output_dir.resolve()
    conds = conditions.split(",")

    if len(conds) != 2:
        raise ValueError("--conditions must contain exactly two values")

    condA, condB = conds

    # ------------------------
    # Load metadata
    # ------------------------
    meta = pd.read_csv(metadata_tsv, sep="\t")

    required_cols = {"library.ID", "Condition"}
    missing = required_cols - set(meta.columns)
    if missing:
        raise ValueError(f"Metadata missing required columns: {missing}")

    # ------------------------
    # Paths
    # ------------------------
    data_dir = output_dir / f"{target}_data"
    tagdir_root = data_dir / f"{target}_tagdirs"
    peaks_dir = data_dir

    merged_peaks = data_dir / f"merged_peaks_{condA}_vs_{condB}.txt"

    # AUTO-GENERATE OUTPUT NAME IF NOT PROVIDED
    if out_tsv is None:
        out_tsv = data_dir / f"raw_counts_from_merged_{condA}_{condB}.tsv"

    # ------------------------
    # Select samples for conditions
    # ------------------------
    meta_sub = meta[meta["Condition"].isin([condA, condB])]

    if meta_sub.empty:
        raise RuntimeError("No samples found for specified conditions")

    # ------------------------
    # Collect peak files
    # ------------------------
    peak_files = []
    for lib_id in meta_sub["library.ID"].unique():
        peak_file = peaks_dir / f"{lib_id}.regions.txt"
        if peak_file.exists():
            peak_files.append(str(peak_file))

    if not peak_files:
        raise RuntimeError("No peak files found for specified conditions")

    # ------------------------
    # mergePeaks
    # ------------------------
    merge_cmd = ["mergePeaks"] + peak_files

    print("Running mergePeaks:")
    print("CMD:", " ".join(merge_cmd))

    with open(merged_peaks, "w") as fh:
        subprocess.run(merge_cmd, stdout=fh, check=True)

    # ------------------------
    # Collect tag directories
    # ------------------------
    tagdirs = []
    sample_cols = []
    for lib_id, ip_type in zip(meta_sub["library.ID"], meta_sub["IP"]):
        if ip_type == "input":
            continue  # skip input
        tagdir = tagdir_root / f"{lib_id}.{target}-tagdir"
        if not tagdir.exists():
            raise FileNotFoundError(f"Missing tagdir: {tagdir}")
        tagdirs.append(str(tagdir))
        sample_cols.append(lib_id)
    
    # Use only IP samples for annotatePeaks output
    annotate_cmd = [
        "annotatePeaks.pl",
        str(merged_peaks),
        target,
        "-size", str(annot_size),
        "-raw",
        "-d",
    ] + tagdirs

    print("Running annotatePeaks:")
    print("CMD:", " ".join(annotate_cmd))

    with open(out_tsv, "w") as fh:
        subprocess.run(annotate_cmd, stdout=fh, check=True)

    print(f"\n✔ Counts matrix written to {out_tsv}\n")

    return out_tsv


# ============================================================
# CLI wrapper (thin)
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description="Create HOMER counts matrix for DESeq2"
    )

    parser.add_argument("--target", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--conditions", required=True)

    parser.add_argument("--out_tsv", default=None, 
                        help="Optional output counts file (auto-named if not provided)")

    parser.add_argument("--annot_size", default="given",
        help="Peak size used for annotatePeaks")

    parser.add_argument("--output_dir", type=Path, default=Path.cwd(),
        help="Path to the ChIP-wrangler project directory")

    args = parser.parse_args()

    create_counts_matrix(
        target=args.target,
        metadata_tsv=Path(args.metadata),
        conditions=args.conditions,
        output_dir=args.output_dir,
        annot_size=args.annot_size,
        out_tsv=Path(args.out_tsv) if args.out_tsv else None,
    )


if __name__ == "__main__":
    main()
