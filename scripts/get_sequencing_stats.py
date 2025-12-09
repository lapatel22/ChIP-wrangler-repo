#!/usr/bin/env python3
"""
05_get_sequencing_stats.py
Update sample_metadata.tsv with sequencing stats, spike-in ratios, QC flags,
and normalized values based on per-species SAM files ending with .nosuffx2.sam.
"""

import pandas as pd
import subprocess
import re
import argparse
from pathlib import Path


def update_sample_metadata(
    output_dir: str,
    target_genome: str,
    spike1_genome: str,
    spike2_genome: str,
    samtools_path: str = "samtools",
    metadata: str = "./sample_names.tsv"
):
    """
    Update sample_metadata.tsv with sequencing stats, spike-in ratios, QC flags,
    and normalization values.

    All required SAMs must follow naming:

        <sample>.<species>.nosuffx2.sam

    The <sample> portion must match:
        Cell_sync_inter_biorep_IP_techrep
    """
    output_dir = Path(output_dir)

    # ---------------- Required input file ----------------
    metadata_file = output_dir / "sample_names.tsv"
    if not metadata_file.exists():
        raise FileNotFoundError(f"Error: sample_names.tsv not found at {metadata_file}")

    df_meta = pd.read_csv(metadata_file, sep="\t")
    if df_meta.empty:
        raise RuntimeError(f"ERROR: sample_names.tsv at {metadata_file} is EMPTY.")

    print(f"Reading metadata from: {metadata_file}")

    # ---------------- Species folders ----------------
    target_dir = output_dir / f"{target_genome}_data" / f"{target_genome}_aligned"
    spike1_dir = output_dir / f"{spike1_genome}_data" / f"{spike1_genome}_aligned"
    spike2_dir = output_dir / f"{spike2_genome}_data" / f"{spike2_genome}_aligned"


    # ---------------- Gather SAMs ----------------
    sam_files = list(target_dir.glob("*.nosuffx2.sam"))
    if not sam_files:
        raise FileNotFoundError(f"Warning: No .nosuffx2.sam files found in {target_dir}")

    records = []

    # ---------------- Process each sample ----------------
    for sam in sam_files:
        sample_full = sam.stem
        # remove trailing .<genome>.nosuffx2 from stem
        sample = re.sub(rf"\.{target_genome}\.nosuffx2$", "", sample_full)
        print(sample)
        # Expected naming pattern
        match = re.match(
            r"^([A-Za-z0-9]+_[A-Za-z0-9]+_[A-Za-z0-9]+)_(\d+)_(\w+)_([0-9]+)$",
            sample
        )

        if not match:
            print(f"Warning: Skipping unrecognized sample name: {sample_full}")
            continue

        condition_base, biorep, ip, techrep = match.groups()

        # Build corresponding spike-in SAM paths
      #  print(sam)
        sam_target = target_dir / f"{sample}.{target_genome}.bam"
        sam_spike1 = spike1_dir / f"{sample}.{spike1_genome}.bam"
        sam_spike2 = spike2_dir / f"{sample}.{spike2_genome}.bam"

        if not (sam_target.exists() and sam_spike1.exists() and sam_spike2.exists()):
            print(f" Missing SAM(s) for sample {sample} — skipping.")
            continue

        # Count reads
        def count_reads(path):
            result = subprocess.run(
                [samtools_path, "view", "-c", str(path)],
                capture_output=True, text=True
            )
            return int(result.stdout.strip() or 0)

        count_target = count_reads(sam_target)
        count_s1 = count_reads(sam_spike1)
        count_s2 = count_reads(sam_spike2)

        records.append({
            "library.ID": sample,
            "Condition": condition_base,
            "Biorep": biorep,
            "IP": ip,
            "TechRep": techrep,
            f"{target_genome}_reads": count_target,
            f"{spike1_genome}_reads": count_s1,
            f"{spike2_genome}_reads": count_s2,
            f"{spike1_genome}/{target_genome}": count_s1 / count_target if count_target > 0 else 0,
            f"{spike2_genome}/{target_genome}": count_s2 / count_target if count_target > 0 else 0
        })

    seqstats = pd.DataFrame(records)
 #   print(seqstats)
    if seqstats.empty:
        raise RuntimeError(f"No valid samples parsed in {target_dir}")

    # ================== QC Flags ==================
    s1_ratio = f"{spike1_genome}/{target_genome}"
    s2_ratio = f"{spike2_genome}/{target_genome}"

    def flag_low_high(row):
        flags = []
        if row["IP"] == "input":
            if row[s1_ratio] < 0.001 or row[s2_ratio] < 0.001:
                flags.append("low_spike_ratio")
            if row[s1_ratio] > 0.25 or row[s2_ratio] > 0.25:
                flags.append("high_spike_ratio")
        return ";".join(flags) if flags else ""

    seqstats["QC_flag"] = seqstats.apply(flag_low_high, axis=1)

    # ================== Input Ratios ==================
    inputs = seqstats.query("IP == 'input'")[["Condition", "Biorep", s1_ratio, s2_ratio]]

    inputs = inputs.rename(columns={
        s1_ratio: f"{spike1_genome}_input",
        s2_ratio: f"{spike2_genome}_input",
    })

    # Merge with overwrite of duplicate columns
    seqstats = seqstats.merge(inputs, on=["Condition", "Biorep"], how="left", suffixes=("", "_old"))

    # Drop any old (duplicate) columns from previous runs
    for col in seqstats.columns:
        if col.endswith("_old"):
            seqstats.drop(columns=[col], inplace=True)

    # Compute ratios
    seqstats[f"{spike1_genome} IP/input"] = (
        seqstats[s1_ratio] / seqstats[f"{spike1_genome}_input"]
    )

    seqstats[f"{spike2_genome} IP/input"] = (
        seqstats[s2_ratio] / seqstats[f"{spike2_genome}_input"]
    )

    # ---------------------------------------
    # Determine control samples
    # ---------------------------------------
    
    if "Control" in df_meta.columns:
        # Use control column from metadata
        df_meta["Control"] = df_meta["Control"].astype(bool)
        control_samples = set(df_meta.loc[df_meta["Control"], "library.ID"])
        print(f"Identified {len(control_samples)} control samples from metadata.")
    
        # restrict seqstats to those samples
        control_subset = seqstats[seqstats["library.ID"].isin(control_samples)]
    
    else:
        # No "Control" column → use ALL IP samples as control group
        print("\nNo 'Control' column found; using ALL IP samples as control.")
        control_subset = seqstats.query("IP != 'input'")


    # Compute control means
    control_means = (
        control_subset
        .groupby("IP")[
            [f"{spike1_genome} IP/input", f"{spike2_genome} IP/input"]
        ]
        .mean()
        .rename(columns={
            f"{spike1_genome} IP/input": f"{spike1_genome}_control_avg",
            f"{spike2_genome} IP/input": f"{spike2_genome}_control_avg"
        })
        .reset_index()
    )

    # Merge back
    seqstats = seqstats.merge(control_means, on="IP", how="left")


    # Normalized values
    seqstats[f"{spike1_genome} IP/input control averaged"] = (
        seqstats[f"{spike1_genome} IP/input"] / seqstats[f"{spike1_genome}_control_avg"]
    )
    seqstats[f"{spike2_genome} IP/input control averaged"] = (
        seqstats[f"{spike2_genome} IP/input"] / seqstats[f"{spike2_genome}_control_avg"]
    )

    # ================== Merge Back to Metadata ==================

    # ---------------- Clean seqstats before merging to prevent duplicate columns ----------------
    merge_key = "library.ID"
    
    # Columns to exclude from seqstats merge because they exist in metadata
    cols_to_exclude = ["Biorep", "IP", "TechRep", "Condition"]

    # Keep library.ID + all other columns not in cols_to_exclude
    cols_to_merge = [merge_key] + [c for c in seqstats.columns if c not in cols_to_exclude]
    
    # Remove duplicates while preserving order
    seen = set()
    cols_to_merge = [x for x in cols_to_merge if not (x in seen or seen.add(x))]
    
    # ---------------- now merge ----------------
    df_updated = df_meta.merge(seqstats[cols_to_merge], on=merge_key, how="left")

    df_updated.to_csv("sample_metadata.tsv", sep="\t", index=False)
    print("\nUpdated sample_metadata.tsv written.")

    print("\n=== Updated metadata preview ===")
    print(df_updated.head().to_string(index=False))


# ============================== CLI ==============================
def main():
    parser = argparse.ArgumentParser(
        description="Update sample_metadata.tsv with sequencing stats and spike-in ratios."
    )
    parser.add_argument("--target_genome", required=True,
                        help="Primary genome (e.g., hg38)")
    parser.add_argument("--spike1_genome", required=True,
                        help="First spike-in genome (e.g., dm6)")
    parser.add_argument("--spike2_genome", required=True,
                        help="Second spike-in genome (e.g., sac3)")
    parser.add_argument("--output_dir", type=Path, default=Path.cwd(),
                        help="Working directory (default: current working directory)")
    parser.add_argument("--samtools_path", default="samtools",
                        help="Path to samtools executable")
    parser.add_argument("--metadata", default="./sample_names.tsv",
                        help="Path to metadata file")
    parser.add_argument("--force_overwrite", action="store_true",
    help="Force all steps to run, even if output files already exist.")

    args = parser.parse_args()

    update_sample_metadata(
        target_genome=args.target_genome,
        spike1_genome=args.spike1_genome,
        spike2_genome=args.spike2_genome,
        output_dir=args.output_dir,
        samtools_path=args.samtools_path,
        metadata=args.metadata,
        force_overwrite=args.force_overwrite,
    )


if __name__ == "__main__":
    main()
