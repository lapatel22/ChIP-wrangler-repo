#!/usr/bin/env python3
"""
05_get_sequencing_stats.py
Update sample_metadata.tsv with sequencing stats, spike-in ratios, QC flags,
and normalized values based on BAM files in hg38, dm6, and sac3 directories.
"""

import pandas as pd
import subprocess
import re
from pathlib import Path

def update_sample_metadata(
    metadata_file: str = None,
    target_dir: str = None,
    dm6_dir: str = None,
    sac3_dir: str = None,
    samtools_path: str = "samtools",
    control_conditions: list = None
):
    """
    Update sample_metadata.tsv with sequencing stats, spike-in ratios, QC flags,
    and normalized values.

    Parameters:
        metadata_file: path to existing sample_metadata.tsv
        target_dir: path to BAMs for the main genome (hg38)
        dm6_dir: path to dm6 BAMs
        sac3_dir: path to sac3 BAMs
        samtools_path: path to samtools executable
        control_conditions: list of control condition names for normalization
    """
    # Base directory is repo root
    base_dir = Path(__file__).parent.parent

    # Build paths if not provided
    if metadata_file is None:
        metadata_file = base_dir / "../sample_metadata.tsv"
    if target_dir is None:
        target_dir = base_dir / "../hg38_data/hg38_aligned"
    if dm6_dir is None:
        dm6_dir = base_dir / "../dm6_data/dm6_aligned"
    if sac3_dir is None:
        sac3_dir = base_dir / "../sac3_data/sac3_aligned"

    metadata_file = Path(metadata_file)
    target_dir = Path(target_dir)
    dm6_dir = Path(dm6_dir)
    sac3_dir = Path(sac3_dir)

    if not metadata_file.exists():
        raise FileNotFoundError(f"Cannot find sample_metadata.tsv at {metadata_file}")

    if control_conditions is None:
        control_conditions = ["HelaS3_100sync_0inter"]

    print(f"Reading metadata from: {metadata_file}")
    df_meta = pd.read_csv(metadata_file, sep="\t")

    bam_files = list(target_dir.glob("*.bam"))
    if not bam_files:
        raise FileNotFoundError(f"No BAM files found in {target_dir}")

    records = []

    for bam in bam_files:
        sample_full = bam.stem
        sample = re.sub(r"\.(hg38|dm6|sac3)$", "", sample_full)

        # expected pattern: Cell_sync_inter_biorep_IP_techrep
        match = re.match(r"^([A-Za-z0-9]+_[0-9]+sync_[0-9]+inter)_(\d+)_(\w+)_([0-9]+)$", sample)
        if not match:
            print(f"⚠ Skipping unrecognized sample name: {sample_full}")
            continue

        condition_base, biorep, ip, techrep = match.groups()
        bam_dm6 = dm6_dir / f"{sample}.dm6.bam"
        bam_sac3 = sac3_dir / f"{sample}.sac3.bam"

        if not (bam.exists() and bam_dm6.exists() and bam_sac3.exists()):
            print(f"⚠ Missing BAM(s) for sample {sample} — skipping.")
            continue

        # Count aligned reads
        def count_reads(path):
            result = subprocess.run([samtools_path, "view", "-c", str(path)],
                                    capture_output=True, text=True)
            return int(result.stdout.strip() or 0)

        count_target = count_reads(bam)
        count_dm6 = count_reads(bam_dm6)
        count_sac3 = count_reads(bam_sac3)

        records.append({
            "library.ID": sample,
            "Condition": condition_base,
            "Biorep": biorep,
            "IP": ip,
            "TechRep": techrep,
            "hg38_reads": count_target,
            "dm6_reads": count_dm6,
            "sac3_reads": count_sac3,
            "dm6/hg38": count_dm6 / count_target if count_target > 0 else 0,
            "sac3/hg38": count_sac3 / count_target if count_target > 0 else 0
        })

    stats_df = pd.DataFrame(records)

    if stats_df.empty:
        raise RuntimeError(f"No valid samples parsed in {target_dir}")

    # QC flag
    def flag_low_high(row):
        flags = []
        if row["IP"] == "input":
            if row["dm6/hg38"] < 0.001 or row["sac3/hg38"] < 0.001:
                flags.append("low_spike_ratio")
            if row["dm6/hg38"] > 0.25 or row["sac3/hg38"] > 0.25:
                flags.append("high_spike_ratio")
        return ";".join(flags) if flags else ""

    stats_df["QC_flag"] = stats_df.apply(flag_low_high, axis=1)

    # IP/input ratios
    inputs = stats_df.query("IP == 'input'")[["Condition", "Biorep", "dm6/hg38", "sac3/hg38"]]
    inputs = inputs.rename(columns={"dm6/hg38": "dm6_input", "sac3/hg38": "sac3_input"})
    stats_df = stats_df.merge(inputs, on=["Condition", "Biorep"], how="left")
    stats_df["dm6 IP/input"] = stats_df["dm6/hg38"] / stats_df["dm6_input"]
    stats_df["sac3 IP/input"] = stats_df["sac3/hg38"] / stats_df["sac3_input"]

    # Control averages
    control_means = (
        stats_df
        .query("Condition in @control_conditions and IP != 'input'")
        .groupby("IP")[["dm6 IP/input", "sac3 IP/input"]]
        .mean()
        .rename(columns={
            "dm6 IP/input": "dm6_control_avg",
            "sac3 IP/input": "sac3_control_avg"
        })
        .reset_index()
    )

    # Merge back to all samples based on IP type only
    stats_df = stats_df.merge(control_means, on="IP", how="left")

    # Compute normalization to control averages
    stats_df["dm6 IP/input control averaged"] = (
        stats_df["dm6 IP/input"] / stats_df["dm6_control_avg"]
    )
    stats_df["sac3 IP/input control averaged"] = (
        stats_df["sac3 IP/input"] / stats_df["sac3_control_avg"]
    )

    # Merge into metadata
    df_updated = df_meta.merge(stats_df, on="library.ID", how="left")

    # Save updated metadata
    df_updated.to_csv(metadata_file, sep="\t", index=False)
    print(f" Updated sample_metadata.tsv written with new columns.")

    # Optional preview
    print("\n=== Updated metadata preview ===")
    print(df_updated.head().to_string(index=False))

    return df_updated


if __name__ == "__main__":
    update_sample_metadata()
