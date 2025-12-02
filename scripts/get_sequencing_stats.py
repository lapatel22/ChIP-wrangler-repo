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
    user_dir: str,
    target_species: str,
    spike1_species: str,
    spike2_species: str,
    samtools_path: str = "samtools",
    control_conditions: list = None
):
    """
    Update sample_metadata.tsv with sequencing stats, spike-in ratios, QC flags,
    and normalization values.

    All required SAMs must follow naming:

        <sample>.<species>.nosuffx2.sam

    The <sample> portion must match:
        Cell_sync_inter_biorep_IP_techrep
    """
    user_dir = Path(user_dir)

    # ---------------- Required input file ----------------
    metadata_file = user_dir / "sample_names.tsv"
    if not metadata_file.exists():
        raise FileNotFoundError(f"Error: sample_names.tsv not found at {metadata_file}")

    df_meta = pd.read_csv(metadata_file, sep="\t")
    if df_meta.empty:
        raise RuntimeError(f"ERROR: sample_names.tsv at {metadata_file} is EMPTY.")

    print(f"Reading metadata from: {metadata_file}")

    # ---------------- Species folders ----------------
    target_dir = user_dir / f"{target_species}_data" / f"{target_species}_aligned"
    spike1_dir = user_dir / f"{spike1_species}_data" / f"{spike1_species}_aligned"
    spike2_dir = user_dir / f"{spike2_species}_data" / f"{spike2_species}_aligned"


    # ---------------- Gather SAMs ----------------
    sam_files = list(target_dir.glob("*.nosuffx2.sam"))
    if not sam_files:
        raise FileNotFoundError(f"Warning: No .nosuffx2.sam files found in {target_dir}")

    records = []

    # ---------------- Process each sample ----------------
    for sam in sam_files:
        sample_full = sam.stem  # includes .<species>
        # remove .species part only for parsing metadata
        # ex: Hela_sync_inter_1_input_1.hg38 -> drop ".hg38"
        sample = re.sub(rf"\.{target_species}.nosuffx2$", "", sample_full)
   #     print(sample)
        # Expected naming pattern
        match = re.match(
            r"^([A-Za-z0-9]+_[0-9]+sync_[0-9]+inter)_(\d+)_(\w+)_([0-9]+)$",
            sample
        )

        if not match:
            print(f"Warning: Skipping unrecognized sample name: {sample_full}")
            continue

        condition_base, biorep, ip, techrep = match.groups()

        # Build corresponding spike-in SAM paths
      #  print(sam)
        sam_target = target_dir / f"{sample}.{target_species}.bam"
        sam_spike1 = spike1_dir / f"{sample}.{spike1_species}.bam"
        sam_spike2 = spike2_dir / f"{sample}.{spike2_species}.bam"

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
            f"{target_species}_reads": count_target,
            f"{spike1_species}_reads": count_s1,
            f"{spike2_species}_reads": count_s2,
            f"{spike1_species}/{target_species}": count_s1 / count_target if count_target > 0 else 0,
            f"{spike2_species}/{target_species}": count_s2 / count_target if count_target > 0 else 0
        })

    seqstats = pd.DataFrame(records)
 #   print(seqstats)
    if seqstats.empty:
        raise RuntimeError(f"No valid samples parsed in {target_dir}")

    # ================== QC Flags ==================
    s1_ratio = f"{spike1_species}/{target_species}"
    s2_ratio = f"{spike2_species}/{target_species}"

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
        s1_ratio: f"{spike1_species}_input",
        s2_ratio: f"{spike2_species}_input",
    })

    # Merge with overwrite of duplicate columns
    seqstats = seqstats.merge(inputs, on=["Condition", "Biorep"], how="left", suffixes=("", "_old"))

    # Drop any old (duplicate) columns from previous runs
    for col in seqstats.columns:
        if col.endswith("_old"):
            seqstats.drop(columns=[col], inplace=True)

    # Compute ratios
    seqstats[f"{spike1_species} IP/input"] = (
        seqstats[s1_ratio] / seqstats[f"{spike1_species}_input"]
    )

    seqstats[f"{spike2_species} IP/input"] = (
        seqstats[s2_ratio] / seqstats[f"{spike2_species}_input"]
    )

    # ================== Control Averages (FROM METADATA COLUMN) ==================

    print("\nSelecting control samples based on metadata column 'Control' == TRUE")

    # Control column may be bool or string ("TRUE"/"FALSE")
    def is_true(x):
        if isinstance(x, bool):
            return x
        if isinstance(x, str):
            return x.strip().upper() == "TRUE"
        return False

    df_updated["Control_flag"] = df_updated["Control"].apply(is_true)

    # Only IP samples should be used for control mean (exclude inputs)
    control_subset = df_updated.query("Control_flag == True and IP != 'input'")

    if control_subset.empty:
        print("\nWARNING: No control samples found in metadata! Using ALL non-input samples instead.\n")
        control_subset = df_updated.query("IP != 'input'")

    # Compute control means per antibody/IP group
    control_means = (
        control_subset
        .groupby("IP")[
            [f"{spike1_species} IP/input", f"{spike2_species} IP/input"]
        ]
        .mean()
        .rename(columns={
            f"{spike1_species} IP/input": f"{spike1_species}_control_avg",
            f"{spike2_species} IP/input": f"{spike2_species}_control_avg"
        })
        .reset_index()
    )

    # Merge back into df_updated
    df_updated = df_updated.merge(control_means, on="IP", how="left")

    # Normalized values
    df_updated[f"{spike1_species} IP/input control averaged"] = (
        df_updated[f"{spike1_species} IP/input"] /
        df_updated[f"{spike1_species}_control_avg"]
    )

    df_updated[f"{spike2_species} IP/input control averaged"] = (
        df_updated[f"{spike2_species} IP/input"] /
        df_updated[f"{spike2_species}_control_avg"]
    )

    # Merge back
    seqstats = seqstats.merge(control_means, on="IP", how="left")


    # Normalized values
    seqstats[f"{spike1_species} IP/input control averaged"] = (
        seqstats[f"{spike1_species} IP/input"] / seqstats[f"{spike1_species}_control_avg"]
    )
    seqstats[f"{spike2_species} IP/input control averaged"] = (
        seqstats[f"{spike2_species} IP/input"] / seqstats[f"{spike2_species}_control_avg"]
    )

    # ================== Merge Back to Metadata ==================
    df_updated = df_meta.merge(seqstats, on="library.ID", how="left")

    df_updated.to_csv("sample_metadata.tsv", sep="\t", index=False)
    print("\nUpdated sample_metadata.tsv written.")

    print("\n=== Updated metadata preview ===")
    print(df_updated.head().to_string(index=False))


# ============================== CLI ==============================
def main():
    parser = argparse.ArgumentParser(
        description="Update sample_metadata.tsv with sequencing stats and spike-in ratios."
    )
    parser.add_argument("--target_species", required=True,
                        help="Primary genome (e.g., hg38)")
    parser.add_argument("--spike1_species", required=True,
                        help="First spike-in genome (e.g., dm6)")
    parser.add_argument("--spike2_species", required=True,
                        help="Second spike-in genome (e.g., sac3)")
    parser.add_argument("--user_dir", type=Path, default=Path.cwd(),
                        help="Working directory (default: current working directory)")
    parser.add_argument("--samtools_path", default="samtools",
                        help="Path to samtools executable")

    args = parser.parse_args()

    update_sample_metadata(
        target_species=args.target_species,
        spike1_species=args.spike1_species,
        spike2_species=args.spike2_species,
        user_dir=args.user_dir,
        samtools_path=args.samtools_path,
    )


if __name__ == "__main__":
    main()
