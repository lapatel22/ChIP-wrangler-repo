#!/usr/bin/env python3
"""
get_sequencing_stats.py (auto-detect spike-ins)

Auto-detect species under <user_dir>/*_data/*_aligned, pick the species
with the most reads as the target and treat others as spike-ins (0-2).
Produces sample_metadata.tsv in user_dir.

Usage (wrapper will call programmatically):
    python get_sequencing_stats.py --user_dir . --samtools_path samtools
"""
from pathlib import Path
import subprocess
import re
import argparse
import pandas as pd
import numpy as np
from typing import List, Dict


def find_species_dirs(user_dir: Path) -> Dict[str, Path]:
    """
    Find species names that have <user_dir>/<species>_data/<species>_aligned.
    Return dict: {species_name: aligned_dir (Path)}
    """
    species_dirs = {}
    for data_dir in user_dir.glob("*_data"):
        if not data_dir.is_dir():
            continue
        # data_dir name is like "hg38_data" or "dm6_data"
        name = data_dir.name
        if not name.endswith("_data"):
            continue
        species = name[:-5]  # strip "_data"
        aligned_dir = data_dir / f"{species}_aligned"
        if aligned_dir.exists() and aligned_dir.is_dir():
            species_dirs[species] = aligned_dir
    return species_dirs


def list_candidate_files(aligned_dir: Path, species: str) -> List[Path]:
    """
    Prefer .nosuffx2.sam files; otherwise prefer *.{species}.bam; otherwise any .bam.
    """
    files = list(aligned_dir.glob("*.nosuffx2.sam"))
    if files:
        return files
    files = list(aligned_dir.glob(f"*.{species}.bam"))
    if files:
        return files
    files = list(aligned_dir.glob("*.bam"))
    return files


def total_reads_for_species(aligned_dir: Path, species: str, samtools_path: str) -> int:
    """
    Sum samtools view -c over candidate files in this species dir.
    """
    files = list_candidate_files(aligned_dir, species)
    total = 0
    for f in files:
        try:
            res = subprocess.run([samtools_path, "view", "-c", str(f)],
                                 capture_output=True, text=True, check=True)
            total += int(res.stdout.strip() or 0)
        except Exception:
            # if samtools fails, skip file
            continue
    return total


def detect_target_and_spikes(user_dir: Path, samtools_path: str):
    species_dirs = find_species_dirs(user_dir)
    if not species_dirs:
        raise RuntimeError(f"No species alignment directories found under {user_dir}")

    totals = {}
    for sp, adir in species_dirs.items():
        totals[sp] = total_reads_for_species(adir, sp, samtools_path)

    # Choose species with max reads as target
    target = max(totals.items(), key=lambda x: x[1])[0]
    spikes = [s for s in species_dirs.keys() if s != target]

    # Keep at most 2 spike-ins (ordered by total reads desc)
    spikes_sorted = sorted(spikes, key=lambda s: totals[s], reverse=True)
    spike1 = spikes_sorted[0] if len(spikes_sorted) >= 1 else None
    spike2 = spikes_sorted[1] if len(spikes_sorted) >= 2 else None

    return target, spike1, spike2, species_dirs


def build_seqstats(user_dir: Path, target: str, spike1: str, spike2: str, species_dirs: Dict[str, Path],
                   samtools_path: str) -> pd.DataFrame:
    """
    Build seqstats DataFrame with per-sample read counts and ratio columns.
    """
    target_dir = species_dirs[target]

    sample_files = list(target_dir.glob(f"*.{target}.bam"))
    if not sample_files:
        raise RuntimeError(f"No sample files found for target {target} in {target_dir}")

    records = []
    for sf in sample_files:
        stem = sf.stem
        sample = re.sub(rf"\.{re.escape(target)}$", "", stem)

        def find_file(species):
            if not species:
                return None
            adir = species_dirs[species]
            for pattern in [f"{sample}.{species}.bam", f"{sample}*.bam"]:
                candidates = list(adir.glob(pattern))
                if candidates:
                    return candidates[0]
            return None

        f_target = find_file(target)
        f_s1 = find_file(spike1)
        f_s2 = find_file(spike2)

        if f_target is None or (spike1 and f_s1 is None) or (spike2 and f_s2 is None):
            print(f"Missing files for {sample} — skipping")
            continue

        def count_reads(f):
            try:
                out = subprocess.run([samtools_path, "view", "-c", str(f)],
                                     capture_output=True, text=True, check=True)
                return int(out.stdout.strip())
            except Exception:
                return 0

        c_t = count_reads(f_target)
        c_s1 = count_reads(f_s1) if f_s1 else 0
        c_s2 = count_reads(f_s2) if f_s2 else 0

        rec = {"library.ID": sample,
               f"{target}_reads": c_t}
        if spike1:
            rec[f"{spike1}_reads"] = c_s1
            rec[f"{spike1}/{target}"] = (c_s1 / c_t) if c_t else np.nan
        if spike2:
            rec[f"{spike2}_reads"] = c_s2
            rec[f"{spike2}/{target}"] = (c_s2 / c_t) if c_t else np.nan

        records.append(rec)

    df = pd.DataFrame(records)
    if df.empty:
        raise RuntimeError("No samples parsed into seqstats (check file naming).")
    return df


def compute_input_and_ip_input(seqstats: pd.DataFrame, spike1: str, spike2: str, target: str) -> pd.DataFrame:
    """
    Compute input group ratios and IP/input ratios.
    """
    s1_col = f"{spike1}/{target}" if spike1 else None
    s2_col = f"{spike2}/{target}" if spike2 else None

    # get inputs (IP == 'input')
    cols = ["Condition", "Biorep"]
    if s1_col:
        cols.append(s1_col)
    if s2_col:
        cols.append(s2_col)

    inputs = seqstats.query("IP == 'input'")[cols].copy()
    rename_map = {}
    if s1_col:
        rename_map[s1_col] = f"{spike1}_input"
    if s2_col:
        rename_map[s2_col] = f"{spike2}_input"
    if rename_map:
        inputs = inputs.rename(columns=rename_map)

    # merge inputs back into seqstats on Condition+Biorep
    seqstats = seqstats.merge(inputs, on=["Condition", "Biorep"], how="left", suffixes=("", "_old"))
    # drop any leftover _old cols if present
    seqstats = seqstats.loc[:, [c for c in seqstats.columns if not c.endswith("_old")]]

    # compute IP/input columns
    if spike1:
        seqstats[f"{spike1} IP/input"] = seqstats[f"{spike1}/{target}"].div(seqstats[f"{spike1}_input"]).replace([np.inf, -np.inf], np.nan)
    if spike2:
        seqstats[f"{spike2} IP/input"] = seqstats[f"{spike2}/{target}"].div(seqstats[f"{spike2}_input"]).replace([np.inf, -np.inf], np.nan)

    return seqstats


def compute_control_averages_and_normalize(df_updated: pd.DataFrame, spike1: str, spike2: str) -> pd.DataFrame:
    """
    Given df_updated (metadata merged with seqstats), compute control averages for each IP
    and append normalized control-averaged columns.
    """
    def is_true(x):
        if isinstance(x, bool):
            return x
        if isinstance(x, str):
            return x.strip().upper() == "TRUE"
        return False

    if "Control" not in df_updated.columns:
        # create Control column default False
        df_updated["Control"] = False

    df_updated["IP"] = df_updated["IP"].astype(str).str.strip()
    df_updated["Control"] = df_updated["Control"].astype(str).str.strip()
    
    control_subset = df_updated[(df_updated["Control"] == True) & (df_updated["IP"] != "input")]
    
    print("Unique values in IP column:", df_updated["IP"].unique())

    inputs = df_updated[df_updated["IP"] == "input"]
    print("Input rows:", len(inputs))
    print(inputs[["library.ID", "IP"]])
    print(df_updated.filter(regex="reads$").head())
    
    # only IPs (exclude input) for control computation
    if control_subset.empty:
        print("No Control==TRUE rows found; falling back to all non-input IP samples for control averaging.")
        control_subset = df_updated[df_updated["IP"] != "input"]

    print("Control subset is:", control_subset[["library.ID", "Control", "IP"]])
    
    # identify which IP/input columns exist
    avg_cols = []
    if spike1 and f"{spike1} IP/input" in df_updated.columns:
        avg_cols.append(f"{spike1} IP/input")
    if spike2 and f"{spike2} IP/input" in df_updated.columns:
        avg_cols.append(f"{spike2} IP/input")

    if not avg_cols:
        # nothing to average, return unchanged
        print("No IP/input columns present to compute control averages.")
        return df_updated

    control_means = control_subset.groupby("IP")[avg_cols].mean().rename(
        columns={c: c.replace(" IP/input", "_control_avg") for c in avg_cols}
    ).reset_index()

    df_updated = df_updated.merge(control_means, on="IP", how="left")

    # compute normalized columns
    if spike1 and f"{spike1} IP/input" in df_updated.columns:
        df_updated[f"{spike1} IP/input control averaged"] = df_updated[f"{spike1} IP/input"].div(df_updated[f"{spike1}_control_avg"]).replace([np.inf, -np.inf], np.nan)
    if spike2 and f"{spike2} IP/input" in df_updated.columns:
        df_updated[f"{spike2} IP/input control averaged"] = df_updated[f"{spike2} IP/input"].div(df_updated[f"{spike2}_control_avg"]).replace([np.inf, -np.inf], np.nan)

    return df_updated


def update_sample_metadata(
    user_dir: str,
    samtools_path: str = "samtools",
    target_species: str = None,
    spike1_species: str = None,
    spike2_species: str = None
):
    user_dir = Path(user_dir)
    species_dirs = find_species_dirs(user_dir)
    if not species_dirs:
        raise RuntimeError("No species alignment directories found under user_dir")

    # ----------------- Determine target and spike-ins -----------------
    if target_species is None:
        target, spike1, spike2, species_dirs = detect_target_and_spikes(user_dir, samtools_path)
        print(f"Auto-detected target: {target}; spike1: {spike1}; spike2: {spike2}")
    else:
        target = target_species
        spike1 = spike1_species
        spike2 = spike2_species
        print(f"Using specified target: {target}; spike1: {spike1}; spike2: {spike2}")
        
    # build per-sample stats
    # Step 1: Build seqstats with library.ID, target_reads, spike_reads, spike/target ratios
    seqstats = build_seqstats(user_dir, target, spike1, spike2, species_dirs, samtools_path)

    # Step 2: Load metadata
    metadata_file = user_dir / "sample_names.tsv"
    df_meta = pd.read_csv(metadata_file, sep="\t")

    # Strip any whitespace in headers
    df_meta.columns = df_meta.columns.str.strip()

    # Now merge works
    df_merged = df_meta.merge(seqstats, on="library.ID", how="left")

    # Step 4: Compute input ratios on merged df
    df_updated = compute_input_and_ip_input(df_merged, spike1, spike2, target)

    # Normalize Control column to real booleans
    # Normalize Control column → real booleans
    df_updated["Control"] = (
        df_updated["Control"]
        .astype(str)
        .str.strip()
        .str.upper()
        .map({"TRUE": True, "FALSE": False})
    )
    
    # Safety: if any unrecognized values exist → raise clear error
    if df_updated["Control"].isna().any():
        bad = df_updated[df_updated["Control"].isna()]["Control"]
        raise ValueError(
            f"ERROR: Control column contains non-boolean values: {bad.tolist()}\n"
            "Valid values are: TRUE or FALSE"
        )
    
    # If metadata has IP column, refine QC_flag to apply only to true inputs
    if "IP" in df_updated.columns:
        df_updated.loc[df_updated["IP"] != "input", "QC_flag"] = ""

    duplicate_columns = [c for c in df_updated.columns if c.endswith("_y")]

    for col_y in duplicate_columns:
        col_base = col_y[:-2]  # remove "_y"
        col_x = f"{col_base}_x"
    
        if col_x in df_updated.columns:
            # check if values are identical where both exist
            mask = df_updated[col_x].notna() & df_updated[col_y].notna()
            if (df_updated.loc[mask, col_x] == df_updated.loc[mask, col_y]).all():
                # values identical → drop _y
                df_updated.drop(columns=[col_y], inplace=True)
            else:
                # keep _y, but rename _x to avoid confusion
                df_updated.rename(columns={col_x: f"{col_base}_meta"}, inplace=True)
        else:
            # no _x, keep _y as base
            df_updated.rename(columns={col_y: col_base}, inplace=True)

    # IP column is now either IP_x or IP_y, pick metadata if available
    if "IP_x" in df_updated.columns:
        df_updated["IP"] = df_updated["IP_x"]
    elif "IP_y" in df_updated.columns:
        df_updated["IP"] = df_updated["IP_y"]

    # drop any leftover _x/_y
    df_updated = df_updated.drop(columns=[c for c in df_updated.columns if c.endswith(("_x","_y"))])

    # df_updated["IP"] = df_updated["IP_x"]
    # df_updated = df_updated.drop(columns=["IP_x", "IP_y"])

    # compute control averages and normalized values
    df_updated = compute_control_averages_and_normalize(df_updated, spike1, spike2)

    # write out
    outp = user_dir / "sample_metadata.tsv"
    df_updated.to_csv(outp, sep="\t", index=False)
    print(f"Wrote updated metadata to: {outp}")


# -------------------- CLI --------------------
def main():
    parser = argparse.ArgumentParser(description="Auto-detect species and update sample metadata with sequencing stats.")
    parser.add_argument("--user_dir", type=str, default=".", help="Project base directory where *_data/*_aligned exist and sample_names.tsv lives.")
    parser.add_argument("--samtools_path", type=str, default="samtools", help="Path to samtools binary")
    args = parser.parse_args()

    update_sample_metadata(user_dir=args.user_dir, samtools_path=args.samtools_path)


if __name__ == "__main__":
    main()
