#!/usr/bin/env python3
import subprocess
from pathlib import Path
import argparse
import pandas as pd
import numpy as np

# ======================== Wrapper Function ========================

def estimate_spikein(
    user_dir: Path = None,
    metadata_file: Path = None,
    spike1_species: str = "dm6",
    spike2_species: str = "sac3",
    target_species: str = "hg38",
    SNR_region: str = "tss",
    frag_length: int = 150,
    hist_size: int = 4000,
    hist_bin: int = 25,
    start_position: int = -100,
    end_position: int = 700,
    genome_dirs: list = None,
    save_file: Path = None  # <--- new argument
):
    if user_dir is None:
        user_dir = Path.cwd()
    else:
        user_dir = Path(user_dir)

    if metadata_file is None:
        metadata_file = user_dir / "sample_metadata.tsv"

    if save_file is None:
        save_file = user_dir / "sample_metadata.norm.tsv"

    # --- Step 1: tagdirs & histograms ---
    make_tagdirs(
        user_dir=user_dir,
        spike1_species=spike1_species,
        spike2_species=spike2_species,
        target_species=target_species,
        SNR_region=SNR_region,
        frag_length=frag_length,
        hist_size=hist_size,
        hist_bin=hist_bin
    )

    # --- Step 2: Read histograms & metadata ---
    dm6_hist_file = user_dir / f"{spike1_species}_data/{spike1_species}_histograms/hist_{SNR_region}_{spike1_species}.txt"
    sac3_hist_file = user_dir / f"{spike2_species}_data/{spike2_species}_histograms/hist_{SNR_region}_{spike2_species}.txt"

    hist_dm6 = pd.read_csv(dm6_hist_file, sep="\t", comment="#")
    hist_sac3 = pd.read_csv(sac3_hist_file, sep="\t", comment="#")
    seqstats = pd.read_csv(metadata_file, sep="\t")

    # --- Step 3: Process histograms & compute AUC ---
    dm6_signal = compute_auc(process_histograms(hist_dm6, spike1_species),
                             start_position=start_position,
                             end_position=end_position)
    sac3_signal = compute_auc(process_histograms(hist_sac3, spike2_species),
                              start_position=start_position,
                              end_position=end_position)
    dm6_signal["library.ID"] = dm6_signal["library.ID"].apply(lambda x: Path(x).name)
    sac3_signal["library.ID"] = sac3_signal["library.ID"].apply(lambda x: Path(x).name)

    print("dm6_signal library.IDs:", dm6_signal['library.ID'].tolist())
    print("sac3_signal library.IDs:", sac3_signal['library.ID'].tolist())

    # --- Step 4: Merge with metadata & normalize ---
    seqstats_norm = normalize_and_merge(seqstats, dm6_signal, sac3_signal)

    # --- Step 5: Always save ---
    try:
        print("\n=== Inspect seqstats_norm before saving ===")
        print("Columns:", seqstats_norm.columns.tolist())
        print("Number of rows:", seqstats_norm.shape[0])
        print("First 10 rows:\n", seqstats_norm.head(10))
        seqstats_norm.to_csv(save_file, sep="\t", index=False)
        print(f"Done! Saved normalized sample metadata to {save_file}")
    except Exception as e:
        print("ERROR saving seqstats_norm:", e)
        raise

    return seqstats_norm

    
# ======================== Core Functions ========================

def make_tagdirs(
    user_dir: Path,
    spike1_species: str,
    spike2_species: str,
    target_species: str,
    SNR_region: str = "tss",
    homer_path: str = "HOMER",
    frag_length: int = 150,
    hist_size: int = 4000,
    hist_bin: int = 25,
):
    """
    Make HOMER tag directories and generate TSS histograms for spike-in species.
    """
    spike_names = [spike1_species, spike2_species]
    homer_genome_dict = {spike1_species: "dm6", spike2_species: "sacCer3"}

    for spike in spike_names:
        aligned_dir = user_dir / f"{spike}_data" / f"{spike}_aligned"
        out_tagdir = user_dir / f"{spike}_data" / f"{spike}_tagdirs"
        out_hist = user_dir / f"{spike}_data" / f"{spike}_histograms"

        out_tagdir.mkdir(parents=True, exist_ok=True)
        out_hist.mkdir(parents=True, exist_ok=True)

        sam_files = sorted(aligned_dir.glob("*.nosuffx2.sam"))
        if not sam_files:
            print(f"No SAM files found for {spike} in {aligned_dir}, skipping.")
            continue

        for sam in sam_files:
            base = sam.stem.replace(".nosuffx2", "")
            tagdir = out_tagdir / f"{base}-tagdir"
            if not tagdir.exists():
                cmd = [
                    "makeTagDirectory", str(tagdir), str(sam),
                    "-genome", homer_genome_dict[spike],
                    "-fragLength", str(frag_length),
                    "-checkGC"
                ]
                print(f"Creating tag directory for {base} ({spike})")
                subprocess.run(cmd, check=True)
            else:
                print(f"Skipping existing tagdir: {tagdir}")

        # Create histogram
        tagdirs = sorted(out_tagdir.glob("*-tagdir"))
        if tagdirs:
            output_hist = out_hist / f"hist_{SNR_region}_{spike}.txt"
            if not output_hist.exists():
                cmd = [
                    "annotatePeaks.pl", SNR_region, homer_genome_dict[spike],
                    "-size", str(hist_size), "-hist", str(hist_bin), "-d"
                ] + [str(td) for td in tagdirs]
                print(f"Generating histogram for {spike}")
                with open(output_hist, "w") as out_f:
                    subprocess.run(cmd, stdout=out_f, check=True)
            else:
                print(f"Skipping existing histogram: {output_hist}")

def process_histograms(df: pd.DataFrame, species: str):
    """Convert HOMER histogram to long format and parse sample metadata."""
    df = df.rename(columns={df.columns[0]: "Distance_from_tss"})
    df.columns = (
        df.columns
        .str.replace(r"\.(dm6|sac3|hg38)-tagdir.*$", "", regex=True)
        .str.replace(r"\.Coverage$", "", regex=True)
        .str.replace(r"\.+", "_", regex=True)
    )
    df_long = df.melt(id_vars="Distance_from_tss", var_name="Sample", value_name="Coverage")
    df_long["Sample"] = df_long["Sample"].str.replace(r"\..*$", "", regex=True)

    # Extract sample metadata (naive split)
    parts = df_long["Sample"].str.split("_", expand=True)
    df_long["cell_type"] = parts[0]
    df_long["treatment"] = parts[1] if parts.shape[1] > 1 else "unknown"
    df_long["timepoint"] = parts[2] if parts.shape[1] > 2 else "unknown"
    df_long["biorep"] = parts[3] if parts.shape[1] > 3 else "unknown"
    df_long["antibody"] = parts[4] if parts.shape[1] > 4 else "unknown"

    print(f"[{species}] Parsed {df_long['Sample'].nunique()} samples")
    return df_long

def compute_auc(df: pd.DataFrame, start_position=-100, end_position=700):
    """Compute AUC per sample in the specified region."""
    df = df[(df["Distance_from_tss"] >= start_position) & (df["Distance_from_tss"] <= end_position)]
    auc_list = []
    for sample, group in df.groupby("Sample"):
        if group["Coverage"].notna().any():
            auc = np.trapz(group["Coverage"].astype(float).values,
                           group["Distance_from_tss"].astype(float).values)
        else:
            auc = np.nan
        auc_list.append({"library.ID": sample, "AUC_peaks": auc})
    return pd.DataFrame(auc_list)


def normalize_and_merge(seqstats, 
                        dm6_signal, 
                        sac3_signal):
    """
    Merge metadata + per-genome AUC signal data, harmonize IDs, 
    compute IP efficiency using automatic Input mapping,
    normalize to control samples, and compute adjusted normalization factors.
    """

    import numpy as np
    import pandas as pd
    from pathlib import Path

    # -------------------------
    # Clean up library.IDs
    # -------------------------
    seqstats["library.ID"] = seqstats["library.ID"].astype(str).str.strip()
    dm6_signal["library.ID"] = dm6_signal["library.ID"].astype(str).apply(lambda x: Path(x).name)
    sac3_signal["library.ID"] = sac3_signal["library.ID"].astype(str).apply(lambda x: Path(x).name)

    # -------------------------
    # Merge AUC signals
    # -------------------------
    seqstats = seqstats.merge(dm6_signal, on="library.ID", how="left", suffixes=("", "_dm6"))
    seqstats = seqstats.merge(sac3_signal, on="library.ID", how="left", suffixes=("", "_sac3"))

    # Rename AUC columns
    seqstats = seqstats.rename(columns={
        "AUC_peaks": "fly_K9ac_signal",
        "AUC_peaks_sac3": "yeast_K9ac_signal"
    })

    # -------------------------
    # Automatic Input mapping
    # -------------------------
    seqstats["IP"] = seqstats["IP"].astype(str).str.lower()
    if "Condition" not in seqstats.columns:
        seqstats["Condition"] = "All"  # fallback if no Condition column

    fly_input_map = (
        seqstats[seqstats["IP"] == "input"]
        .groupby("Condition")["fly_K9ac_signal"]
        .mean()
        .to_dict()
    )
    yeast_input_map = (
        seqstats[seqstats["IP"] == "input"]
        .groupby("Condition")["yeast_K9ac_signal"]
        .mean()
        .to_dict()
    )

    seqstats["fly_input_value"] = seqstats["Condition"].map(fly_input_map)
    seqstats["yeast_input_value"] = seqstats["Condition"].map(yeast_input_map)

    # -------------------------
    # Compute IP efficiency
    # -------------------------
    seqstats["fly_ip_efficiency"] = seqstats["fly_K9ac_signal"] - seqstats["fly_input_value"]
    seqstats["yeast_ip_efficiency"] = seqstats["yeast_K9ac_signal"] / seqstats["yeast_input_value"]

    # -------------------------
    # Normalize within antibody groups using control samples
    # -------------------------
    seqstats["antibody"] = seqstats["library.ID"].str.extract(r"_(H3[^_]+)_")[0]

    for ab, group in seqstats.groupby("antibody"):
        print(f"\nProcessing antibody group: {ab}")

        # Identify control samples from Control column
        control_mask = group["Control"].astype(str).str.upper() == "TRUE"
        control_samples = group.loc[control_mask]

        if control_mask.any():
            fly_ctrl_mean = np.nanmean(control_samples["fly_ip_efficiency"])
            yeast_ctrl_mean = np.nanmean(control_samples["yeast_ip_efficiency"])
            print(f"{ab}: Using Control column TRUE samples as control.")
        else:
            fly_ctrl_mean = np.nanmean(group["fly_ip_efficiency"])
            yeast_ctrl_mean = np.nanmean(group["yeast_ip_efficiency"])
            print(f"{ab}: No Control==TRUE samples found. Using antibody group average.")

        # Safety check
        if fly_ctrl_mean == 0 or yeast_ctrl_mean == 0:
            print(f"{ab}: Warning: control mean is zero, skipping normalization.")
            continue

        # Apply normalization
        seqstats.loc[group.index, "fly_ip_efficiency_norm"] = group["fly_ip_efficiency"] / fly_ctrl_mean
        seqstats.loc[group.index, "yeast_ip_efficiency_norm"] = group["yeast_ip_efficiency"] / yeast_ctrl_mean

        print(f"{ab}: Normalized fly/yeast IP efficiencies (fly_ctrl_mean={fly_ctrl_mean:.4f}, yeast_ctrl_mean={yeast_ctrl_mean:.4f})")

    # -------------------------
    # Compute adjusted normalization factors
    # -------------------------
    seqstats["dm6.normfactor.ipeff.adj"] = seqstats["fly_ip_efficiency_norm"] * seqstats.get("dm6 IP/input control averaged", 1)
    seqstats["sac3.normfactor.ipeff.adj"] = seqstats["yeast_ip_efficiency_norm"] * seqstats.get("sac3 IP/input control averaged", 1)
    seqstats["dual.normfactor.ipeff.adj"] = seqstats[["dm6.normfactor.ipeff.adj", "sac3.normfactor.ipeff.adj"]].mean(axis=1)

    print("Success: Added normalized IP efficiency and adjusted spike-in columns.")

    return seqstats


# ======================== Main Workflow ========================
def main():
    parser = argparse.ArgumentParser(
        description="Estimate spike-in IP efficiency using HOMER tag directories and histograms."
    )
    parser.add_argument("--user_dir", type=Path, default=Path.cwd())
    parser.add_argument("--target_species", required=True)
    parser.add_argument("--spike1_species", default="dm6")
    parser.add_argument("--spike2_species", default="sac3")
    parser.add_argument("--SNR_region", default="tss")
    parser.add_argument("--frag_length", type=int, default=150)
    parser.add_argument("--hist_size", type=int, default=4000)
    parser.add_argument("--hist_bin", type=int, default=25)
    parser.add_argument("--start_position", type=int, default=-100)
    parser.add_argument("--end_position", type=int, default=700)
    args = parser.parse_args()

    print("=== SCRIPT START ===", flush=True)

    # Path to input metadata
    metadata_file = args.user_dir / "sample_metadata.tsv"
    if not metadata_file.exists():
        raise FileNotFoundError(f"{metadata_file} not found.")

    # Path to output normalized metadata (optional override)
    save_file = args.user_dir / "sample_metadata.norm.tsv"

    # Run the wrapper (always saves internally)
    print("STEP 6: Estimate spike-in IP efficiency and normalize...", flush=True)
    try:
        estimate_spikein(
            user_dir=args.user_dir,
            metadata_file=metadata_file,
            spike1_species=args.spike1_species,
            spike2_species=args.spike2_species,
            target_species=args.target_species,
            SNR_region=args.SNR_region,
            frag_length=args.frag_length,
            hist_size=args.hist_size,
            hist_bin=args.hist_bin,
            start_position=args.start_position,
            end_position=args.end_position,
            save_file=save_file
        )
    except Exception as e:
        print("ERROR in estimate_spikein:", e, flush=True)
        raise

    print("=== SCRIPT END ===", flush=True)

if __name__ == "__main__":
    main()