#!/usr/bin/env python3
import subprocess
from pathlib import Path
import argparse
import pandas as pd
import numpy as np

# ======================== Wrapper Function ========================

def estimate_spikein(
    output_dir: Path = None,
    metadata_file: Path = None,
    target_species: str = None,
    spike1_species: str=None,
    spike2_species: str=None,
    SNR_region: str = "tss",
    frag_length: int = 150,
    hist_size: int = 4000,
    hist_bin: int = 25,
    start_position: int = -100,
    end_position: int = 700,
    save_file: Path = None,
    force_overwrite: bool = False
):
    if output_dir is None:
        output_dir = Path.cwd()
    else:
        output_dir = Path(output_dir)

    if metadata_file is None:
        metadata_file = output_dir / "sample_metadata.tsv"

    if save_file is None:
        save_file = output_dir / "sample_metadata.norm.tsv"

    # --- Step 1: tagdirs & histograms ---
    make_tagdirs(
        output_dir=output_dir,
        spike1_species=spike1_species,
        spike2_species=spike2_species,
        target_species=target_species,
        SNR_region=SNR_region,
        frag_length=frag_length,
        hist_size=hist_size,
        hist_bin=hist_bin,
        force_overwrite=force_overwrite
    )

    # --- Step 2: Read histograms & metadata ---
    spike1_hist_file = output_dir / f"{spike1_species}_data/{spike1_species}_histograms/hist_{SNR_region}_{spike1_species}.txt"
    spike2_hist_file = output_dir / f"{spike2_species}_data/{spike2_species}_histograms/hist_{SNR_region}_{spike2_species}.txt"

    hist_spike1 = pd.read_csv(spike1_hist_file, sep="\t", comment="#")
    hist_spike2 = pd.read_csv(spike2_hist_file, sep="\t", comment="#")
    seqstats = pd.read_csv(metadata_file, sep="\t")

    # --- Step 3: Process histograms & compute AUC ---
    spike1_signal = compute_auc(process_histograms(hist_spike1, spike1_species),
                                start_position=start_position,
                                end_position=end_position)
    spike2_signal = compute_auc(process_histograms(hist_spike2, spike2_species),
                                start_position=start_position,
                                end_position=end_position)

    # --- Step 4: Merge with metadata & normalize ---
    seqstats_norm = normalize_and_merge(seqstats, spike1_signal, spike2_signal,
                                        spike1_species=spike1_species,
                                        spike2_species=spike2_species)

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
    output_dir: Path,
    target_species: str = None,
    spike1_species: str = None,
    spike2_species: str = None,
    SNR_region: str = "tss",
    homer_path: str = "HOMER",
    frag_length: int = 150,
    hist_size: int = 4000,
    hist_bin: int = 25,
    force_overwrite: bool = False
):
    """
    Make HOMER tag directories and generate TSS histograms for spike-in species.
    """
    spike_names = [spike1_species, spike2_species]
    homer_genome_dict = {spike1_species: spike1_species, spike2_species: spike2_species}

    for spike in spike_names:
        aligned_dir = output_dir / f"{spike}_data" / f"{spike}_aligned"
        out_tagdir = output_dir / f"{spike}_data" / f"{spike}_tagdirs"
        out_hist = output_dir / f"{spike}_data" / f"{spike}_histograms"

        out_tagdir.mkdir(parents=True, exist_ok=True)
        out_hist.mkdir(parents=True, exist_ok=True)

        sam_files = sorted(aligned_dir.glob("*.nosuffx2.sam"))
        if not sam_files:
            print(f"No SAM files found for {spike} in {aligned_dir}, skipping.")
            continue

        for sam in sam_files:
            base = sam.stem.replace(".nosuffx2", "")
            tagdir = out_tagdir / f"{base}-tagdir"
            if not tagdir.exists() or force_overwrite:
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
            if not output_hist.exists() or force_overwrite:
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
    """Convert HOMER histogram to long format and parse sample metadata from metadata file."""
    df = df.rename(columns={df.columns[0]: "Distance_from_tss"})
    df.columns = (
        df.columns
        .str.replace(rf"\.{species}-tagdir.*$", "", regex=True)
        .str.replace(r"\.Coverage$", "", regex=True)
        .str.replace(r"\.+", "_", regex=True)
    )
    df_long = df.melt(id_vars="Distance_from_tss", var_name="Sample", value_name="Coverage")
    df_long["Sample"] = df_long["Sample"].str.replace(r"\..*$", "", regex=True)

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
                        spike1_signal, 
                        spike2_signal,
                        spike1_species,
                        spike2_species):
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
    spike1_signal["library.ID"] = spike1_signal["library.ID"].astype(str).apply(lambda x: Path(x).name)
    spike2_signal["library.ID"] = spike2_signal["library.ID"].astype(str).apply(lambda x: Path(x).name)

    # -------------------------
    # Merge AUC signals
    # -------------------------
    seqstats = seqstats.merge(spike1_signal, on="library.ID", how="left", suffixes=("", f"_{spike1_species}"))
    seqstats = seqstats.merge(spike2_signal, on="library.ID", how="left", suffixes=("", f"_{spike2_species}"))

    # Rename AUC columns to standardized spike-in names
    seqstats = seqstats.rename(columns={
        "AUC_peaks": f"{spike1_species}_signal",
        f"AUC_peaks_{spike2_species}": f"{spike2_species}_signal"
    })

    # -------------------------
    # Clean IP and Condition columns
    # -------------------------
    seqstats["IP"] = seqstats["IP"].astype(str).str.lower().str.strip()
    seqstats["Condition"] = seqstats["Condition"].astype(str).str.strip()

    # -------------------------
    # Automatic Input mapping
    # -------------------------
    spike1_input_map = (
        seqstats.loc[seqstats["IP"] == "input", [f"{spike1_species}_signal", "Condition"]]
        .groupby("Condition")[f"{spike1_species}_signal"]
        .mean()
        .to_dict()
    )
    spike2_input_map = (
        seqstats.loc[seqstats["IP"] == "input", [f"{spike2_species}_signal", "Condition"]]
        .groupby("Condition")[f"{spike2_species}_signal"]
        .mean()
        .to_dict()
    )

    # Map input values to all samples; fill missing with 0 (for subtraction) or 1 (for division)
    seqstats[f"{spike1_species}_input_value"] = seqstats["Condition"].map(spike1_input_map).fillna(0)
    seqstats[f"{spike2_species}_input_value"] = seqstats["Condition"].map(spike2_input_map).fillna(1)

    # -------------------------
    # Compute IP efficiency
    # -------------------------
    seqstats[f"{spike1_species}_ip_efficiency"] = seqstats[f"{spike1_species}_signal"] - seqstats[f"{spike1_species}_input_value"]
    seqstats[f"{spike2_species}_ip_efficiency"] = seqstats[f"{spike2_species}_signal"] / seqstats[f"{spike2_species}_input_value"]

    # -------------------------
    # Normalize within IP (antibody) groups using Control column
    # -------------------------
    # Here "antibody" is taken from IP column
    seqstats["antibody"] = seqstats["IP"]

    for ab, group in seqstats.groupby("antibody"):
        print(f"\nProcessing IP group: {ab}")

        # Identify control samples from Control column
        control_mask = group["Control"].astype(str).str.upper() == "TRUE"
        control_samples = group.loc[control_mask]

        if control_mask.any():
            spike1_ctrl_mean = np.nanmean(control_samples[f"{spike1_species}_ip_efficiency"])
            spike2_ctrl_mean = np.nanmean(control_samples[f"{spike2_species}_ip_efficiency"])
            print(f"{ab}: Using Control==TRUE samples for normalization.")
        else:
            spike1_ctrl_mean = np.nanmean(group[f"{spike1_species}_ip_efficiency"])
            spike2_ctrl_mean = np.nanmean(group[f"{spike2_species}_ip_efficiency"])
            print(f"{ab}: No Control==TRUE samples found. Using antibody group average.")

        # Safety check
        if spike1_ctrl_mean == 0 or spike2_ctrl_mean == 0:
            print(f"{ab}: Warning: control mean is zero, skipping normalization.")
            continue

        # Apply normalization
        seqstats.loc[group.index, f"{spike1_species}_ip_efficiency_norm"] = group[f"{spike1_species}_ip_efficiency"] / spike1_ctrl_mean
        seqstats.loc[group.index, f"{spike2_species}_ip_efficiency_norm"] = group[f"{spike2_species}_ip_efficiency"] / spike2_ctrl_mean
        print(f"{ab}: Normalized {spike1_species}/{spike2_species} IP efficiencies "
              f"({spike1_species}_ctrl_mean={spike1_ctrl_mean:.4f}, {spike2_species}_ctrl_mean={spike2_ctrl_mean:.4f})")

    # -------------------------
    # Compute adjusted normalization factors
    # -------------------------
    seqstats[f"{spike1_species}.normfactor.ipeff.adj"] = seqstats[f"{spike1_species}_ip_efficiency_norm"] * seqstats.get(f"{spike1_species} IP/input control averaged", 1)
    seqstats[f"{spike2_species}.normfactor.ipeff.adj"] = seqstats[f"{spike2_species}_ip_efficiency_norm"] * seqstats.get(f"{spike2_species} IP/input control averaged", 1)
    seqstats["dual.normfactor.ipeff.adj"] = seqstats[[f"{spike1_species}.normfactor.ipeff.adj", f"{spike2_species}.normfactor.ipeff.adj"]].mean(axis=1)

    print("Success: Added normalized IP efficiency and adjusted spike-in columns.")

    return seqstats

# ======================== Main Workflow ========================
def main():
    parser = argparse.ArgumentParser(
        description="Estimate spike-in IP efficiency using HOMER tag directories and histograms."
    )
    parser.add_argument("--output_dir", type=Path, default=Path.cwd())
    parser.add_argument("--target_species", required=True)
    parser.add_argument("--spike1_species", default="dm6")
    parser.add_argument("--spike2_species", default="sacCer3")
    parser.add_argument("--SNR_region", default="tss")
    parser.add_argument("--frag_length", type=int, default=150)
    parser.add_argument("--hist_size", type=int, default=4000)
    parser.add_argument("--hist_bin", type=int, default=25)
    parser.add_argument("--start_position", type=int, default=-100)
    parser.add_argument("--end_position", type=int, default=700)
    parser.add_argument("--force_overwrite", action="store_true")
    args = parser.parse_args()

    print("=== SCRIPT START ===", flush=True)

    metadata_file = args.output_dir / "sample_metadata.tsv"
    if not metadata_file.exists():
        raise FileNotFoundError(f"{metadata_file} not found.")

    save_file = args.output_dir / "sample_metadata.norm.tsv"

    estimate_spikein(
        output_dir=args.output_dir,
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
        save_file=save_file,
        force_overwrite=args.force_overwrite
    )

    print("=== SCRIPT END ===", flush=True)


if __name__ == "__main__":
    main()
