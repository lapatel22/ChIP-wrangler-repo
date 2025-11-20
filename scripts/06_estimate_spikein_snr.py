#!/usr/bin/env python3
import subprocess
from pathlib import Path
import argparse
import pandas as pd
import numpy as np

# ======================== Config ========================

def make_tagdirs(
    user_dir: str,
    target_species: str,
    spike1_species: str,
    spike2_species: str,
    SNR_region: str = "tss",
    homer_path: str = "HOMER",
    frag_length: int = 150,
    hist_size: int = 4000,
    hist_bin: int = 25,
    experiment_id: str = "LP78"
):

    spike1_outdir = user_dir / f"{spike1_species}_data" / f"{spike1_species}_tagdirs"
    spike2_outdir = user_dir / f"{spike2_species}_data" / f"{spike2_species}_tagdirs"
    spike1_histogramdir = user_dir / f"{spike1_species}_data" / f"{spike1_species}_histograms"
    spike2_histogramdir = user_dir / f"{spike2_species}_data" / f"{spike2_species}_histograms"
    
    spike_names = [spike1_species, spike2_species]
    spike_outdirs = [spike1_outdir, spike2_outdir]
    dictionary = {spike1_species: "dm6", spike2_species: "sacCer3"}

    # HOMER binaries
    make_tag_dir_bin = "makeTagDirectory"
    annotate_peaks_bin = "annotatePeaks.pl"

    # ======================== Step 1: Make tag directories ========================
    for spike_name, spike_outdir in zip(spike_names, spike_outdirs):
    
        aligned_dir = user_dir / f"{spike_name}_data" / f"{spike_name}_aligned"    
        spike_outdir.mkdir(parents=True, exist_ok=True)

        sam_files = sorted(aligned_dir.glob("*.nosuffx2.sam"))
        if not sam_files:
            print(f"No SAM files found in {aligned_dir}, skipping {spike_name}.")
            continue

        for sam in sam_files:
            base = sam.stem.replace(".nosuffx2", "")
            tagdir = spike_outdir / f"{base}-tagdir"
            if tagdir.exists():
                print(f"Skipping {tagdir} — already exists.")
                continue

            cmd = [
                make_tag_dir_bin,
                str(tagdir),
                str(sam),
                "-genome", dictionary[spike_names],
                "-fragLength", str(frag_length),
                "-checkGC"
            ]

            print(f"Success: Creating tag directory for {base}")
            subprocess.run(cmd, check=True)

    # ======================== Step 2: Generate TSS metagene histogram ========================
    for spike_name, spike_outdir in zip(spike_names, spike_outdirs):
        tagdirs = sorted(spike_outdir.glob("*-tagdir"))

        if not tagdirs:
            print(f"Warning: No tag directories found for {spike_name}, skipping TSS histogram.")
            continue

        hist_dir = user_dir / f"{spike_name}_data" / f"{spike_name}_histograms"
        hist_dir.mkdir(parents=True, exist_ok=True)

        output_hist = hist_dir / f"hist_{SNR_region}_{spike_name}_{experiment_id}.txt"
        if output_hist.exists():
            print(f"Skipping TSS histogram — {output_hist} already exists.")
            continue

        cmd = [
            annotate_peaks_bin,
            str(SNR_region),
            dictionary[spike_name],
            "-size", str(hist_size),
            "-hist", str(hist_bin),
            "-d"
        ] + [str(td) for td in tagdirs]

        print(f"Success: Generating TSS histogram for {base}")
        with open(output_hist, "w") as out_f:
            subprocess.run(cmd, stdout=out_f, check=True)

    print("All tag directories and per-genome TSS histograms completed.")

# ======================== Helper Functions ========================

def process_histograms(df, species):
    """
    Convert HOMER histogram to long format and extract sample metadata.
    """
    df = df.rename(columns={df.columns[0]: "Distance_from_tss"})

    # Remove trailing tagdir suffix and any coverage suffix
    df.columns = (
        df.columns
        .str.replace(r"\.(dm6|sac3|hg38)-tagdir.*$", "", regex=True)
        .str.replace(r"\.Coverage$", "", regex=True)
        .str.replace(r"\.+", "_", regex=True)
    )

    # Melt to long format
    df_long = df.melt(id_vars="Distance_from_tss", var_name="Sample", value_name="Coverage")

    # Clean up sample names: keep only what’s before the first '.' or any extra text
    df_long["Sample"] = df_long["Sample"].str.replace(r"\..*$", "", regex=True)

    # Extract sample metadata
    parts = df_long["Sample"].str.split("_", expand=True)
    if parts.shape[1] >= 5:
        df_long["cell_type"] = parts[0]
        df_long["treatment"] = parts[1]
        df_long["timepoint"] = parts[2]
        df_long["biorep"] = parts[3]
        df_long["antibody"] = parts[4]
    elif parts.shape[1] >= 4:
        df_long["cell_type"] = parts[0]
        df_long["treatment"] = parts[1]
        df_long["timepoint"] = parts[2]
        df_long["biorep"] = parts[3]
        df_long["antibody"] = "unknown"
    else:
        print("Warning: Unexpected sample name format in histogram data.")

    print(f"\n[{species}] Parsed {df_long['Sample'].nunique()} samples from histogram.")
    return df_long

def compute_auc(df):
    """
    Compute AUC within -100 to +700 bp from TSS for each sample.
    """
    df = df[(df["Distance_from_tss"] >= -100) & (df["Distance_from_tss"] <= 700)]

    auc_list = []
    for sample, group in df.groupby("Sample"):
        if group["Coverage"].notna().any():
            x = group["Distance_from_tss"].astype(float).values
            y = group["Coverage"].astype(float).values
            auc = np.trapz(y, x)
        else:
            auc = np.nan
        auc_list.append({"library.ID": sample, "AUC_peaks": auc})

    auc_df = pd.DataFrame(auc_list)
    print(f"Success: Computed AUC for {auc_df.shape[0]} samples (non-NA: {(auc_df['AUC_peaks'].notna()).sum()})")
    return auc_df
    

def normalize_and_merge(seqstats, dm6_signal, sac3_signal, control_pattern="HelaS3_100sync_0inter_"):
    """
    Merge metadata + per-genome AUC signal data, harmonize IDs, 
    compute IP efficiency (fly: subtract input, yeast: divide by input),
    normalize to control means, and compute adjusted normalization factors.
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
    # Map input values to all IP samples
    # -------------------------
    # Build mapping from condition -> input value
    fly_input_map = seqstats[seqstats["IP"].str.lower() == "input"].set_index("Condition")["fly_K9ac_signal"].to_dict()
    yeast_input_map = seqstats[seqstats["IP"].str.lower() == "input"].set_index("Condition")["yeast_K9ac_signal"].to_dict()

    # Map input values to all rows
    seqstats["fly_input_value"] = seqstats["Condition"].map(fly_input_map)
    seqstats["yeast_input_value"] = seqstats["Condition"].map(yeast_input_map)

    # -------------------------
    # Compute IP efficiency
    # -------------------------
    seqstats["fly_ip_efficiency"] = seqstats["fly_K9ac_signal"] - seqstats["fly_input_value"]
    seqstats["yeast_ip_efficiency"] = seqstats["yeast_K9ac_signal"] / seqstats["yeast_input_value"]

    # Initialize normalized columns
    seqstats["fly_ip_efficiency_norm"] = np.nan
    seqstats["yeast_ip_efficiency_norm"] = np.nan

    # -------------------------
    # Normalize within antibody groups
    # -------------------------
    # Extract antibody name
    seqstats["antibody"] = seqstats["library.ID"].str.extract(r"_(H3[^_]+)_")

    for ab, group in seqstats.groupby("antibody"):
        control_mask = (
            group["library.ID"].str.contains(control_pattern, regex=True, na=False)
            & group["library.ID"].str.contains(ab, regex=True, na=False)
        )

        control_candidates = group.loc[control_mask, "library.ID"].tolist()
        print(f"{ab}: control candidates = {control_candidates}")

        if not control_mask.any():
            print(f"{ab}: Warning: No controls matching '{control_pattern}' found. Skipping normalization.\n")
            continue

        fly_ctrl_mean = np.nanmean(group.loc[control_mask, "fly_ip_efficiency"])
        yeast_ctrl_mean = np.nanmean(group.loc[control_mask, "yeast_ip_efficiency"])

        if np.isnan(fly_ctrl_mean) or np.isnan(yeast_ctrl_mean) or fly_ctrl_mean == 0 or yeast_ctrl_mean == 0:
            print(f"{ab}: Warning: Invalid control mean(s). Skipping normalization.\n")
            continue

        seqstats.loc[group.index, "fly_ip_efficiency_norm"] = group["fly_ip_efficiency"] / fly_ctrl_mean
        seqstats.loc[group.index, "yeast_ip_efficiency_norm"] = group["yeast_ip_efficiency"] / yeast_ctrl_mean

        print(f"{ab}: Success: Normalized fly and yeast IP efficiencies (control mean fly={fly_ctrl_mean:.4f}, yeast={yeast_ctrl_mean:.4f})")

    # -------------------------
    # Compute adjusted normalization factors
    # -------------------------
    seqstats["dm6.normfactor.ipeff.adj"] = seqstats["fly_ip_efficiency_norm"] * seqstats["dm6 IP/input control averaged"]
    seqstats["sac3.normfactor.ipeff.adj"] = seqstats["yeast_ip_efficiency_norm"] * seqstats["sac3 IP/input control averaged"]

    seqstats["dual.normfactor.ipeff.adj"] = seqstats[["dm6.normfactor.ipeff.adj", "sac3.normfactor.ipeff.adj"]].mean(axis=1)

    print("Success: Added columns: fly_ip_efficiency, yeast_ip_efficiency,")
    print("                 fly_ip_efficiency_norm, yeast_ip_efficiency_norm,")
    print("                 dm6.normfactor.ipeff.adj, sac3.normfactor.ipeff.adj,")
    print("                 dual.normfactor.ipeff.adj")

    return seqstats




# --- Clean up seqstats_TSS_histone_signal.tsv columns ---
def clean_seqstats_columns(df):
    # Drop redundant or empty columns
    drop_cols = [
        "antibody", "fly_snr_ratio", "yeast_snr_ratio",
        "fly_snr_ratio.adj", "yeast_snr_ratio.adj",
        "condition"
    ]
    df = df.drop(columns=[c for c in drop_cols if c in df.columns], errors="ignore")

    # Fix malformed concatenated columns
    df.columns = df.columns.str.replace(" ", "_")

    # Ensure consistent column order
    desired_order = [
        "library.ID", "Condition", "Biorep", "IP", "TechRep",
        "hg38_reads", "dm6_reads", "sac3_reads", "dm6/hg38", "sac3/hg38", "QC_flag",
        "dm6_input", "sac3_input",
        "dm6_IP_input", "sac3_IP_input",
        "dm6_control_avg", "sac3_control_avg",
        "dm6_IP_input_control_avgd", "sac3_IP_input_control_avgd",
        "fly_K9ac_signal", "yeast_K9ac_signal",
        "fly_K9ac_signal.adj", "yeast_K9ac_signal.adj"
    ]
    # Keep any missing columns at the end
    df = df[[c for c in desired_order if c in df.columns] + 
            [c for c in df.columns if c not in desired_order]]

    return df

# ======================== Argparse ========================

def get_args():
    parser = argparse.ArgumentParser(
        description="Generate tag directories, histograms, and normalization factors for spike-in ChIP-seq."
    )

    parser.add_argument("--user_dir", required=True, type=Path,
                        help="Base directory where species_data folders live.")
    parser.add_argument("--target_species", required=True, type=str,
                        help="Primary genome (e.g. hg38).")
    parser.add_argument("--spike1_species", required=True, type=str,
                        help="First spike-in species (e.g. dm6).")
    parser.add_argument("--spike2_species", required=True, type=str,
                        help="Second spike-in species (e.g. sacCer3).")

    # Optional parameters
    parser.add_argument("--SNR_region", default="tss", type=str,
                        help="Region used for SNR estimation.")
    parser.add_argument("--homer_path", default="HOMER", type=str,
                        help="Path to HOMER installation.")
    parser.add_argument("--frag_length", default=150, type=int)
    parser.add_argument("--hist_size", default=4000, type=int)
    parser.add_argument("--hist_bin", default=25, type=int)
    parser.add_argument("--experiment_id", default="LP78", type=str)

    return parser.parse_args()

# ======================== Main ========================

def main():
    args = get_args()

    user_dir = args.user_dir
    target_species = args.target_species
    spike1_species = args.spike1_species
    spike2_species = args.spike2_species
    SNR_region = args.SNR_region
    homer_path = args.homer_path
    frag_length = args.frag_length
    hist_size = args.hist_size
    hist_bin = args.hist_bin
    experiment_id = args.experiment_id

    # ---- Run tagdir + histogram creation ----
    make_tagdirs(
        user_dir=user_dir,
        target_species=target_species,
        spike1_species=spike1_species,
        spike2_species=spike2_species,
        SNR_region=SNR_region,
        homer_path=homer_path,
        frag_length=frag_length,
        hist_size=hist_size,
        hist_bin=hist_bin,
        experiment_id=experiment_id
    )

    # ---- Continue with your existing processing code ----
    metadata_file = user_dir / "sample_metadata.tsv"
    dm6_hist_file = user_dir / f"{spike1_species}_data/{spike1_species}_histograms/hist_{SNR_region}_{spike1_species}_{experiment_id}.txt"
    sac3_hist_file = user_dir / f"{spike2_species}_data/{spike2_species}_histograms/hist_{SNR_region}_{spike2_species}_{experiment_id}.txt"

    if not metadata_file.exists():
        raise FileNotFoundError(f"{metadata_file} not found.")

    print("Reading HOMER histograms...")
    hist_dm6 = pd.read_csv(dm6_hist_file, sep="\t", comment="#")
    hist_sac3 = pd.read_csv(sac3_hist_file, sep="\t", comment="#")
    seqstats = pd.read_csv(metadata_file, sep="\t")

    print("Processing dm6 histograms...")
    hist_dm6_long = process_histograms(hist_dm6, spike1_species)
    dm6_signal = compute_auc(hist_dm6_long)

    print("Processing sacCer3 histograms...")
    hist_sac3_long = process_histograms(hist_sac3, spike2_species)
    sac3_signal = compute_auc(hist_sac3_long)

    print("Merging with metadata and normalizing...")
    seqstats_norm = normalize_and_merge(
        seqstats, dm6_signal, sac3_signal,
        control_pattern="HelaS3_100sync_0inter_"
    )

    dm6_signal.to_csv("dm6_H3K9ac_signal.tsv", sep="\t", index=False)
    sac3_signal.to_csv("sac3_H3K9ac_signal.tsv", sep="\t", index=False)

    seqstats_df = clean_seqstats_columns(seqstats_norm)
    seqstats_df.to_csv("seqstats_TSS_histone_signal.tsv", sep="\t", index=False)
    print(" Cleaned and saved seqstats_TSS_histone_signal.tsv")

    updated_metadata = seqstats_norm[
        [
            "library.ID",
            "fly_ip_efficiency_norm",
            "yeast_ip_efficiency_norm",
            "dm6.normfactor.ipeff.adj",
            "sac3.normfactor.ipeff.adj",
            "dual.normfactor.ipeff.adj",
        ]
    ].copy()

    seqstats["library.ID"] = seqstats["library.ID"].astype(str).str.strip()
    updated_metadata["library.ID"] = updated_metadata["library.ID"].astype(str).str.strip()

    final_metadata = seqstats.merge(updated_metadata, on="library.ID", how="left")
    final_metadata.to_csv(f"{user_dir}/sample_metadata.norm.tsv", sep="\t", index=False)
    print("Success: Updated sample_metadata.norm.tsv with new normalization factors")
    
if __name__ == "__main__":
    main()
