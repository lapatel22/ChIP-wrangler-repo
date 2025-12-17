import subprocess
import pandas as pd
from pathlib import Path
from typing import Optional


def run_homer_peak_calling(
    output_dir: Path,
    target: str,
    metadata_tsv: Path,
    style: str,
    size: Optional[int] = None,
    min_dist: Optional[int] = None,
    force_overwrite: bool = False,
):
    """
    Run HOMER findPeaks using metadata-driven IP ↔ input matching.

    Parameters
    ----------
    output_dir : Path
        Base working directory (user --output_dir)
    metadata_tsv : Path
        Sample metadata TSV
    style : str
        HOMER peak style (e.g. histone, factor)
    size : int, optional
        Peak size (style-dependent default if None)
    min_dist : int, optional
        Minimum distance between peaks (style-dependent default if None)
    force_overwrite : bool
        Overwrite existing peak files
    """

    tagdir_root = output_dir / f"{target}_data" / f"{target}_tagdirs"
    peaks_outdir = output_dir / f"{target}_data"

    peaks_outdir.mkdir(parents=True, exist_ok=True)

    # ------------------------
    # Load metadata
    # ------------------------
    meta = pd.read_csv(metadata_tsv, sep="\t")

    required_cols = {"library.ID", "IP", "Condition", "Biorep"}
    missing = required_cols - set(meta.columns)
    if missing:
        raise ValueError(f"Metadata missing required columns: {missing}")

    # ------------------------
    # Style-based defaults
    # ------------------------
    if style == "histone":
        size = 1000 if size is None else size
        min_dist = 2500 if min_dist is None else min_dist
    elif style == "factor":
        # HOMER defaults → do not pass size/minDist
        pass
    else:
        raise ValueError(f"Unsupported HOMER style: {style}")

    # ------------------------
    # Split IP vs input
    # ------------------------
    ip_meta = meta[meta["IP"] != "input"]
    input_meta = meta[meta["IP"] == "input"]

    if input_meta.empty:
        raise RuntimeError("No input samples found in metadata (IP == 'input').")

    # ------------------------
    # Iterate over IP samples
    # ------------------------
    for _, ip_row in ip_meta.iterrows():

        lib_id = ip_row["library.ID"]
        condition = ip_row["Condition"]
        biorep = ip_row["Biorep"]

        # Match input
        input_match = input_meta[
            (input_meta["Condition"] == condition) &
            (input_meta["Biorep"] == biorep)
        ]

        if len(input_match) != 1:
            raise RuntimeError(
                f"Expected exactly 1 input for "
                f"Condition={condition}, Biorep={biorep}; "
                f"found {len(input_match)}"
            )

        input_lib_id = input_match.iloc[0]["library.ID"]

        # ------------------------
        # Build paths
        # ------------------------
        ip_tagdir = tagdir_root / f"{lib_id}.{target}-tagdir"
        input_tagdir = tagdir_root / f"{input_lib_id}.{target}-tagdir"

        if not ip_tagdir.exists():
            raise FileNotFoundError(f"Missing IP tagdir: {ip_tagdir}")
        if not input_tagdir.exists():
            raise FileNotFoundError(f"Missing input tagdir: {input_tagdir}")

        output_file = peaks_outdir / f"{lib_id}.regions.txt"

        if output_file.exists() and not force_overwrite:
            print(f"Skipping {lib_id} — output exists")
            continue

        # ------------------------
        # Build HOMER command
        # ------------------------
        cmd = [
            "findPeaks",
            str(ip_tagdir),
            "-style", style,
            "-i", str(input_tagdir),
        ]

        if size is not None:
            cmd.extend(["-size", str(size)])
        if min_dist is not None:
            cmd.extend(["-minDist", str(min_dist)])

        print(f"Calling peaks for {lib_id}")
        print("CMD:", " ".join(cmd))

        with open(output_file, "w") as fh:
            subprocess.run(cmd, stdout=fh, check=True)
