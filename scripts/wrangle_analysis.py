#!/usr/bin/env python3

import argparse
import subprocess
from pathlib import Path

from peak_calling import run_homer_peak_calling
from create_counts_matrix import create_counts_matrix

### helper function to run the R script in the python wrapper #### 


### Note: adding a way to check where the R script is based on python scripts ###
# need to add this variable, which finds the python scripts, then strip out the actual python file, and use this as a variable when calling R script
import os

def get_Rscript_directory():

    script_path = os.path.realpath(__file__)
    script_dir = os.path.dirname(script_path)
    return script_dir

print(f"Script directory: {get_Rscript_directory()}")

#######

def run_deseq2_with_chipwrangler(
    r_script: Path = None,
    counts_tsv: Path = None,
    metadata_tsv: Path = None,
    conditions: str = "",
    outprefix: Path = None,
    spike_genomes: list[str] = []
):
    # Use your helper if r_script is not provided
    if r_script is None:
        r_script = Path(get_Rscript_directory()) / "deseq2_with_chipwrangler.R"
        
    """
    Run DESeq2 + ChIP-wrangler R script via Rscript.
    """

    cmd = [
        "Rscript",
        str(r_script),
        "--counts", str(counts_tsv),
        "--metadata", str(metadata_tsv),
        "--conditions", conditions,
        "--outprefix", str(outprefix)
    ] + ["--spike_genomes"] + spike_genomes

    print("Running DESeq2:")
    print("CMD:", " ".join(cmd))

    subprocess.run(cmd, check=True)


def main():

    parser = argparse.ArgumentParser(
        description="Unified ChIP-wrangler analysis pipeline, doing peak finding and differential peak analysis"
    )

    # General options
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--target", required=True)

    # Options for peak calling
    parser.add_argument("--style", required=True, choices=["histone", "factor"])
    parser.add_argument("--size", type=int)
    parser.add_argument("--minDist", type=int)

    parser.add_argument("--annot_size", type=str, default="given")

    # Options for DESeq2
    parser.add_argument("--conditions", required=True)
    parser.add_argument(
    "--r_script", default=Path(get_Rscript_directory()) /"deseq2_with_chipwrangler.R",)
    parser.add_argument("--spike_genomes", required=True, nargs="+",
    help="List of spike-in species (1 or 2), space-separated"
)

    parser.add_argument("--force_overwrite", action="store_true")

    args = parser.parse_args()
    output_dir = Path(args.output_dir).resolve()

    print("\nSTEP 1: Calling peaks with HOMER\n")
    run_homer_peak_calling(
        output_dir=output_dir,
        target=args.target,
        metadata_tsv=Path(args.metadata),
        style=args.style,
        size=args.size,
        min_dist=args.minDist,
        force_overwrite=args.force_overwrite,
    )
    
    import sys

    # ------------------------
    # STEP 2: Generating counts matrix
    # ------------------------
    
    data_dir = output_dir / f"{args.target}_data"
    condA, condB = args.conditions.split(",")
    
    default_counts_file = data_dir / f"raw_counts_from_merged_{condA}_{condB}.tsv"
    
    if default_counts_file.exists() and not args.force_overwrite:
        print(f" Counts matrix already exists: {default_counts_file}")
        counts_matrix = default_counts_file
    else:
        print("\nSTEP 2: Generating counts matrix\n")
        counts_matrix = create_counts_matrix(
            output_dir=output_dir,
            target=args.target,
            metadata_tsv=Path(args.metadata),
            conditions=args.conditions,
            annot_size=args.annot_size,
            out_tsv=default_counts_file,
        )


    print(f"DEBUG: Using counts matrix: {counts_matrix}")

    print("\nSTEP 3: Running DESeq2\n")
    run_deseq2_with_chipwrangler(
        r_script=Path(args.r_script),  # this can still override the default
        counts_tsv=counts_matrix,
        metadata_tsv=Path(args.metadata),
        conditions=args.conditions,
        outprefix=output_dir / f"{args.target}_DESeq2_results",
        spike_genomes=args.spike_genomes
    )


    print("\n Pipeline finished successfully\n")


if __name__ == "__main__":
    main()
