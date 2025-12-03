#!/usr/bin/env python3
"""
Wrapper script to run the full ChIP-wrangler workflow from raw fastq to QC metrics.
"""

import argparse
import subprocess
import sys
from pathlib import Path

# Make sure the scripts directory is on the PYTHONPATH
SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.append(str(SCRIPT_DIR))

# Now import the modules by filename (without .py)
import preprocessing as pp
import trimming as trim
import alignment as align
import remove_duplicates as pcr
import generate_species_bams as species_bams
import get_sequencing_stats as stats
import estimate_spikein_ipeff as ipeff
import normalize_tagdirs as norm
import qc_data as qc


def wrangle_all(
    fastq_dir: str,
    genome_names: list,
    output_dir: str,
    paired_end: bool = False,
    threads: int = 4,
    skip_trimming: bool = False,
    skip_alignment: bool = False,
    indexed_genome_dir: str = None,
    fasta_paths: list = None
):
    """
    Runs the full ChIP-seq wrangling workflow with options to skip trimming or alignment.
    """

    fastq_dir = Path(fastq_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)

    # -------------------------------------------
    # STEP 0: Decide whether to build custom genome
    # -------------------------------------------

    if indexed_genome_dir:
        print("STEP 0: Skipping custom genome creation (using pre-indexed genome)")
        genome_dir = Path(indexed_genome_dir)

    else:
        print("STEP 0: Preprocessing (custom genome creation)")

        # New: pass ABSOLUTE FASTA FILES directly into your preprocessing function
        pp.create_custom_genome(
            genome_names=genome_names,
            output_dir=output_dir,
            fasta_paths=fasta_paths      # <–– NEW optional override
        )

        pp.generate_sample_metadata(fastq_dir=fastq_dir, output_dir=output_dir)

        genome_dir = output_dir / "custom_genome"

    # -------------------------------------------
    # STEP 1: Trimming
    # -------------------------------------------
    if not skip_trimming:
        print("STEP 1: Trimming")
        # STEP 1: Trimming
        trim.trim_fastq(
            fastq_dir=fastq_dir,          # raw fastq
            output_dir=output_dir,        # base project dir
            paired_end=paired_end
        )
    else:
        print("STEP 1: Trimming skipped!")

    # -------------------------------------------
    # STEP 2: Alignment
    # -------------------------------------------
    if not skip_alignment:
        print("STEP 2: Alignment")
        align.align_fastq(
            fastq_dir=output_dir/"fastq_trimmed",
            genome_dir=genome_dir,
            output_dir=output_dir,   # BASE PROJECT DIR
            paired_end=paired_end,
            threads=threads
        )
    else:
        print("STEP 2: Alignment skipped!")

    # -------------------------------------------
    # STEP 3: Remove PCR duplicates
    # -------------------------------------------
    print("STEP 3: Remove PCR duplicates")

    pcr.remove_duplicates(
        output_dir=output_dir,    # BASE PROJECT DIR
        paired=False,        # or however you detect paired-end
        umis=False,            # True/False
        threads=threads
    )
    
    # -------------------------------------------
    # STEP 4: Generate species-specific BAMs
    # -------------------------------------------
    
    print("STEP 4: Generate species-specific BAMs")
    from generate_species_bams import generate_species_bams

    # Assign target and spike species
    target_species = genome_names[0]           # e.g., "hg38"
    spike1_species = genome_names[1] if len(genome_names) > 1 else "none"
    spike2_species = genome_names[2] if len(genome_names) > 2 else "none"

    species_bams.generate_species_bams(
        spike1_species=spike1_species,
        spike2_species=spike2_species,
        target_species=target_species,
        user_dir=output_dir,       # this is the base project dir
        mapq_threshold=50,         # optional: change if needed
        threads=threads
    )

    # -------------------------------------------
    # STEP 5: Update sequencing stats
    # -------------------------------------------
    print("STEP 5: Update sequencing stats")
    
    stats.update_sample_metadata(
        user_dir=output_dir,
        target_species=genome_names[0],  # hg38
        spike1_species=genome_names[1] if len(genome_names) > 1 else None,  # dm6
        spike2_species=genome_names[2] if len(genome_names) > 2 else None,  # sac3
        samtools_path="samtools"
    )

    # -------------------------------------------
    # STEP 6: Estimate spike-in SNR
    # -------------------------------------------
    print("STEP 6: Estimate spike-in IP efficiency")

    # Define where to save normalized metadata
    save_file = output_dir / "sample_metadata.norm.tsv"

    ipeff.estimate_spikein(
        metadata_file=output_dir / "sample_metadata.tsv",  # input metadata
        user_dir=output_dir,                                # working directory
        spike1_species="dm6",                               # adjust if needed
        spike2_species="sac3",                              # adjust if needed
        target_species="hg38",                              # adjust if needed
        SNR_region="tss",                                   # optional
        frag_length=150,                                    # optional
        hist_size=4000,                                     # optional
        hist_bin=25,                                        # optional
        save_file=save_file                                 # always save
    )

    # -------------------------------------------
    # STEP 7: Normalize tag directories
    # -------------------------------------------

    print("STEP 7: Normalize tag directories")
    norm.normalize_tagdirs(
        metadata_file=output_dir / "sample_metadata.norm.tsv",
        genome_dirs=genome_dirs,
        target_species="hg38",
        spike_species=["dm6", "sac3"]
    )

    # -------------------------------------------
    # STEP 8: Generate QC data
    # -------------------------------------------
    print("STEP 8: Generate QC data")
    qc.generate_qc(
        metadata_file=output_dir/"sample_metadata.norm.tsv",
        output_dir=output_dir/"qc_outputs",
        target_species="hg38",
        spike1_species="dm6",
        spike2_species="sac3"
    )

    print("Workflow complete!")


def main():
    parser = argparse.ArgumentParser(description="Run full ChIP-seq wrangling workflow.")
    parser.add_argument("--fastq_dir", required=True)
    parser.add_argument("--genomes", nargs="+", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--paired_end", action="store_true")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--skip_trimming", action="store_true")
    parser.add_argument("--skip_alignment", action="store_true")

    # NEW ARGUMENTS:
    parser.add_argument(
        "--indexed_genome_dir",
        help="If provided, skip custom genome creation and use this directory with BWA index files."
    )
    parser.add_argument(
        "--fasta_paths",
        nargs="+",
        help="Absolute paths to FASTA files for custom combined genome construction."
    )

    args = parser.parse_args()

    wrangle_all(
        fastq_dir=args.fastq_dir,
        genome_names=args.genomes,
        output_dir=args.output_dir,
        paired_end=args.paired_end,
        threads=args.threads,
        skip_trimming=args.skip_trimming,
        skip_alignment=args.skip_alignment,
        indexed_genome_dir=args.indexed_genome_dir,
        fasta_paths=args.fasta_paths
    )


if __name__ == "__main__":
    main()
