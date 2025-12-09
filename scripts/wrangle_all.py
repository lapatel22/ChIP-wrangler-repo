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

# Import modules
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
    output_dir: str,
    target_genome: str,
    target_fasta: str,
    spike_genomes: list[str],
    spike_fastas: list[str],
    metadata: str,
    indexed_genome_dir: str = None,
    paired_end: bool = False,
    umis: bool = False,
    threads: int = 4,
    skip_trimming: bool = False,
    skip_alignment: bool = False,
    force_overwrite: bool = False
):
    """
    Runs the full ChIP-seq wrangling workflow with options to skip trimming or alignment.
    """

    ############ Define global variables used inside the functions #############

    fastq_dir = Path(fastq_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)

    target_fasta = Path(target_fasta)
    spike_fastas = [Path(sf) for sf in spike_fastas]

    target_species = target_genome
    spike1_species = spike_genomes[0] if len(spike_genomes) > 0 else None
    spike2_species = spike_genomes[1] if len(spike_genomes) > 1 else None

    # -------------------------------------------
    # STEP 0: Decide whether to build custom genome
    # -------------------------------------------

    if indexed_genome_dir:
        print("STEP 0: Skipping custom genome creation (using pre-indexed genome)")
        genome_dir = Path(indexed_genome_dir)
        genome_names = [target_genome] + spike_genomes
    else:
        print("STEP 0: Preprocessing (custom genome creation)")
        genome_dir, genome_names = pp.create_custom_genome(
            output_dir=output_dir,
            target_genome=target_genome,
            target_fasta=target_fasta,
            spike_genomes=spike_genomes,
            spike_fastas=spike_fastas
        )

    # -------------------------------------------
    # STEP 1: Trimming
    # -------------------------------------------

    trimmed_dir = output_dir / "fastq_trimmed"

    if skip_trimming and trimmed_dir.exists() and not force_overwrite:
        print("STEP 1: Trimming skipped!")
    else:
        print("STEP 1: Trimming")
        trim.trim_fastq(
            fastq_dir=fastq_dir,
            output_dir=output_dir,
            paired_end=paired_end,
            force_overwrite=force_overwrite
        )

    # -------------------------------------------
    # STEP 2: Alignment
    # -------------------------------------------

    alignment_dir = output_dir / "alignment"

    if skip_alignment and alignment_dir.exists() and not force_overwrite:
        print("STEP 2: Alignment skipped!")
    else:
        print("STEP 2: Alignment")
        align.align_fastq(
            fastq_dir=output_dir / "fastq_trimmed",
            genome_dir=genome_dir,
            output_dir=output_dir,
            paired_end=paired_end,
            threads=threads,
            force_overwrite=force_overwrite
        )

    # -------------------------------------------
    # STEP 3: Remove PCR duplicates
    # -------------------------------------------

    print("STEP 3: Remove PCR duplicates")

    pcr.remove_duplicates(
        output_dir=output_dir,
        paired=paired_end,
        umis=umis,
        threads=threads,
        force_overwrite=force_overwrite
    )

    # -------------------------------------------
    # STEP 4: Generate species-specific BAMs
    # -------------------------------------------

    print("STEP 4: Generate species-specific BAMs")

    species_bams.generate_species_bams(
        output_dir=output_dir,
        target_species=target_species,
        spike1_species=spike1_species,
        spike2_species=spike2_species,
        mapq_threshold=50,
        threads=threads,
        force_overwrite=force_overwrite
    )

    # -------------------------------------------
    # STEP 5: Update sequencing stats
    # -------------------------------------------

    print("STEP 5: Update sequencing stats")

    stats.update_sample_metadata(
        output_dir=output_dir,
        metadata=metadata,
        target_genome=target_genome,
        spike1_genome=spike1_species,
        spike2_genome=spike2_species,
        samtools_path="samtools"
    )

    # -------------------------------------------
    # STEP 6: Estimate spike-in IP efficiency
    # -------------------------------------------

    print("STEP 6: Estimate spike-in IP efficiency")

    save_file = output_dir / "sample_metadata.norm.tsv"

    ipeff.estimate_spikein(
        metadata_file=output_dir / "sample_metadata.tsv",
        output_dir=output_dir,
        spike1_species=spike1_species,
        spike2_species=spike2_species,
        target_species=target_species,
        SNR_region="tss",
        frag_length=150,
        hist_size=4000,
        hist_bin=25,
        save_file=save_file,
        force_overwrite=force_overwrite
    )

    # -------------------------------------------
    # STEP 7: Normalize tag directories
    # -------------------------------------------

    print("STEP 7: Normalize tag directories")

    genome_dirs = {genome: output_dir / f"{genome}_data" for genome in genome_names}

    norm.normalize_tagdirs(
        metadata_file=save_file,
        output_dir=output_dir,
        target_species=target_species,
        spike1_species=spike1_species,
        spike2_species=spike2_species,
        force_overwrite=force_overwrite
    )

    # -------------------------------------------
    # STEP 8: Generate QC data
    # -------------------------------------------

    print("STEP 8: Generate QC data")

    qc_dir = output_dir / "basedist_and_gc_outputs"
    qc_dir.mkdir(exist_ok=True, parents=True)
    log_file = qc_dir / "GC_bias_QC_report.log"

    with qc.tee_stdout(log_file):
        qc.generate_qc(
            metadata_file=save_file,
            output_dir=output_dir,
            target_species=target_species,
            spike1_species=spike1_species,
            spike2_species=spike2_species
        )

    print("Workflow complete!")


def main():
    parser = argparse.ArgumentParser(description="Run full ChIP-seq wrangling workflow.")
    parser.add_argument("--fastq_dir", required=True)
    parser.add_argument("--target_genome", required=True)
    parser.add_argument("--target_fasta", required=True)
    parser.add_argument("--spike_genomes", nargs="+", required=True)
    parser.add_argument("--spike_fastas", nargs="+", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--paired_end", action="store_true")
    parser.add_argument("--umis", action="store_true")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--skip_trimming", action="store_true")
    parser.add_argument("--skip_alignment", action="store_true")
    parser.add_argument("--force_overwrite", action="store_true")
    parser.add_argument("--indexed_genome_dir")

    args = parser.parse_args()

    wrangle_all(
        fastq_dir=args.fastq_dir,
        target_genome=args.target_genome,
        target_fasta=args.target_fasta,
        spike_genomes=args.spike_genomes,
        spike_fastas=args.spike_fastas,
        metadata=args.metadata,
        output_dir=args.output_dir,
        paired_end=args.paired_end,
        umis=args.umis,
        threads=args.threads,
        skip_trimming=args.skip_trimming,
        skip_alignment=args.skip_alignment,
        force_overwrite=args.force_overwrite,
        indexed_genome_dir=args.indexed_genome_dir
    )


if __name__ == "__main__":
    main()
