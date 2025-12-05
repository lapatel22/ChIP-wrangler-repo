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
    output_dir: str,
    target_genome: str,
    target_fasta: str,
    spike_genomes: str,
    spike_fastas: str,
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

    fastq_dir = Path(fastq_dir)
    target_fasta: Path(target_fasta)
    spike_fastas: Path(spike_fastas)
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
        genome_dir = pp.create_custom_genome(
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
    
    if skip_trimming and not force_overwrite and trimmed_dir.exists():
        print("STEP 1: Trimming skipped!")
    else:
        print("STEP 1: Trimming")
        # STEP 1: Trimming
        trim.trim_fastq(
            fastq_dir=fastq_dir,          # raw fastq
            output_dir=output_dir,        # base project dir
            paired_end=paired_end,
            force_overwrite=force_overwrite
        )
        
    # -------------------------------------------
    # STEP 2: Alignment
    # -------------------------------------------
    
    alignment_dir = output_dir / "alignment"

    if skip_alignment and not force_overwrite and alignment_dir.exists():
        print("STEP 2: Alignment skipped!")
    else:
        print("STEP 2: Alignment")
        align.align_fastq(
            fastq_dir=output_dir/"fastq_trimmed",
            genome_dir=genome_dir,
            output_dir=output_dir,   # BASE PROJECT DIR
            paired_end=paired_end,
            threads=threads,
            force_overwrite=force_overwrite
        )


    # -------------------------------------------
    # STEP 3: Remove PCR duplicates
    # -------------------------------------------
    print("STEP 3: Remove PCR duplicates")

    pcr.remove_duplicates(
        output_dir=output_dir,    # BASE PROJECT DIR
        paired=False,        # or however you detect paired-end
        umis=False,            # True/False
        threads=threads,
        force_overwrite=force_overwrite
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
        threads=threads,
        force_overwrite=force_overwrite
    )

    # -------------------------------------------
    # STEP 5: Update sequencing stats
    # -------------------------------------------
    print("STEP 5: Update sequencing stats")
    
    stats.update_sample_metadata(
        user_dir=output_dir,
        metadata=metadata,
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
        save_file=save_file,                                # always save
        force_overwrite=force_overwrite
    )

    # -------------------------------------------
    # STEP 7: Normalize tag directories
    # -------------------------------------------

    print("STEP 7: Normalize tag directories")
    norm.normalize_tagdirs(
        metadata_file=output_dir/"sample_metadata.norm.tsv",
        genome_dirs={genome: output_dir/f"{genome}_data" for genome in genome_names}, force_overwrite=force_overwrite
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
            metadata_file=output_dir / "sample_metadata.norm.tsv",
            output_dir=output_dir,  
            target_species="hg38",
            spike1_species="dm6",
            spike2_species="sac3"
        )
    
    print("Workflow complete!")



def main():
    parser = argparse.ArgumentParser(description="Run full ChIP-seq wrangling workflow.")
    parser.add_argument("--fastq_dir", required=True)
    parser.add_argument("--target_genome", required=True, 
                       help = "Name of target genome")
    parser.add_argument("--target_fasta", required=True,
                        help="Absolute paths to FASTA files for custom combined genome construction.")
    parser.add_argument("--spike_genomes", nargs="+", required=True, 
                       help = "Name of spike-in genomes (list of strings, example: dm6 sac3)")
    parser.add_argument("--spike_fastas", nargs="+", required=True,
                        help="Absolute paths to spike-in FASTA files for custom combined genome construction, example: dm6_genome.fa sacCer3_genome.fa")
    parser.add_argument("--metadata", required=True, 
                        help = "Path to the sample metadata tsv file")
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--paired_end", action="store_true")
    parser.add_argument("--umis", action="store_true")
    parser.add_argument("--threads", type=int, default=1, help = "Specifies number of cpus to use when multithreding is possible")
    parser.add_argument("--skip_trimming", action="store_true")
    parser.add_argument("--skip_alignment", action="store_true")
    parser.add_argument("--force_overwrite", action="store_true",
    help="Global option to force overwriting of every step."
)
    parser.add_argument("--indexed_genome_dir",
                        help="If provided, skip custom genome creation and use this directory containing BWA index files.")

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
