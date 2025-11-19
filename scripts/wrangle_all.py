#!/usr/bin/env python3
"""
Wrapper script to run the full ChIP-wrangler workflow from raw fastq to QC metrics.
"""

### Example of how to run: 
# python run_wrangle_all.py \
#    --fastq_dir ./fastqfiles \
#    --genomes hg38 dm6 sacCer3 \
#    --output_dir ./chip_wrangled_outputs \
#    --paired_end \
#    --threads 8
####

import argparse
import subprocess
from pathlib import Path

# Import your functions
from your_package import (
    pre_processing as pp,           
    trimming as trim,               
    alignment as align,             
    remove_pcr_dups as pcr,        
    generate_species_bams as species_bams,  
    sequencing_stats as stats,      
    estimate_spikein_snr as snr,   
    normalize_tagdirs as norm,      
    qc_data as qc                   
)


def wrangle_all(
    fastq_dir: str,
    genome_names: list,
    output_dir: str,
    paired_end: bool = False,
    threads: int = 4,
    skip_trimming: bool = False,
    skip_alignment: bool = False
):
    """
    Runs the full ChIP-seq wrangling workflow with options to skip trimming or alignment.
    """

    fastq_dir = Path(fastq_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)

    print("STEP 0: Preprocessing")
    pp.create_custom_genome(genome_names=genome_names, output_dir=output_dir)
    pp.generate_sample_metadata(fastq_dir=fastq_dir, output_dir=output_dir)

    if not skip_trimming:
        print("STEP 1: Trimming")
        trim.trim_fastq(
            fastq_dir=fastq_dir,
            output_dir=output_dir/"fastq_trimmed",
            paired_end=paired_end
        )
    else:
        print("STEP 1: Trimming skipped!")

    if not skip_alignment:
        print("STEP 2: Alignment")
        align.align_fastq(
            fastq_dir=output_dir/"fastq_trimmed",
            genome_dir=output_dir/"custom_genome",
            output_dir=output_dir/"concat_align",
            paired_end=paired_end,
            threads=threads
        )
    else:
        print("STEP 2: Alignment skipped!")

    print("STEP 3: Remove PCR duplicates")
    pcr.remove_duplicates(
        input_dir=output_dir/"concat_align",
        output_dir=output_dir/"concat_align/dedup_out"
    )

    print("STEP 4: Generate species-specific BAMs")
    species_bams.split_species_bams(
        dedup_dir=output_dir/"concat_align/dedup_out",
        output_dirs={genome: output_dir/f"{genome}_data"/f"{genome}_aligned" for genome in genome_names}
    )

    print("STEP 5: Update sequencing stats")
    stats.update_sample_metadata(metadata_file=output_dir/"sample_metadata.tsv")

    print("STEP 6: Estimate spike-in SNR")
    snr.estimate_spikein(
        metadata_file=output_dir/"sample_metadata.tsv",
        genome_dirs={genome: output_dir/f"{genome}_data" for genome in genome_names}
    )

    print("STEP 7: Normalize tag directories")
    norm.normalize_tagdirs(
        metadata_file=output_dir/"sample_metadata.tsv",
        genome_dirs={genome: output_dir/f"{genome}_data" for genome in genome_names}
    )

    print("STEP 8: Generate QC data")
    qc.generate_qc(
        metadata_file=output_dir/"sample_metadata.tsv",
        output_dir=output_dir/"qc_outputs"
    )

    print("Workflow complete!")


def main():
    parser = argparse.ArgumentParser(description="Run full ChIP-seq wrangling workflow.")
    parser.add_argument("--fastq_dir", required=True, help="Directory with raw fastq files")
    parser.add_argument("--genomes", nargs="+", required=True, help="List of genomes, e.g. hg38 dm6 sacCer3")
    parser.add_argument("--output_dir", required=True, help="Directory to store outputs")
    parser.add_argument("--paired_end", action="store_true", help="Use if files are paired-end")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads for alignment/trimming")
    parser.add_argument("--skip_trimming", action="store_true", help="Skip the trimming step")
    parser.add_argument("--skip_alignment", action="store_true", help="Skip the alignment step")

    args = parser.parse_args()

    wrangle_all(
        fastq_dir=args.fastq_dir,
        genome_names=args.genomes,
        output_dir=args.output_dir,
        paired_end=args.paired_end,
        threads=args.threads,
        skip_trimming=args.skip_trimming,
        skip_alignment=args.skip_alignment
    )


if __name__ == "__main__":
    main()