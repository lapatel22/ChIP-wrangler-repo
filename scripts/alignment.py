#!/usr/bin/env python3
"""
Alignment module for ChIP-wrangler
"""

import subprocess
from pathlib import Path
import re
import argparse
from typing import List, Tuple

# -------------------- Helper functions --------------------

def clean_sample_name(fq_name: str) -> str:
    name = fq_name
    name = re.sub(r"_R[12]", "", name)
    name = re.sub(r"_trimmed", "", name)
    name = re.sub(r"\.fastq$|\.fq$|\.fastq\.gz|\.fq\.gz", "", name)
    return name


def detect_fastq_pairs(fastq_dir: Path):
    r1_files = sorted(fastq_dir.glob("*_R1*.fastq.gz"))
    pairs = []
    for r1 in r1_files:
        r2 = Path(str(r1).replace("_R1", "_R2"))
        if r2.exists():
            pairs.append((r1, r2))
    return pairs


def detect_single_end_fastqs(fastq_dir: Path):
    r1 = set(fastq_dir.glob("*_R1*.fastq.gz"))
    pairs_r1 = {p[0] for p in detect_fastq_pairs(fastq_dir)}
    singles = sorted(list(r1 - pairs_r1))
    return singles


# ============================================================
# ============ UPDATED FUNCTION (WRAPPER-FRIENDLY) ============
# ============================================================

def align_fastq(
    fastq_dir: Path = None,
    target_genome: str = None,
    spike_genomes: list[str] = [],
    genome_dir: Path = None,
    output_dir: Path = None,
    paired_end: bool = False,
    threads: int = 1,
    force_overwrite: bool = False
):
    """
    Align FASTQ files using BWA MEM.

    Parameters
    ----------
    genome_dir : Path
        Directory containing BWA index files (.amb, .ann, etc.)
    fastq_dir : Path
        Directory containing trimmed FASTQ files.
        Default: <output_dir>/fastq_trimmed
    output_dir : Path
        Directory where BAMs will be written.
        Default: <output_dir>/concat_align
    """

    output_dir = Path(output_dir)

    # ------------------------------
    # Determine genome_dir automatically if not provided
    # ------------------------------
    if genome_dir is None:
        if target_genome is None:
            raise ValueError("target_genome must be provided if genome_dir is not set")
        parts = [target_genome] + spike_genomes
        genome_dir = output_dir / ("_".join(parts) + "_indexed")

    genome_dir = Path(genome_dir)
    if not genome_dir.exists():
        raise FileNotFoundError(f"Genome directory does not exist: {genome_dir}")

    print(f"Using genome directory: {genome_dir}")

    # ------------------------------
    # Resolve fastq_dir and bam_dir
    # ------------------------------
    if fastq_dir is None:
        fastq_dir = output_dir / "fastq_trimmed"
    else:
        fastq_dir = Path(fastq_dir)

    bam_dir = output_dir / "concat_align"
    bam_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------
    # Check for trimmed FASTQs
    # ------------------------------
    trimmed_fastqs = list(fastq_dir.glob("*.fastq.gz"))
    if len(trimmed_fastqs) == 0:
        raise FileNotFoundError(
            f"No trimmed FASTQ files found in {fastq_dir}.\n"
            f"Did you run trimming?"
        )

    # ------------------------------
    # Find the BWA index prefix
    # ------------------------------

    bwa_prefix = None
    for file in genome_dir.glob("*"):
        if file.suffix == ".amb":
            bwa_prefix = file.with_suffix("")  # strip .amb
            break

    if bwa_prefix is None:
        raise RuntimeError(f"No BWA index (.amb) found in {genome_dir}")

    print(f"Using BWA index prefix: {bwa_prefix}")

    # ------------------------------
    # Pair detection
    # ------------------------------
    if paired_end:
        pairs = detect_fastq_pairs(fastq_dir)
        singles = []
    else:
        pairs = []
        singles = detect_single_end_fastqs(fastq_dir)

    # ------------------------------
    # Paired-end alignment
    # ------------------------------
    for r1, r2 in pairs:
        sample = clean_sample_name(r1.stem)
        bam_out = bam_dir / f"{sample}.bam"
        log_file = bam_dir / f"{sample}.bwa.log"

        if bam_out.exists() and not force_overwrite:
            print(f"[Skipping] BAM exists for {sample}")
            continue

        print(f"Aligning paired-end: {sample}")

        with open(bam_out, "wb") as bam, open(log_file, "w") as log:
            subprocess.run(
                f"bwa mem -t {threads} {bwa_prefix} {r1} {r2} | samtools view -bS -",
                shell=True,
                stdout=bam,
                stderr=log,
                check=True
            )

    # ------------------------------
    # Single-end alignment
    # ------------------------------
    for fq in singles:
        sample = clean_sample_name(fq.stem)
        bam_out = bam_dir / f"{sample}.bam"
        log_file = bam_dir / f"{sample}.bwa.log"

        if bam_out.exists() and not force_overwrite:
            print(f"[Skipping] BAM exists for {sample}")
            continue

        print(f"Aligning single-end: {sample}")

        with open(bam_out, "wb") as bam, open(log_file, "w") as log:
            subprocess.run(
                f"bwa mem -t {threads} {bwa_prefix} {fq} | samtools view -bS -",
                shell=True,
                stdout=bam,
                stderr=log,
                check=True
            )

    print("Alignment completed.")


# -------------------- CLI --------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Align FASTQs with BWA MEM")
    parser.add_argument("--target_genome", type=str, help="Nickname of target genome, e.g., hg38")
    parser.add_argument("--spike_genomes", nargs="*", default=[], help="Spike-in genome nicknames")
    parser.add_argument("--genome_dir", required=True, type=Path)
    parser.add_argument("--output_dir", required=True, type=Path,
                        help="Base directory containing fastq_trimmed/ and where concat_align/ will be written")
    parser.add_argument("--fastq_dir", type=Path,
                        help="Directory containing trimmed FASTQs; default <output_dir>/fastq_trimmed")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--paired", action="store_true")
    parser.add_argument("force_overwrite", action="store_true", 
                       help = "Global option to force overwriting of every step")

    args = parser.parse_args()

    align_fastq(
        genome_dir=args.genome_dir,
        target_genome=args.target_genome,
        spike_genomes=args.spike_genomes,
        fastq_dir=args.fastq_dir,
        output_dir=args.output_dir,
        paired_end=args.paired,
        threads=args.threads,
        force_overwrite=args.force_overwrite
    )
