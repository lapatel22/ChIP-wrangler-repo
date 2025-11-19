#!/usr/bin/env python3
"""
Alignment module for ChIP-wrangler
- Detects paired-end or single-end FASTQs
- Aligns to custom genome using BWA MEM
- Creates logs
- Skips samples already aligned
"""

import subprocess
from pathlib import Path
import re
import argparse
from typing import List, Tuple, Dict


# ---------------------------------------------------------
# Helper functions
# ---------------------------------------------------------

def run(cmd: List[str], log_file: Path = None):
    """Run a command with optional logging."""
    print("Running:", " ".join(cmd))

    if log_file:
        with open(log_file, "w") as log:
            subprocess.run(cmd, stdout=log, stderr=log, check=True)
    else:
        subprocess.run(cmd, check=True)


def detect_fastq_pairs(fastq_dir: Path) -> List[Tuple[Path, Path]]:
    """
    Find paired-end FASTQ files.
    Supports patterns:
        sample_R1.fastq.gz / sample_R2.fastq.gz
        sample_1.fastq.gz  / sample_2.fastq.gz
        sample.1.fq.gz     / sample.2.fq.gz
    """
    r1_files = sorted(fastq_dir.glob("*R1*.fastq.gz")) + \
               sorted(fastq_dir.glob("*_1*.fastq.gz")) + \
               sorted(fastq_dir.glob("*.1.fq.gz"))

    pairs = []

    for r1 in r1_files:
        r2 = None

        # try to identify R2 file using common substitutions
        candidates = [
            Path(str(r1).replace("R1", "R2")),
            Path(str(r1).replace("_1", "_2")),
            Path(str(r1).replace(".1.trim.fastq.gz", ".2.trim.fastq.gz"))
        ]

        for c in candidates:
            if c.exists():
                r2 = c
                break

        if r2:
            pairs.append((r1, r2))

    return pairs


def detect_single_end_fastqs(fastq_dir: Path) -> List[Path]:
    """Return FASTQs that are not part of paired-end sets."""
    all_fastqs = set(fastq_dir.glob("*.fastq.gz"))
    pe_fastqs = set()

    # these are excluded from SE list
    for r1, r2 in detect_fastq_pairs(fastq_dir):
        pe_fastqs.add(r1)
        pe_fastqs.add(r2)

    return sorted(list(all_fastqs - pe_fastqs))


# ---------------------------------------------------------
# Main alignment function
# ---------------------------------------------------------

def align_fastq(
    fastq_dir: Path,
    genome_dir: Path,
    output_dir: Path,
    paired_end: bool = False,
    threads: int = 4
):
    """
    Align FASTQ files using BWA MEM.
    
    fastq_dir: directory containing trimmed FASTQs
    genome_dir: directory containing BWA index prefix ("<genome>_<spike1>_<spike2>")
    output_dir: directory for alignment output BAMs
    paired_end: user indicates that reads *should* be paired-end
    threads: CPU threads for BWA

    Skips samples where output BAM already exists.
    """

    output_dir.mkdir(parents=True, exist_ok=True)

    # detect BWA index prefix
    bwa_prefix = None
    for file in genome_dir.glob("*"):
        if file.suffix == ".amb":  # BWA index marker
            bwa_prefix = str(file.with_suffix(""))
            break

    if bwa_prefix is None:
        raise RuntimeError(f"Could not find BWA index prefix in: {genome_dir}")

    print(f"Using BWA index: {bwa_prefix}")

    # determine paired / single
    if paired_end:
        pairs = detect_fastq_pairs(fastq_dir)
        singles = []
    else:
        pairs = []
        singles = detect_single_end_fastqs(fastq_dir)

    # ----------- ALIGNMENT FOR PAIRED-END -----------
    for r1, r2 in pairs:
        sample = re.sub(r"_R1.*|_1.*|\.1\.fq\.gz", "", r1.stem)
        bam_out = output_dir / f"{sample}.bam"
        log_file = output_dir / f"{sample}.bwa.log"

        if bam_out.exists():
            print(f"Skipping {sample}: BAM already exists.")
            continue

        cmd = [
            "bwa", "mem",
            "-t", str(threads),
            bwa_prefix,
            str(r1), str(r2)
        ]

        # pipe into samtools to convert directly to BAM
        cmd_samtools = [
            "samtools", "view", "-bS", "-"
        ]

        print(f"Aligning paired-end sample: {sample}")

        with open(bam_out, "wb") as bam, open(log_file, "w") as log:
            p1 = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=log)
            p2 = subprocess.Popen(cmd_samtools, stdin=p1.stdout, stdout=bam, stderr=log)
            p1.stdout.close()
            p2.communicate()

    # ----------- ALIGNMENT FOR SINGLE-END -----------
    for fq in singles:
        sample = fq.stem
        bam_out = output_dir / f"{sample}.bam"
        log_file = output_dir / f"{sample}.bwa.log"

        if bam_out.exists():
            print(f"Skipping {sample}: BAM already exists.")
            continue

        cmd = [
            "bwa", "mem",
            "-t", str(threads),
            bwa_prefix,
            str(fq)
        ]

        print(f"Aligning single-end sample: {sample}")

        cmd_samtools = ["samtools", "view", "-bS", "-"]

        with open(bam_out, "wb") as bam, open(log_file, "w") as log:
            p1 = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=log)
            p2 = subprocess.Popen(cmd_samtools, stdin=p1.stdout, stdout=bam, stderr=log)
            p1.stdout.close()
