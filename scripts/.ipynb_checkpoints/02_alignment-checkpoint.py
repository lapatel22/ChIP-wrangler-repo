#!/usr/bin/env python3
"""
Alignment module for ChIP-wrangler
- Auto-builds directories from user_dir
- Detects paired-end or single-end FASTQs
- Aligns to user-specified genome using BWA MEM
- Creates logs
- Skips samples already aligned
"""

import subprocess
from pathlib import Path
import re
import argparse
from typing import List, Tuple

# -------------------- Helper functions --------------------

def clean_sample_name(fq_name: str) -> str:
    """Remove only _R1, _R2, and _trimmed suffixes from FASTQ file name."""
    name = fq_name
    name = re.sub(r"_R[12]", "", name)
    name = re.sub(r"_trimmed", "", name)
    name = re.sub(r"\.fastq$|\.fq$|\.fastq\.gz|\.fq\.gz", "", name)
    return name


def detect_fastq_pairs(fastq_dir: Path) -> List[Tuple[Path, Path]]:
    r1_files = sorted(fastq_dir.glob("*_R1*.fastq.gz"))
    pairs = []
    for r1 in r1_files:
        r2 = Path(str(r1).replace("_R1", "_R2"))
        if r2.exists():
            pairs.append((r1, r2))
    return pairs


def detect_single_end_fastqs(fastq_dir: Path) -> List[Path]:
    r1 = set(fastq_dir.glob("*_R1*.fastq.gz"))
    pairs_r1 = {p[0] for p in detect_fastq_pairs(fastq_dir)}
    singles = sorted(list(r1 - pairs_r1))
    return singles


# ============================================================
# ================  UPDATED FUNCTION HERE  ====================
# ============================================================

def align_fastq(user_dir: Path, genome_dir: Path, paired_end: bool = False, threads: int = 4):
    """Align FASTQ files using BWA MEM, using user_dir to define input/output dirs."""

    # Build directories automatically
    fastq_dir = user_dir / "fastq_trimmed"
    output_dir = user_dir / "concat_align"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Check for trimmed FASTQs
    trimmed_fastqs = list(fastq_dir.glob("*.fastq.gz"))
    if len(trimmed_fastqs) == 0:
        raise FileNotFoundError(
            f"No trimmed FASTQ files found in: {fastq_dir}\n"
            f"Expected files matching *.fastq.gz.\n"
            f"Did you run trimming?"
        )

    # Detect BWA index prefix inside genome_dir
    bwa_prefix = None
    for file in genome_dir.glob("*"):
        if file.suffix == ".amb":
            bwa_prefix = str(file.with_suffix(""))
            break
    if bwa_prefix is None:
        raise RuntimeError(f"Could not find BWA index prefix in: {genome_dir}")
    print(f"Using BWA index: {bwa_prefix}")

    # Detect paired/single
    if paired_end:
        pairs = detect_fastq_pairs(fastq_dir)
        singles = []
    else:
        pairs = []
        singles = detect_single_end_fastqs(fastq_dir)

    # Paired-end alignment
    for r1, r2 in pairs:
        sample = clean_sample_name(r1.stem)
        bam_out = output_dir / f"{sample}.bam"
        log_file = output_dir / f"{sample}.bwa.log"

        if bam_out.exists():
            print(f"Skipping {sample}: BAM already exists.")
            continue

        print(f"Aligning paired-end sample: {sample}")
        with open(bam_out, "wb") as bam, open(log_file, "w") as log:
            subprocess.run(
                f"bwa mem -t {threads} {bwa_prefix} {r1} {r2} | samtools view -bS -",
                shell=True,
                stdout=bam,
                stderr=log,
                check=True
            )

    # Single-end alignment
    for fq in singles:
        sample = clean_sample_name(fq.stem)
        bam_out = output_dir / f"{sample}.bam"
        log_file = output_dir / f"{sample}.bwa.log"

        if bam_out.exists():
            print(f"Skipping {sample}: BAM already exists.")
            continue

        print(f"Aligning single-end sample: {sample}")
        with open(bam_out, "wb") as bam, open(log_file, "w") as log:
            subprocess.run(
                f"bwa mem -t {threads} {bwa_prefix} {fq} | samtools view -bS -",
                shell=True,
                stdout=bam,
                stderr=log,
                check=True
            )


# -------------------- CLI --------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Align FASTQs with BWA MEM")

    # REQUIRED
    parser.add_argument("--genome_dir", required=True, type=Path,
                        help="Path to BWA index directory containing *.amb/*.ann/etc")

    # OPTIONAL user_dir
    parser.add_argument("--user_dir", type=Path, default=Path.cwd(),
                        help="Working directory (default: current working directory)")

    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--paired", action="store_true")

    args = parser.parse_args()

    # Run alignment with auto-constructed dirs
    align_fastq(
        user_dir=args.user_dir,
        genome_dir=args.genome_dir,
        paired_end=args.paired,
        threads=args.threads
    )
