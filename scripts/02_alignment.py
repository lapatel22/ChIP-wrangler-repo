#!/usr/bin/env python3
"""
Alignment module for ChIP-wrangler
- Auto-builds directories from user_dir, target, and spike-ins
- Detects paired-end or single-end FASTQs
- Aligns to custom genome using BWA MEM
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
    name = re.sub(r"_R[12]", "", name)       # remove _R1 or _R2 at the end
    name = re.sub(r"_trimmed", "", name)     # remove _trimmed at the end
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


def align_fastq(fastq_dir: Path, genome_dir: Path, output_dir: Path, paired_end: bool = False, threads: int = 4):
    """Align FASTQ files using BWA MEM"""
    output_dir.mkdir(parents=True, exist_ok=True)

    # Detect BWA index prefix
    bwa_prefix = None
    for file in genome_dir.glob("*"):
        if file.suffix == ".amb":
            bwa_prefix = str(file.with_suffix(""))
            break
    if bwa_prefix is None:
        raise RuntimeError(f"Could not find BWA index prefix in: {genome_dir}")
    print(f"Using BWA index: {bwa_prefix}")

    # Determine paired / single
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
    parser.add_argument("--user_dir", required=True, type=Path, help="Base directory for project")
    parser.add_argument("--target_genome", required=True, type=str, help="Target genome, e.g., hg38")
    parser.add_argument("--spikein_genomes", required=True, nargs="+", help="Space-separated spike-ins, e.g., dm6 sac3")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--paired", action="store_true")
    args = parser.parse_args()

    # Convert sac3 → sacCer3 in genome dir name
    spikein_clean = [g if g != "sac3" else "sacCer3" for g in args.spikein_genomes]

    # Build directories
    fastq_dir = args.user_dir / "fastq_trimmed"
    genome_name = "_".join([args.target_genome] + spikein_clean)
    genome_dir = args.user_dir / "genomes" / genome_name
    output_dir = args.user_dir / "concat_align"

    align_fastq(
        fastq_dir=fastq_dir,
        genome_dir=genome_dir,
        output_dir=output_dir,
        paired_end=args.paired,
        threads=args.threads
    )
