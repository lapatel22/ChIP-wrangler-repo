#!/usr/bin/env python3
"""
Deduplication script for ChIP-wrangler outputs

Three modes:
1) Single-end, no UMIs
2) Paired-end, no UMIs
3) Any reads with UMIs
"""

import subprocess
from pathlib import Path
import argparse
import sys
from contextlib import contextmanager

# -------------------- Helper -------------------- #

class Tee:
    """Duplicate stdout/stderr to console + file."""
    def __init__(self, logfile):
        self.terminal = sys.stdout
        self.log = open(logfile, "w")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        self.terminal.flush()
        self.log.flush()

def run(cmd, log_file=None):
    """Run a shell command, optionally tee to log_file."""
    print(f"\nRunning: {' '.join(cmd)}\n")
    if log_file:
        tee = Tee(log_file)
        sys.stdout = tee
        sys.stderr = tee
    try:
        subprocess.run(cmd, check=True)
    finally:
        if log_file:
            sys.stdout = sys.__stdout__
            sys.stderr = sys.__stderr__

# -------------------- Deduplication -------------------- #

def dedup_single_end_no_umis(bam_dir, intermed_dir, dedup_dir, threads=4):
    bam_files = sorted(Path(bam_dir).glob("*.bam"))
    for bam in bam_files:
        base = bam.stem
        fixmate_bam = intermed_dir / f"{base}.fixmate.bam"
        sorted_bam = intermed_dir / f"{base}.sorted.bam"
        nodup_bam = dedup_dir / f"{base}.nodup.bam"
        log_file = dedup_dir / f"{base}.log"

        if nodup_bam.exists():
            print(f"Skipping {base}: output already exists.")
            continue

        run(["samtools", "fixmate", "-m", str(bam), str(fixmate_bam)], log_file)
        run(["samtools", "sort", str(fixmate_bam), "-o", str(sorted_bam)], log_file)
        run(["samtools", "markdup", "-r", "-s", str(sorted_bam), str(nodup_bam)], log_file)


def dedup_paired_end_no_umis(bam_dir, intermed_dir, dedup_dir, threads=4):
    bam_files = sorted(Path(bam_dir).glob("*.sam"))
    for sam in bam_files:
        base = sam.stem
        bam_file = intermed_dir / f"{base}.bam"
        collated_bam = intermed_dir / f"{base}.collated.bam"
        fixmate_bam = intermed_dir / f"{base}.fixmate.bam"
        sorted_bam = intermed_dir / f"{base}.sorted.bam"
        nodup_bam = dedup_dir / f"{base}.nodup.bam"
        log_file = dedup_dir / f"{base}.log"

        if nodup_bam.exists():
            print(f"Skipping {base}: output already exists.")
            continue

        run(["samtools", "view", "-bS", str(bam), "-o", str(bam_file)], log_file)
        run(["samtools", "collate", "-o", str(collated_bam), str(bam_file)], log_file)
        run(["samtools", "fixmate", "-m", str(collated_bam), str(fixmate_bam)], log_file)
        run(["samtools", "sort", str(fixmate_bam), "-o", str(sorted_bam)], log_file)
        run(["samtools", "markdup", "-r", "-s", str(sorted_bam), str(nodup_bam)], log_file)


def dedup_with_umis(intermed_dir, dedup_dir, threads=4):
    import glob

    sorted_bams = glob.glob(str(intermed_dir / "*.sorted.bam"))
    for bam in sorted_bams:
        base = Path(bam).stem.replace(".sorted", "")
        output_bam = dedup_dir / f"{base}.dedup.new.unique.bam"
        log_file = dedup_dir / f"{base}.log"

        if output_bam.exists():
            print(f"Skipping {base}: output already exists.")
            continue

        print(f"Processing {base} with umi_tools...")
        run([
            "umi_tools", "dedup",
            "--stdin", str(bam),
            "--stdout", str(output_bam),
            "--extract-umi-method", "read_id",
            "--umi-separator", ":",
            "--method", "unique",
            "--log", str(log_file),
            "--output-stats", str(dedup_dir / base)
        ], log_file)


# -------------------- CLI -------------------- #

def main():
    parser = argparse.ArgumentParser(description="Deduplicate aligned reads")
    parser.add_argument("--user_dir", type=Path, default=Path.cwd(),
                        help="Working directory (default: current working directory)")
    parser.add_argument("--paired", action="store_true", help="Indicate paired-end reads")
    parser.add_argument("--umis", required=True, choices=["TRUE","FALSE"], help="Whether UMIs are present in reads")
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()

    bam_dir = args.user_dir / "concat_align"
    intermed_dir = args.user_dir / "bams_pcrdup_intermed"
    dedup_dir = args.user_dir / "concat_align" / "dedup_out"

    dedup_dir.mkdir(exist_ok=True, parents=True)
    intermed_dir.mkdir(exist_ok=True, parents=True)

    if args.umis == "TRUE":
        dedup_with_umis(intermed_dir, dedup_dir, threads=args.threads)
    else:
        if args.paired:
            dedup_paired_end_no_umis(bam_dir, intermed_dir, dedup_dir, threads=args.threads)
        else:
            dedup_single_end_no_umis(bam_dir, intermed_dir, dedup_dir, threads=args.threads)


if __name__ == "__main__":
    main()
