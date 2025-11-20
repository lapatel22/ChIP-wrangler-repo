#!/usr/bin/env python3
"""
trimming.py — Run fastp on FASTQ files and output trimmed reads.

Features:
- Auto-detect paired vs single-end FASTQs
- Uses fastp adapter autodetect
- Sliding-window trimming: 4bp window, mean Q < 20
- Minimum read length: 15bp
- Skips trimming if output FASTQs already exist
"""

import subprocess
from pathlib import Path
import argparse
import sys
from contextlib import contextmanager

class Tee(object):
    def __init__(self, *streams):
        self.streams = streams
    
    def write(self, data):
        for s in self.streams:
            s.write(data)
    
    def flush(self):
        for s in self.streams:
            s.flush()

@contextmanager
def tee_stdout(log_path):
    with open(log_path, "w") as logfile:
        tee = Tee(sys.stdout, logfile)
        old_stdout = sys.stdout
        sys.stdout = tee
        try:
            yield
        finally:
            sys.stdout = old_stdout

def detect_fastq_pairs(fastq_dir: Path):
    """
    Detect FASTQ files and group them into pairs or singles.

    Returns:
        list of dicts:
        [
          { "sample": "A", "R1": "A_R1.fastq.gz", "R2": "A_R2.fastq.gz" },
          { "sample": "B", "R1": "B.fastq.gz",     "R2": None }
        ]
    """
    fq_files = sorted(fastq_dir.glob("*.fastq.gz")) + sorted(fastq_dir.glob("*.fq.gz"))
    samples = {}

    for fq in fq_files:
        name = fq.name

        if "_R1" in name:
            key = name.replace("_R1", "")
            samples.setdefault(key, {"sample": key, "R1": None, "R2": None})
            samples[key]["R1"] = fq

        elif "_R2" in name:
            key = name.replace("_R2", "")
            samples.setdefault(key, {"sample": key, "R1": None, "R2": None})
            samples[key]["R2"] = fq

        else:
            # Single-end
            key = name.replace(".fastq.gz", "").replace(".fq.gz", "")
            samples[key] = {"sample": key, "R1": fq, "R2": None}

    return list(samples.values())


def trim_fastq(fastq_dir: Path, output_dir: Path, paired_end: bool = False, threads: int = 4):
    """
    Trim FASTQs using fastp.

    Args:
        fastq_dir: directory with raw FASTQs
        output_dir: directory for trimmed FASTQs
        paired_end: force paired-end? If False, auto-detect per-file
        threads: number of threads for fastp
    """

    fastq_dir = Path(fastq_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)

    sample_entries = detect_fastq_pairs(fastq_dir)

    for entry in sample_entries:
        sample = entry["sample"]
        R1 = entry["R1"]
        R2 = entry["R2"]

        # Override detection if paired_end=True
        is_pe = paired_end or (R1 is not None and R2 is not None)

        out_R1 = output_dir / f"{sample}_trimmed_R1.fastq.gz"
        out_R2 = output_dir / f"{sample}_trimmed_R2.fastq.gz" if is_pe else None

        json_report = output_dir / f"{sample}.fastp.json"
        html_report = output_dir / f"{sample}.fastp.html"

        # Skip if trimmed file already exists (R1 is enough to decide)
        if out_R1.exists():
            print(f"[Skipping] {sample} already trimmed → {out_R1}")
            continue

        print(f"Trimming sample: {sample} (paired-end={is_pe})")

        # Build fastp command
        cmd = [
            "fastp",
            "--thread", str(threads),
            "--trim_front1", "0",
            "--trim_tail1", "0",
            # Enable 3' end quality trimming
            "--cut_tail",
            # Quality-cut parameters
            "--cut_window_size", "4",
            "--cut_mean_quality", "20",
            "--length_required", "15",
            "--json", str(json_report),
            "--html", str(html_report),
            "--in1", str(R1),
            "--out1", str(out_R1)
        ]

        if is_pe:
            cmd += [
                "--in2", str(R2),
                "--out2", str(out_R2)
            ]

        # fastp auto-detects adapters by default → no adapter args required

        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            print(f"ERROR running fastp on sample {sample}: {e}")
            continue

    print("Trimming complete.")


# -------------------- CLI SUPPORT -------------------- #

def main():
    parser = argparse.ArgumentParser(description="Trim FASTQ files using fastp.")
    parser.add_argument("--fastq_dir", required=True, help="Input directory with raw FASTQs")
    parser.add_argument("--output_dir", required=True, help="Output directory for trimmed FASTQs")
    parser.add_argument("--paired_end", action="store_true", help="Force paired-end mode")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads for fastp")

    args = parser.parse_args()

    trim_fastq(
        fastq_dir=Path(args.fastq_dir),
        output_dir=Path(args.output_dir),
        paired_end=args.paired_end,
        threads=args.threads,
    )


if __name__ == "__main__":
    log_file = "trimming.log"
    with tee_stdout(log_file):
        main()
