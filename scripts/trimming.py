#!/usr/bin/env python3
"""
trimming.py — Run fastp on FASTQ files and output trimmed reads.
"""

import subprocess
from pathlib import Path
import argparse
import sys
from contextlib import contextmanager

# --------------- Logging Helpers --------------- #

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

# --------------- FASTQ Detection --------------- #

def detect_fastq_pairs(fastq_dir: Path):
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
            key = name.replace(".fastq.gz", "").replace(".fq.gz", "")
            samples[key] = {"sample": key, "R1": fq, "R2": None}

    return list(samples.values())


# --------------- Main Logic --------------- #

def trim_fastq(
    fastq_dir: Path,
    output_dir: Path,
    paired_end: bool = False,
    threads: int = 4,
    phred_cutoff: int = 20
):
    """
    Trim FASTQs using fastp.
    fastq_dir: where raw FASTQ files live
    output_dir: where trimmed FASTQs will be written
    """
    fastq_dir = Path(fastq_dir)
    output_dir = Path(output_dir) / "fastq_trimmed"
    output_dir.mkdir(exist_ok=True, parents=True)

    sample_entries = detect_fastq_pairs(fastq_dir)

    for entry in sample_entries:
        sample = entry["sample"]
        R1 = entry["R1"]
        R2 = entry["R2"]

        is_pe = paired_end or (R1 is not None and R2 is not None)

        out_R1 = output_dir / f"{sample}_trimmed_R1.fastq.gz"
        out_R2 = output_dir / f"{sample}_trimmed_R2.fastq.gz" if is_pe else None

        json_report = output_dir / f"{sample}.fastp.json"
        html_report = output_dir / f"{sample}.fastp.html"

        if out_R1.exists():
            print(f"[Skipping] {sample} already trimmed → {out_R1}")
            continue

        print(f"Trimming sample: {sample} (paired-end={is_pe})")

        cmd = [
            "fastp",
            "--thread", str(threads),
            "--cut_tail",
            "--cut_window_size", "4",
            "--cut_mean_quality", str(phred_cutoff),
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

        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            print(f"ERROR running fastp on sample {sample}: {e}")
            continue

    print("Trimming complete.")


# --------------- CLI --------------- #

def main():
    parser = argparse.ArgumentParser(description="Trim FASTQ files using fastp.")
    parser.add_argument("--fastq_dir", required=True,
                        help="Directory containing raw FASTQ files")
    parser.add_argument("--output_dir", required=True,
                        help="Directory where trimmed FASTQs will be written")
    parser.add_argument("--paired_end", action="store_true")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--phred_cutoff", type=int, default=20)

    args = parser.parse_args()

    trim_fastq(
        fastq_dir=Path(args.fastq_dir),
        output_dir=Path(args.output_dir),
        paired_end=args.paired_end,
        threads=args.threads,
        phred_cutoff=args.phred_cutoff
    )


if __name__ == "__main__":
    with tee_stdout("trimming.log"):
        main()
