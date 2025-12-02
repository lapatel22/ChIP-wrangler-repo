#!/usr/bin/env python3
"""
00_prepare_combined_genome.py

Flexible preprocessing script:
- Accepts absolute paths to target genome and up to 2 spike-in genomes
- OR accepts a pre-indexed genome directory and skips rebuilding
- If building, combines FASTAs, applies safe renaming to spike-in headers,
  outputs combined FASTA in <output_dir>, and BWA indexes it.
- Outputs the final indexed genome directory path for downstream script 02.
"""

import argparse
from pathlib import Path
import subprocess
import re
import sys


def safe_header(header: str, prefix: str) -> str:
    """Ensure unique, prefix-safe FASTA headers."""
    header = header.strip()
    if header.startswith(">"):
        header = header[1:]
    return f">{prefix}_{header}"


def rewrite_fasta_headers(input_fasta: Path, prefix: str) -> str:
    """Read FASTA and rewrite headers with prefix."""
    out = []
    with open(input_fasta) as f:
        for line in f:
            if line.startswith(">"):
                out.append(safe_header(line[1:], prefix) + "\n")
            else:
                out.append(line)
    return "".join(out)


def build_combined_genome(target: Path, spikes: list[Path], output_dir: Path) -> Path:
    """Combine FASTAs, apply header rules, return combined FASTA path."""
    output_dir.mkdir(parents=True, exist_ok=True)
    combined_fa = output_dir / "combined_genome.fa"

    print(f"→ Building combined FASTA: {combined_fa}")

    with open(combined_fa, "w") as out:
        # target fasta (unchanged)
        with open(target) as f:
            out.write(f.read())

        # spike-in fastas rewritten
        for i, spike in enumerate(spikes, start=1):
            prefix = f"spike{i}"
            print(f"  Adding spike-in: {spike} (prefix={prefix})")
            out.write(rewrite_fasta_headers(spike, prefix))

    return combined_fa


def bwa_index(fasta: Path):
    """Run BWA index on combined genome."""
    print(f"→ Running BWA index on: {fasta}")
    subprocess.run(["bwa", "index", str(fasta)], check=True)


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--target_fasta", type=Path, help="Absolute path to target genome FASTA")
    parser.add_argument("--spikein_fastas", type=Path, nargs="*", default=[],
                        help="Absolute paths to 0–2 spike-in FASTAs")
    parser.add_argument("--output_dir", type=Path, required=True,
                        help="Directory where the combined genome + index will be written")
    parser.add_argument("--indexed_genome_dir", type=Path,
                        help="If provided, skip building/indexing and use this directory instead")

    args = parser.parse_args()

    # -------------------------------------------------------------
    # CASE 1: User provided pre-indexed genome → SKIP building
    # -------------------------------------------------------------
    if args.indexed_genome_dir:
        idx = args.indexed_genome_dir.resolve()
        if not idx.exists():
            sys.exit(f"ERROR: Provided indexed genome directory does not exist: {idx}")

        print(f"✔ Skipping build — using pre-indexed genome at:\n  {idx}")
        print(f"OUTPUT_INDEX_DIR={idx}")
        return

    # -------------------------------------------------------------
    # CASE 2: Build and index combined genome
    # -------------------------------------------------------------
    if not args.target_fasta:
        sys.exit("ERROR: target_fasta is required unless --indexed_genome_dir is provided")

    all_spikes = args.spikein_fastas or []

    if len(all_spikes) > 2:
        sys.exit("ERROR: Maximum of 2 spike-in genomes supported.")

    # Build combined FASTA
    combined_fasta = build_combined_genome(
        args.target_fasta.resolve(),
        [p.resolve() for p in all_spikes],
        args.output_dir.resolve(),
    )

    # Run BWA index
    bwa_index(combined_fasta)

    # Output the final indexed genome directory
    print(f"OUTPUT_INDEX_DIR={args.output_dir.resolve()}")


if __name__ == "__main__":
    main()
