#!/usr/bin/env python3
"""
preprocessing.py

Provides:
- create_custom_genome(): Python function for wrapper use
- CLI usage preserving old behavior
"""

import argparse
from pathlib import Path
import subprocess
import sys

# ------------------------- Helper functions -------------------------

def safe_header(header: str, suffix: str) -> str:
    """Append a suffix to FASTA headers (e.g. chr2L_dm6)."""
    header = header.strip()
    if header.startswith(">"):
        header = header[1:]
    return f">{header}_{suffix}"


def rewrite_fasta_headers(input_fasta: Path, suffix: str) -> str:
    """Read FASTA and rewrite headers with a genome suffix."""
    out = []
    with open(input_fasta) as f:
        for line in f:
            if line.startswith(">"):
                out.append(safe_header(line[1:], suffix) + "\n")
            else:
                out.append(line)
    return "".join(out)


def build_combined_genome(target: Path, spikes: list[Path], spike_genomes: list[str], output_dir: Path) -> Path:

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
            suffix = spike_genomes[i-1]   # use actual genome name
            print(f"  Adding spike-in: {spike} (suffix={suffix})")
            out.write(rewrite_fasta_headers(spike, suffix))

    return combined_fa

def bwa_index(fasta: Path):
    """Run BWA index on combined genome."""
    print(f"→ Running BWA index on: {fasta}")
    subprocess.run(["bwa", "index", str(fasta)], check=True)

# ------------------------- Wrapper function -------------------------

def create_custom_genome(output_dir, target_genome, target_fasta, spike_genomes=[], spike_fastas=[]):
    """
    Build a combined genome from target + spike-ins and BWA index it.
    The indexed genome will be placed in a directory named like:
        <target>_<spike1>_<spike2>_indexed
    inside output_dir.
    """
    output_dir = Path(output_dir)

    # -------------------- build indexed genome dir name --------------------
    parts = [target_genome] + spike_genomes
    index_dir_name = "_".join(parts) + "_indexed"
    index_dir = output_dir / index_dir_name
    index_dir.mkdir(parents=True, exist_ok=True)

    # -------------------- resolve spike FASTAs --------------------
    if spike_fastas:
        if len(spike_fastas) != len(spike_genomes):
            raise ValueError("Number of spike_fastas must match number of spike_genomes")
        spikes_resolved = [Path(f).resolve() for f in spike_fastas]
    else:
        spikes_resolved = []

    # -------------------- build combined FASTA --------------------
    combined_fasta = build_combined_genome(
        target=Path(target_fasta).resolve(),
        spikes=spikes_resolved,
        spike_genomes=spike_genomes,
        output_dir=index_dir
    )

    # -------------------- BWA index --------------------
    bwa_index(combined_fasta)

    return index_dir


# ------------------------- CLI usage -------------------------

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--target_genome", type=str, required=True,
                        help="Nickname of the target genome, e.g., 'hg38'")
    parser.add_argument("--target_fasta", type=Path, required=True,
                        help="Absolute path to target genome FASTA")
    parser.add_argument("--spike_genomes", type=str, nargs="*", default=[],
                        help="0–2 spike-in genome nicknames, e.g., dm6 sacCer3")
    parser.add_argument("--spike_fastas", type=Path, nargs="*", default=[],
                        help="Absolute paths to 0–2 spike-in FASTAs, must match spike_genomes order")
    parser.add_argument("--output_dir", type=Path, required=True,
                        help="Directory where the combined genome + index will be written")
    parser.add_argument("--indexed_genome_dir", type=Path,
                        help="If provided, skip building/indexing and use this directory instead")

    args = parser.parse_args()

    # Use pre-indexed genome if provided
    if args.indexed_genome_dir:
        idx = args.indexed_genome_dir.resolve()
        if not idx.exists():
            sys.exit(f"ERROR: Provided indexed genome directory does not exist: {idx}")
        print(f" Skipping build — using pre-indexed genome at:\n  {idx}")
        print(f"OUTPUT_INDEX_DIR={idx}")
        return

    # Build combined genome
    create_custom_genome(
        output_dir=args.output_dir.resolve(),
        target_genome=args.target_genome,
        target_fasta=args.target_fasta.resolve(),
        spike_genomes=args.spike_genomes,
        spike_fastas=[p.resolve() for p in args.spike_fastas]
    )

    # Print final output
    print(f"OUTPUT_INDEX_DIR={args.output_dir.resolve()}")


if __name__ == "__main__":
    main()
