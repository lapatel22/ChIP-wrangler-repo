#!/usr/bin/env python3
"""
Flexible preprocessing script:
- Accepts target genome and 1–2 spike-in genomes
- Finds genome FASTAs named <genome>_genome.fa inside <user_dir>/genomes/
- Applies safe suffix rules to spike-in FASTA headers
- Combines FASTAs
- BWA indexes combined genome
"""

import subprocess
from pathlib import Path
import re
import argparse


# ======================== Helper functions ========================

def run(cmd, shell=False):
    print(f"Running: {cmd}")
    subprocess.run(cmd, check=True, shell=shell)


def fasta_path(genome: str, genome_dir: Path) -> Path:
    """
    Expected FASTA: <genome_dir>/<genome>_genome.fa
    """
    fp = genome_dir / f"{genome}_genome.fa"
    if not fp.exists():
        raise FileNotFoundError(f"ERROR: FASTA not found: {fp}")
    return fp


# ======================== Suffix Construction Rules ========================

def sanitize_suffix(name: str) -> str:
    return re.sub(r'[^A-Za-z0-9]', '', name)


def strip_prefix(name: str) -> str:
    if "." in name:
        return name.split(".", 1)[1]
    return name


def truncate(name: str, max_len=20) -> str:
    return name[:max_len]


def make_suffix(genome: str) -> str:
    special_map = {"sacCer3": "sac3"}
    if genome in special_map:
        return special_map[genome]

    core = strip_prefix(genome)
    core = sanitize_suffix(core)
    core = truncate(core)
    return core


# ======================== Main ========================

def main():
    parser = argparse.ArgumentParser(description="Combine genomes and create BWA index with safe suffixes.")
    parser.add_argument("--user_dir", type=Path, required=True,
                        help="Base directory containing the 'genomes/' folder")
    parser.add_argument("--target_genome", required=True, help="Primary genome name (e.g., hg38)")
    parser.add_argument("--spikein_genomes", nargs='+', required=True,
                        help="Required spike-in genome(s) (1 or more)")

    args = parser.parse_args()

    genome_dir = args.user_dir / "genomes"
    target_genome = args.target_genome
    spikein_genomes = args.spikein_genomes

    # ======================== Step 1: Add spike-in suffixes ========================

    spike_suffix_fasta_paths = []

    for spike in spikein_genomes:
        orig_fa = fasta_path(spike, genome_dir)
        suffix_label = make_suffix(spike)
        suffix = f"_{suffix_label}"

        suffixed_fa = orig_fa.with_name(f"{spike}_genome{suffix}.fa")

        print(f"Adding suffix '{suffix}' to FASTA headers in {orig_fa}")
        cmd = f"sed 's/>.*/&{suffix}/' {orig_fa} > {suffixed_fa}"
        run(cmd, shell=True)

        spike_suffix_fasta_paths.append(suffixed_fa)

    # ======================== Step 2: Combine FASTAs ========================

    combined_name_parts = [target_genome] + spikein_genomes
    combined_dir = genome_dir / "_".join(combined_name_parts)
    combined_dir.mkdir(parents=True, exist_ok=True)

    combined_fa = combined_dir / "genome_combined.fa"

    files_to_cat = spike_suffix_fasta_paths + [fasta_path(target_genome, genome_dir)]

    print(f"Combining genomes into {combined_fa}")
    run(f"cat {' '.join(map(str, files_to_cat))} > {combined_fa}", shell=True)

    # ======================== Step 3: BWA Index ========================

    bwa_prefix = combined_dir / "_".join(combined_name_parts)

    print(f"Indexing combined genome with BWA prefix: {bwa_prefix}")
    run(["bwa", "index", "-p", str(bwa_prefix), str(combined_fa)])

    print("\nAll steps completed successfully!")


if __name__ == "__main__":
    main()
