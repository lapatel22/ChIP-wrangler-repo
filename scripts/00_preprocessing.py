#!/usr/bin/env python3
"""
Flexible preprocessing script:
- Accepts any target genome and 0–2 spike-in genomes
- Finds genome FASTAs named <genome>_genome.fa inside ./genomes/<genome>/
- Applies safe suffix rules to spike-in FASTA headers
- Combines FASTAs
- BWA indexes combined genome
"""

import subprocess
from pathlib import Path
import re


# ======================== User-defined variables ========================

genome_dir = Path("../genomes")

target_genome = "hg38"
spikein_genomes = ["dm6", "sacCer3"]   # can be ["dm6"] or empty []


# ======================== Helper functions ========================

def run(cmd, shell=False):
    print(f"Running: {cmd}")
    subprocess.run(cmd, check=True, shell=shell)


def fasta_path(genome: str) -> Path:
    """
    Expected FASTA: ../genomes/<genome>_genome.fa
    """
    fp = genome_dir / f"{genome}_genome.fa"
    if not fp.exists():
        raise FileNotFoundError(f"ERROR: FASTA not found: {fp}")
    return fp


# ======================== Suffix Construction Rules ========================

def sanitize_suffix(name: str) -> str:
    """
    Keep only letters and numbers; remove everything else.
    """
    return re.sub(r'[^A-Za-z0-9]', '', name)


def strip_prefix(name: str) -> str:
    """
    If name contains a dot, remove the prefix (e.g., corn.AGPv3 → AGPv3).
    """
    if "." in name:
        return name.split(".", 1)[1]
    return name


def truncate(name: str, max_len=20) -> str:
    return name[:max_len]


def make_suffix(genome: str) -> str:
    """
    Return safe suffix label based on rules:
      - Custom short form for sacCer3 → sac3
      - Remove prefix before dot
      - Strip to alphanumeric only
      - Limit to 20 characters
    """
    special_map = {
        "sacCer3": "sac3"
    }

    if genome in special_map:
        return special_map[genome]

    core = strip_prefix(genome)
    core = sanitize_suffix(core)
    core = truncate(core)
    return core


# ======================== Step 1: Add spike-in suffixes ========================

spike_suffix_fasta_paths = []

for spike in spikein_genomes:

    orig_fa = fasta_path(spike)
    suffix_label = make_suffix(spike)
    suffix = f"_{suffix_label}"

    suffixed_fa = orig_fa.with_name(f"{spike}_genome{suffix}.fa")

    print(f"Adding suffix '{suffix}' to FASTA headers in {orig_fa}")

    cmd = f"sed 's/>.*/&{suffix}/' {orig_fa} > {suffixed_fa}"
    run(cmd, shell=True)

    print(f"Header check for {suffixed_fa}:")
    run(f"perl -ne 'if(/^>(\\S+)/){{print \"$1\\n\"}}' {suffixed_fa}", shell=True)

    spike_suffix_fasta_paths.append(suffixed_fa)


# ======================== Step 2: Combine FASTAs ========================

combined_name_parts = [target_genome] + spikein_genomes
combined_dir = genome_dir / "_".join(combined_name_parts)
combined_dir.mkdir(parents=True, exist_ok=True)

combined_fa = combined_dir / "genome_combined.fa"

files_to_cat = spike_suffix_fasta_paths + [fasta_path(target_genome)]

print(f"Combining genomes into {combined_fa}")
run(f"cat {' '.join(map(str, files_to_cat))} > {combined_fa}", shell=True)


# ======================== Step 3: BWA Index ========================

bwa_prefix = combined_dir / "_".join(combined_name_parts)

print(f"Indexing combined genome with BWA prefix: {bwa_prefix}")
run(["bwa", "index", "-p", str(bwa_prefix), str(combined_fa)])


print("\nAll steps completed successfully!")
