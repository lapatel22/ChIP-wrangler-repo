import os
import glob
import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

# ======================== Configurable Variables ========================

input_dir = Path("../concat_align/dedup_out") ### Make user input ARGUMENT!
output_dir = Path("../bams_chr_sep")
filtered_dir = Path("../filtered_bams")

spike_species_1 = "dm6"
spike_species_2 = "sac3"
target_species = "hg38"

bamutil_path = "bam"  # path to bamUtil binary
samtools_path = "samtools"              # can also set to full path if needed

suffix_to_remove = ".concat.nodup.bam"
mapq_threshold = 50                         ### Make user input!
max_threads = 4                           ### Make user input!

# Create directories if needed
output_dir.mkdir(parents=True, exist_ok=True)
filtered_dir.mkdir(parents=True, exist_ok=True)

# ======================== Step 1: Filter low MAPQ reads ========================

print("=== Step 1: Filtering low-quality reads ===")

bam_files = glob.glob(str(input_dir / f"*{suffix_to_remove}"))

def filter_bam(bam_file, min_mapq):
    bam_file = Path(bam_file)
    base = bam_file.name.replace(suffix_to_remove, "")
    filtered_bam = filtered_dir / f"{base}.filtered.bam"
    
    if filtered_bam.exists():
        return f"Skipping {filtered_bam.name} (already exists)"
    
    cmd = [
        samtools_path, "view",
        "-b", "-q", str(min_mapq), "-F", "4",
        "-o", str(filtered_bam),
        str(bam_file)
    ]
    subprocess.run(cmd, check=True)
    return f"Filtered {bam_file.name} → {filtered_bam.name}"

with ThreadPoolExecutor(max_workers=max_threads) as executor:
    futures = [executor.submit(filter_bam, f, mapq_threshold) for f in bam_files]
    for fut in as_completed(futures):
        print(fut.result())

# ======================== Step 2: Split BAMs by chromosome ========================

print("\n=== Step 2: Splitting BAMs by chromosome ===")

filtered_bams = list(filtered_dir.glob("*.filtered.bam"))

def split_bam(input_bam):
    input_bam = Path(input_bam)
    base = input_bam.stem.replace(".filtered", "")
    output_prefix = output_dir / f"{base}.dedup."

    sentinel_chrX = output_dir / f"{base}.dedup.chrX.bam"
    if sentinel_chrX.exists():
        return f"Skipping {base} (already split)"
    
    cmd = [
        bamutil_path,
        "splitChromosome",
        "--in", str(input_bam),
        "--out", str(output_prefix)
    ]
    subprocess.run(cmd, check=True)
    return f"Finished splitting {base}"

with ThreadPoolExecutor(max_workers=max_threads) as executor:
    futures = [executor.submit(split_bam, bam) for bam in filtered_bams]
    for fut in as_completed(futures):
        print(fut.result())

# ======================== Step 2b: Remove unwanted chromosomes ========================

print("\n=== Step 2b: Removing unwanted contigs (alt, Un, random) ===")

for bam in filtered_dir.glob("*.filtered.bam"):
    base = bam.stem.replace(".filtered", "")
    for bad_bam in output_dir.glob(f"{base}.dedup.chr*"):
        # Match unwanted contigs
        if any(x in bad_bam.name for x in ["_alt", "Un", "random"]):
            print(f"Removing {bad_bam.name} (unwanted contig)")
            bad_bam.unlink()

# ======================== Step 3: Merge by species ========================

print("\n=== Step 3: Merging chromosomes by species ===")

# Create output directories for each species
species_map = {
    spike_species_1: Path(f"../{spike_species_1}_data/{spike_species_1}_aligned"),
    spike_species_2: Path(f"../{spike_species_2}_data/{spike_species_2}_aligned"),
    target_species: Path(f"../{target_species}_data/{target_species}_aligned")
}

for path in species_map.values():
    path.mkdir(parents=True, exist_ok=True)

# Folder to hold spike-in chromosome BAMs
spike_chr_dir = output_dir / "spike_bams_chr_sep"
spike_chr_dir.mkdir(exist_ok=True)

# Loop through all filtered BAMs
for bam in filtered_bams:
    base = Path(bam).stem.replace(".filtered", "")

    # Move spike-in chromosome BAMs to spike_chr_dir
    for pattern in [f"{base}.dedup.chr*_{spike_species_1}.bam",
                    f"{base}.dedup.chr*_{spike_species_2}.bam"]:
        for f in output_dir.glob(pattern):
            f.rename(spike_chr_dir / f.name)

    # Define output BAM paths
    outputs = {
        spike_species_1: species_map[spike_species_1] / f"{base}.{spike_species_1}.bam",
        spike_species_2: species_map[spike_species_2] / f"{base}.{spike_species_2}.bam",
        target_species: species_map[target_species] / f"{base}.{target_species}.bam",
    }

    # Merge spike-in species
    for species in [spike_species_1, spike_species_2]:
        spike_files = sorted(spike_chr_dir.glob(f"{base}.dedup.chr*_{species}.bam"))
        out_bam = outputs[species]
        if out_bam.exists():
            print(f"Skipping {out_bam.name} (already exists)")
            continue
        if not spike_files:
            print(f"No chromosome files found for {base} / {species}")
            continue

        print(f"Merging {species} BAMs for {base}")
        subprocess.run([samtools_path, "merge", str(out_bam)] + [str(f) for f in spike_files], check=True)

    # Merge target species
    target_files = sorted(output_dir.glob(f"{base}.dedup.chr*.bam"))
    out_bam = outputs[target_species]
    if out_bam.exists():
        print(f"Skipping {out_bam.name} (already exists)")
        continue
    if not target_files:
        print(f"No chromosome files found for {base} / {target_species}")
        continue

    print(f"Merging {target_species} BAMs for {base}")
    subprocess.run([samtools_path, "merge", str(out_bam)] + [str(f) for f in target_files], check=True)

# ======================== Step 4: Remove species suffixes (post-merge) ========================

print("\n=== Step 4: Removing species suffixes from merged SAMs ===")

def remove_suffixes(species, dir_path):
    for bam in dir_path.glob(f"*.{species}.bam"):
        base = bam.stem
        nosuff_sam = dir_path / f"{base}.nosuff.sam"
        nosuff2_sam = dir_path / f"{base}.nosuffx2.sam"

        if nosuff2_sam.exists():
            print(f"Skipping {nosuff2_sam.name} (already done)")
            continue

        # Convert BAM to SAM
        subprocess.run([samtools_path, "view", "-h", "-o", str(nosuff_sam), str(bam)], check=True)

        # Remove both suffixes in order
        with open(nosuff_sam) as fin, open(nosuff2_sam, "w") as fout:
            for line in fin:
                line = line.replace(f"_{spike_species_1}", "").replace(f"_{spike_species_2}", "")
                fout.write(line)

        os.remove(nosuff_sam)
        print(f"Cleaned {nosuff2_sam.name}")

for species, path in species_map.items():
    remove_suffixes(species, path)

print("\n All steps complete.")
