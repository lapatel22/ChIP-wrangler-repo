#!/usr/bin/env python3
import os
import glob
import subprocess
import argparse
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed


def generate_species_bams(
    target_species: str,
    output_dir: str,
    spike1_species: str,
    spike2_species: str,
    mapq_threshold: int = 50,
    threads: int = 1,
    force_overwrite: bool = False,
):

    # ----------------------- Setup paths -----------------------
    output_dir = Path(output_dir)

    input_dir = output_dir / "concat_align" / "dedup_out"
    output_dir = output_dir / "bams_chr_sep"
    filtered_dir = output_dir / "filtered_bams"

    bamutil_path = "bam"       # bamUtil binary
    samtools_path = "samtools" # samtools binary
    suffix_to_remove = ".nodup.bam"

    output_dir.mkdir(parents=True, exist_ok=True)
    filtered_dir.mkdir(parents=True, exist_ok=True)

    # ----------------------- Step 1: Filter low MAPQ -----------------------
    print("=== Step 1: Filtering low-quality reads ===")
    bam_files = glob.glob(str(input_dir / f"*{suffix_to_remove}"))

    def filter_bam(bam_file, min_mapq):
        bam_file = Path(bam_file)
        base = bam_file.name.replace(suffix_to_remove, "")
        filtered_bam = filtered_dir / f"{base}.filtered.bam"

        if filtered_bam.exists() and not force_overwrite:
            return f"Skipping {filtered_bam.name} (already exists)"

        cmd = [
            samtools_path, "view",
            "-b", "-q", str(min_mapq), "-F", "4",
            "-o", str(filtered_bam),
            str(bam_file)
        ]
        subprocess.run(cmd, check=True)
        return f"Filtered {bam_file.name} → {filtered_bam.name}"

    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = [executor.submit(filter_bam, f, mapq_threshold) for f in bam_files]
        for fut in as_completed(futures):
            print(fut.result())

# ----------------------- Step 2: Split BAMs by chromosome -----------------------
    print("\n=== Step 2: Splitting BAMs by chromosome ===")

    filtered_bams = list(filtered_dir.glob("*.filtered.bam"))

    def split_bam(input_bam):
        input_bam = Path(input_bam)
        base = input_bam.stem.replace(".filtered", "")
        output_prefix = output_dir / f"{base}."
        log_file = output_dir / f"{base}.splitting.log"

        sentinel_chrX = output_dir / f"{base}.chrX.bam"
        if sentinel_chrX.exists() and not force_overwrite:
            return f"Skipping {base} (already split)"

        cmd = [
            bamutil_path,
            "splitChromosome",
            "--in", str(input_bam),
            "--out", str(output_prefix)
        ]
        with open(log_file, "w") as log:
            subprocess.run(cmd, stdout=subprocess.PIPE, stderr=log, check=True)
        return f"Finished splitting {base}"
        

    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = [executor.submit(split_bam, bam) for bam in filtered_bams]
        for fut in as_completed(futures):
            print(fut.result())
    
    # ----------------------- Step 2b: Remove unwanted chromosomes -----------------------
    print("\n=== Step 2b: Removing unwanted contigs (alt, Un, random) ===")

    for bam in filtered_bams:
        base = bam.stem.replace(".filtered", "")
        for bad_bam in output_dir.glob(f"{base}.chr*"):
            if any(x in bad_bam.name for x in ["_alt", "Un", "random", "chrM"]):
                # dont print the unwanted contigs! too many!
        #        print(f"Removing {bad_bam.name} (unwanted contig)")
                bad_bam.unlink()

    # ----------------------- Step 3: Merge by species -----------------------
    print("\n=== Step 3: Merging chromosomes by species ===")

    species_map = {
        spike1_species: output_dir / f"{spike1_species}_data" / f"{spike1_species}_aligned",
        target_species: output_dir / f"{target_species}_data" / f"{target_species}_aligned",
    }

    # Only include spike2 if it's real
    spike2_present = spike2_species.lower() not in ["none", "na", "null", "0", ""]
    if spike2_present:
        species_map[spike2_species] = output_dir / f"{spike2_species}_data" / f"{spike2_species}_aligned"

    # Create directories
    for path in species_map.values():
        path.mkdir(parents=True, exist_ok=True)

    # spike chromosome folder
    spike_chr_dir = output_dir / "spike_bams_chr_sep"
    spike_chr_dir.mkdir(exist_ok=True)

    # Merge per species
    for bam in filtered_bams:
        base = bam.stem.replace(".filtered", "")

        # Move spike-in chromosome BAMs
        for species in [spike1_species, spike2_species] if spike2_present else [spike1_species]:
            for f in output_dir.glob(f"{base}.chr*_{species}.bam"):
                f.rename(spike_chr_dir / f.name)

        outputs = {
            spike1_species: species_map[spike1_species] / f"{base}.{spike1_species}.bam",
            target_species: species_map[target_species] / f"{base}.{target_species}.bam",
        }
        if spike2_present:
            outputs[spike2_species] = species_map[spike2_species] / f"{base}.{spike2_species}.bam"

        # Merge spike species
        for species in outputs:
            if species == target_species:
                continue

            spike_files = sorted(spike_chr_dir.glob(f"{base}.chr*_{species}.bam"))
            out_bam = outputs[species]
            if out_bam.exists() and not force_overwrite:
                print(f"Skipping {out_bam.name} (already exists)")
                continue
            if not spike_files:
                print(f"No chromosome files found for {base} / {species}")
                continue

            print(f"Merging {species} BAMs for {base}")
            subprocess.run([samtools_path, "merge", str(out_bam)] +
                           [str(f) for f in spike_files], check=True)

        # Merge target species
        target_files = sorted(output_dir.glob(f"{base}.chr*.bam"))
        out_bam = outputs[target_species]
        if out_bam.exists() and not force_overwrite:
            print(f"Skipping {out_bam.name} (already exists)")
            continue
        if not target_files:
            print(f"No chromosome files found for {base} / {target_species}")
            continue

        print(f"Merging {target_species} BAMs for {base}")
        subprocess.run([samtools_path, "merge", str(out_bam)] +
                       [str(f) for f in target_files], check=True)


    # ----------------------- Step 4: Remove suffixes -----------------------
    print("\n=== Step 4: Removing species suffixes in final SAM ===")

    def remove_suffixes(species, dir_path):
        for bam in dir_path.glob(f"*.{species}.bam"):
            base = bam.stem
            nosuff_sam = dir_path / f"{base}.nosuff.sam"
            nosuff2_sam = dir_path / f"{base}.nosuffx2.sam"

            if nosuff2_sam.exists() and not force_overwrite:
                print(f"Skipping {nosuff2_sam.name} (already done)")
                continue

            subprocess.run([samtools_path, "view", "-h", "-o", str(nosuff_sam), str(bam)], check=True)

            with open(nosuff_sam) as fin, open(nosuff2_sam, "w") as fout:
                for line in fin:
                    line = line.replace(f"_{spike1_species}", "")
                    if spike2_present:
                        line = line.replace(f"_{spike2_species}", "")
                    fout.write(line)

            os.remove(nosuff_sam)
            print(f"Cleaned {nosuff2_sam.name}")

    for species, path in species_map.items():
        remove_suffixes(species, path)

    print("\nAll steps complete.")


# ----------------------- CLI WRAPPER -----------------------

def main():
    parser = argparse.ArgumentParser(description="Split and filter BAMs by species.")
    parser.add_argument("--spike1_species", required=True, help="Spike-in species 1 (e.g., dm6)")
    parser.add_argument("--spike2_species", required=True, help="Spike-in species 2 (e.g., sacCer3)")
    parser.add_argument("--target_species", required=True, help="Target species (e.g., hg38)")
    parser.add_argument("--output_dir", type=Path, default=Path.cwd(), help="Working directory (default: current working directory)")
    parser.add_argument("--threads", type=int, default=1, help="Threads for BAM processing")
    parser.add_argument("--mapq", type=int, default=50, help="Minimum MAPQ to keep")
    parser.add_argument("--force_overwrite", action="store_true",
    help="Force all steps to run, even if output files already exist.")

    args = parser.parse_args()

    generate_species_bams(
        spike1_species=args.spike1_species,
        spike2_species=args.spike2_species,
        target_species=args.target_species,
        output_dir=args.output_dir,
        mapq_threshold=args.mapq,
        threads=args.threads,
        force_overwrite=args.force_overwrite,
    )


if __name__ == "__main__":
    main()
