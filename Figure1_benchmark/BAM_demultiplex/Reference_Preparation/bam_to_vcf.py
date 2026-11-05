# -*- coding: utf-8 -*-
import os
import sys
import subprocess
import argparse
import shutil
from concurrent.futures import ThreadPoolExecutor, as_completed

def check_dependency(tool_name):
    """Check if the tool exists in the PATH."""
    if not shutil.which(tool_name):
        print(f"Error: Tool '{tool_name}' not found. Please ensure it is installed and added to PATH.")
        sys.exit(1)

def run_cmd(cmd, description=None):
    """Execute a Shell command."""
    if description:
        print(f"--- {description} ---")
    print(f"Run: {cmd}")
    try:
        subprocess.run('set -o pipefail; ' + cmd, shell=True, check=True, executable='/bin/bash')
    except subprocess.CalledProcessError as e:
        print(f"Error executing command: {cmd}")
        raise e
        

def process_single_sample(bam_path, ref_genome, out_dir, threads, fast_mode):
    """Process a single sample: Index BAM -> Mpileup -> Call -> Index VCF"""
    
    # 1. Determine Sample ID
    parent_dir = os.path.basename(os.path.dirname(bam_path))
    filename = os.path.basename(bam_path)
    
    if filename == "possorted_genome_bam.bam":
        sample_id = parent_dir
    else:
        sample_id = filename.replace(".bam", "")

    print(f"--- Processing Sample: {sample_id} ---")

    output_vcf_gz = os.path.join(out_dir, f"{sample_id}.vcf.gz")
    
    # Skip if output already exists
    if os.path.exists(output_vcf_gz) and os.path.exists(output_vcf_gz + ".tbi"):
        print(f"  Skipping {sample_id}, output already exists.")
        return output_vcf_gz

    # 2. Check BAM index
    if not (os.path.exists(bam_path + ".bai") or os.path.exists(bam_path + ".csi")):
        print(f"  Indexing BAM: {sample_id}")
        run_cmd(f"samtools index -@ {min(4, threads)} '{bam_path}'")

    # 3. Execute pipeline
    mpileup_opts = "-q 10 -Q 20 -a AD,DP,SP,ADF,ADR -d 1000"
    if fast_mode:
        mpileup_opts += " -B"

    cmd = (
        f"bcftools mpileup -f '{ref_genome}' "
        f"--threads {threads} {mpileup_opts} -Ou '{bam_path}' | "
        f"bcftools call -m -v -Oz -o '{output_vcf_gz}'"
    )
    
    try:
        run_cmd(cmd, description=f"Calling variants for {sample_id}")
    except Exception:
        print(f"Error processing {sample_id}, removing partial files.")
        if os.path.exists(output_vcf_gz): os.remove(output_vcf_gz)
        return None

    # 4. Index VCF
    # Note: tabix uses htslib which is sensitive to zlib version
    run_cmd(f"tabix -p vcf '{output_vcf_gz}'", description=f"Indexing VCF for {sample_id}")
    
    return output_vcf_gz

def main():
    parser = argparse.ArgumentParser(description="Batch extract VCFs (Parallel & Fast Mode).")
    
    parser.add_argument("-i", "--inputs", nargs='+', required=True, help="Input BAM files")
    parser.add_argument("-r", "--ref", required=True, help="Reference Genome FASTA")
    parser.add_argument("-o", "--outdir", required=True, help="Output directory")
    
    parser.add_argument("-m", "--merge", action="store_true", help="Merge VCFs")
    parser.add_argument("--rename-map", help="Sample rename map file")
    parser.add_argument("-t", "--threads", type=int, default=5, help="Threads per job")
    parser.add_argument("-j", "--jobs", type=int, default=1, help="Concurrent jobs")
    parser.add_argument("--fast", action="store_true", help="Fast mode (No BAQ)")
    parser.add_argument("--keep-temp", action="store_true", help="Keep temp files")

    args = parser.parse_args()  

    if not os.path.exists(args.ref):
        sys.exit(f"Error: Reference genome not found: {args.ref}")
        
    os.makedirs(args.outdir, exist_ok=True)

    # Dependency check
    check_dependency("bcftools")
    check_dependency("samtools")

    # 1. Index reference genome
    if not os.path.exists(args.ref + ".fai"):
        run_cmd(f"samtools faidx '{args.ref}'", description="Indexing Ref Genome")

    # 2. Batch process samples
    successful_vcfs = []
    print(f"Found {len(args.inputs)} BAM files.")
    print(f"Parallel Config: {args.jobs} jobs, {args.threads} threads/job.")
    
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:

        future_to_bam = {
            executor.submit(process_single_sample, bam, args.ref, args.outdir, args.threads, args.fast): bam 
            for bam in args.inputs
        }
        
        for future in as_completed(future_to_bam):
            bam = future_to_bam[future]
            try:
                vcf_path = future.result()
                if vcf_path:
                    successful_vcfs.append(vcf_path)
            except Exception as exc:
                print(f"Sample {bam} exception: {exc}")

    if not successful_vcfs:
        sys.exit("No VCFs generated. Exiting.")

    # 3. Merge VCFs
    merged_vcf = os.path.join(args.outdir, "all_samples_merged.vcf.gz")
    
    if args.merge:
        print("\n--- Merging VCFs ---")
        vcf_list_str = " ".join([f"'{v}'" for v in successful_vcfs])
        
        merge_cmd = (
            f"bcftools merge -Oz -o '{merged_vcf}' "
            f"--threads {args.threads} {vcf_list_str}"
        )
        run_cmd(merge_cmd, description="Running bcftools merge")
        run_cmd(f"tabix -p vcf '{merged_vcf}'", description="Indexing merged VCF")

        # 4. Rename samples
        if args.rename_map:
            if not os.path.exists(args.rename_map):
                print(f"Warning: Map file {args.rename_map} not found.")
            else:
                renamed_vcf = merged_vcf.replace(".vcf.gz", ".renamed.vcf.gz")
                reheader_cmd = (
                    f"bcftools reheader --samples '{args.rename_map}' "
                    f"-o '{renamed_vcf}' '{merged_vcf}'"
                )
                run_cmd(reheader_cmd)
                run_cmd(f"tabix -p vcf '{renamed_vcf}'", description="Indexing renamed VCF")
                print(f"Renamed VCF: {renamed_vcf}")

        # 5. Cleanup
        if not args.keep_temp:
            print("\nCleaning up individual VCF files...")
            for vcf in successful_vcfs:
                if os.path.exists(vcf): os.remove(vcf)
                if os.path.exists(vcf + ".tbi"): os.remove(vcf + ".tbi")

    print("\nAll tasks finished successfully.")

if __name__ == "__main__":
    main()
