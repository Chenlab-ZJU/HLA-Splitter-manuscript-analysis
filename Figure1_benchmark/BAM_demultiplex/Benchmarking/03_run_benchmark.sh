#!/usr/bin/env bash
# Run selected tools on existing pools; preserve nonempty prediction files.
set -euo pipefail
source "$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)/config.sh"
selected() { [[ " $TOOLS " == *" $1 "* ]]; }
needs_run() { [[ "$RERUN_EXISTING" == 1 || ! -s "$1" ]]; }
read -r -a methods <<< "$TOOLS"
for tool in "${methods[@]}"; do
  case "$tool" in HLA-Splitter|Vireo_free|Vireo_genotype|Demuxlet|Souporcell) ;; *) echo "Unknown tool: $tool" >&2; exit 2 ;; esac
done
read -r -a pools <<< "$CELL_NUMS_STR"
for pool in "${pools[@]}"; do
  dir="$OUT_DIR/$pool"
  bam="$dir/sorted_bam.bam"
  barcodes="$BARCODE_DIR/$pool.txt"
  [[ -s "$bam" && -s "$barcodes" ]] || { echo "Missing pool inputs: $pool" >&2; exit 1; }
  if selected HLA-Splitter && needs_run "$dir/HLA_Splitter_res/HLA_demultiplex.csv"; then
    output="$dir/HLA_Splitter_res"
    mkdir -p "$output"
    gzip -c "$barcodes" > "$barcodes.gz"
    graph_args=()
    [[ "$HLA_GRAPH_CORR" != 1 ]] || graph_args=(--graph_corr)
    /usr/bin/time -v HLA-Splitter -B "$bam" -b "$barcodes.gz" \
      -L "$HLA_LIST" -H "$HLADB" -o "$output" -t "$THREADS" \
      "${graph_args[@]}" > "$output/time.log" 2>&1
  fi
  for mode in free genotype; do
    selected "Vireo_$mode" || continue
    if [[ "$mode" == free ]]; then
      output="$dir/Vireo_res_without_genotypes"
      variants="$COMMON_SNP_VCF"
      donor_args=(-N "$SAMPLE_NUM")
    else
      output="$dir/Vireo_res_with_genotypes"
      variants="$MERGED_DONOR_VCF"
      donor_args=(-d "$MERGED_DONOR_VCF")
    fi
    needs_run "$output/vireo_results/donor_ids.tsv" || continue
    mkdir -p "$output/cellsnp"
    {
      /usr/bin/time -v cellsnp-lite -s "$bam" -b "$barcodes" -O "$output/cellsnp" \
        -R "$variants" -p "$THREADS" --minMAF 0.1 --minCOUNT 20 --cellTAG CB --UMItag UB --gzip
      /usr/bin/time -v vireo -c "$output/cellsnp" "${donor_args[@]}" -o "$output/vireo_results"
    } > "$output/time.log" 2>&1
  done
  if selected Demuxlet && needs_run "$dir/Demuxlet_res/demuxlet_results.best"; then
    output="$dir/Demuxlet_res"
    mkdir -p "$output"
    {
      /usr/bin/time -v popscle dsc-pileup --sam "$bam" --vcf "$MERGED_DONOR_VCF" \
        --tag-group CB --tag-UMI UB --min-MQ 20 --min-BQ 10 --out "$output/pooled.pileup"
      /usr/bin/time -v popscle demuxlet --plp "$output/pooled.pileup" --field GT \
        --vcf "$MERGED_DONOR_VCF" --out "$output/demuxlet_results"
    } > "$output/time.log" 2>&1
  fi
  if selected Souporcell && needs_run "$dir/Souporcell_res/clusters.tsv"; then
    output="$dir/Souporcell_res"
    mkdir -p "$output"
    /usr/bin/time -v "$CONTAINER_CMD" exec "$SOUP_SIF" souporcell_pipeline.py \
      -i "$bam" -b "$barcodes" -f "$REF_GENOME" -t "$THREADS" \
      -o "$output" -k "$SAMPLE_NUM" > "$output/time.log" 2>&1
  fi
  printf 'Selected methods checked/completed: %s\n' "$pool"
done
