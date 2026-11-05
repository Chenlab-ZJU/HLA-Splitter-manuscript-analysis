#!/usr/bin/env bash
# Reuse saved barcode selections; this stage does not randomly sample cells.
set -euo pipefail
source "$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)/config.sh"
read -r -a CELL_NUMS <<< "$CELL_NUMS_STR"
read -r -a SAMPLE_IDS <<< "$SAMPLE_IDS_STR"
for pool in "${CELL_NUMS[@]}"; do
  current_out="$OUT_DIR/$pool"
  barcodes_file="$BARCODE_DIR/$pool.txt"
  if [[ "$RERUN_EXISTING" != 1 && -s "$current_out/sorted_bam.bam" && -s "$current_out/sorted_bam.bam.bai" ]]; then
    printf 'Existing pool retained: %s\n' "$pool"
    continue
  fi
  [[ -s "$barcodes_file" ]] || { echo "Missing barcode list: $barcodes_file" >&2; exit 1; }
  # Check every donor before producing a partial pool.
  for donor in "${SAMPLE_IDS[@]}"; do
    [[ -s "$BAM_DIR/$donor/possorted_genome_bam.bam" ]] || { echo "Missing donor: $donor" >&2; exit 1; }
  done
  mkdir -p "$current_out"
  temp_bams=()
  for donor in "${SAMPLE_IDS[@]}"; do
    subset="$current_out/${donor}_subset.bam"
    subset-bam --bam "$BAM_DIR/$donor/possorted_genome_bam.bam" \
      -c "$barcodes_file" -o "$subset" --cores "$THREADS"
    temp_bams+=("$subset")
  done
  samtools merge -o "$current_out/merged_bam.bam" "${temp_bams[@]}"
  samtools sort "$current_out/merged_bam.bam" -o "$current_out/sorted_bam.bam"
  samtools index "$current_out/sorted_bam.bam"
  rm -- "${temp_bams[@]}" "$current_out/merged_bam.bam"
done
