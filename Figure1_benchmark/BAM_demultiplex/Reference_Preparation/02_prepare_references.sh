#!/usr/bin/env bash
# Review donor names and HLA calls manually before using these references.
set -euo pipefail
source "$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)/config.sh"
shopt -s nullglob
bams=("$BAM_DIR"/*/possorted_genome_bam.bam)
(( ${#bams[@]} > 0 )) || { echo "No donor BAMs found in $BAM_DIR" >&2; exit 1; }
rename_args=()
if [[ -n "$SAMPLE_RENAME_MAP" ]]; then
  [[ -s "$SAMPLE_RENAME_MAP" ]] || { echo 'Missing sample rename map.' >&2; exit 1; }
  rename_args=(--rename-map "$SAMPLE_RENAME_MAP")
fi
mkdir -p "$SNP_OUT_DIR" "$HLA_OUT_DIR"
python "$CODE_DIR/Reference_Preparation/bam_to_vcf.py" \
  -i "${bams[@]}" -r "$REF_GENOME" -o "$SNP_OUT_DIR" \
  -m -j 4 -t 5 --fast --keep-temp "${rename_args[@]}"
python "$CODE_DIR/Reference_Preparation/extract_hla_genotypes.py" \
  --inputs "${bams[@]}" --outdir "$HLA_OUT_DIR" \
  --arcas "$(command -v arcasHLA)" --threads "$THREADS" \
  --genes "A,B,C,DPB1,DQB1,DRB1"
printf '%s\n' 'Inspect HLA_Genotypes_Summary.csv and the donor VCF sample names.' \
  'If renamed VCFs are used, set MERGED_DONOR_VCF to the .renamed.vcf.gz file.' \
  'No demultiplexing is started by this stage.'
