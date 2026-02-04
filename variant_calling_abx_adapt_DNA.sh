#!/usr/bin/env bash
set -euo pipefail

REF="refs/GCF_040687945.1_ASM4068794v1_genomic.fna"
ALN_DIR="aln"
ABX="PEN_DNA"
OUT_DIR="demos/data/vars/variant_calling_${ABX}"

mkdir -p "$OUT_DIR"

# -----------------------------
# Define samples here
# -----------------------------

WT_SAMPLES=(
  "T4-P1-D32"
  "T4-P2-D32"
  "T4-P3-D32"
)

ADAPT_SAMPLES=(
  "T4-N1-D32"
  "T4-N2-D32"
  "T4-N3-D32"
)

# -----------------------------
# Function: BAM → filtered VCF
# -----------------------------

call_variants () {
  SAMPLE=$1
  BAM="${ALN_DIR}/${SAMPLE}.bam"
  PREFIX="${OUT_DIR}/${SAMPLE}"

  echo "Processing ${SAMPLE}"

  bcftools mpileup \
    -f "$REF" \
    -Q 20 -q 20 \
    -a AD,DP,SP \
    --max-depth 10000 \
    "$BAM" \
    -Ov \
    -o "${PREFIX}.vcf"

  bcftools call \
    "${PREFIX}.vcf" \
    --ploidy 1 \
    -c -v \
    -Ov \
    -o "${PREFIX}_calls.vcf"

  bcftools filter \
    -i 'FORMAT/DP>=5 && QUAL>=15' \
    "${PREFIX}_calls.vcf" \
    -Ov \
    -o "${PREFIX}_calls.filtered.vcf"

  bgzip -f "${PREFIX}_calls.filtered.vcf"
  bcftools index -f "${PREFIX}_calls.filtered.vcf.gz"
}

# -----------------------------
# Call variants for all samples
# -----------------------------

for S in "${WT_SAMPLES[@]}"; do
  call_variants "$S"
done

for S in "${ADAPT_SAMPLES[@]}"; do
  call_variants "$S"
done

# -----------------------------
# Merge WT and adapted (if >1)
# -----------------------------

WT_VCFS=()
for S in "${WT_SAMPLES[@]}"; do
  WT_VCFS+=("${OUT_DIR}/${S}_calls.filtered.vcf.gz")
done

ADAPT_VCFS=()
for S in "${ADAPT_SAMPLES[@]}"; do
  ADAPT_VCFS+=("${OUT_DIR}/${S}_calls.filtered.vcf.gz")
done

bcftools merge "${WT_VCFS[@]}" -Ov -o "${OUT_DIR}/WT_merged.vcf"
bcftools merge "${ADAPT_VCFS[@]}" -Ov -o "${OUT_DIR}/ADAPT_merged.vcf"

bgzip -f "${OUT_DIR}/WT_merged.vcf"
bgzip -f "${OUT_DIR}/ADAPT_merged.vcf"
bcftools index -f "${OUT_DIR}/WT_merged.vcf.gz"
bcftools index -f "${OUT_DIR}/ADAPT_merged.vcf.gz"

# -----------------------------
# Subtract WT from adapted
# -----------------------------

bcftools isec \
  -C \
  "${OUT_DIR}/ADAPT_merged.vcf.gz" \
  "${OUT_DIR}/WT_merged.vcf.gz" \
  -p "${OUT_DIR}/adapted_unique"

echo "Done. Final variants:"
echo "${OUT_DIR}/adapted_unique/0000.vcf"
