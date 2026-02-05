#!/usr/bin/env bash

REF="refs/GCF_040687945.1_ASM4068794v1_genomic.fna"
DATA_DIR="demos/data"
ALN_DIR="$DATA_DIR/aln"
ABX="PEN"
OUT_DIR="demos/data/vars/variant_calling_PEN_PRJNA856279"

mkdir -p "$OUT_DIR"

# download SRA files and align to reference

# Bioproject PRJNA856279
# PEN-adapted TIGR4, populations 1-4, day 32 of culture
fasterq-dump SRR20002096 -O $DATA_DIR  # T4PenD32pop1
fasterq-dump SRR20002095 -O $DATA_DIR  # T4PenD32pop2
fasterq-dump SRR20002094 -O $DATA_DIR  # T4PenD32pop3
fasterq-dump SRR20002093 -O $DATA_DIR  # T4PenD32pop4

echo "downloaded PEN adapted populations data"

bwa mem $REF $DATA_DIR/SRR20002096_1.fastq $DATA_DIR/SRR20002096_2.fastq -o $ALN_DIR/T4PenD32pop1.sam -t 10
bwa mem $REF $DATA_DIR/SRR20002095_1.fastq $DATA_DIR/SRR20002095_2.fastq -o $ALN_DIR/T4PenD32pop2.sam -t 10
bwa mem $REF $DATA_DIR/SRR20002094_1.fastq $DATA_DIR/SRR20002094_2.fastq -o $ALN_DIR/T4PenD32pop3.sam -t 10
bwa mem $REF $DATA_DIR/SRR20002093_1.fastq $DATA_DIR/SRR20002093_2.fastq -o $ALN_DIR/T4PenD32pop4.sam -t 10

echo "aligned PEN adapted populations data"

# NDC-adapted TIGR4 (passaged in SDMM), populations 1-4, day 32 of culture
fasterq-dump SRR20002092 -O $DATA_DIR  # T4SDMMPop1
fasterq-dump SRR20002091 -O $DATA_DIR  # T4SDMMPop2
fasterq-dump SRR20002090 -O $DATA_DIR  # T4SDMMPop3
fasterq-dump SRR20002089 -O $DATA_DIR  # T4SDMMPop4

echo "downloaded NDC adapted populations data"

bwa mem $REF $DATA_DIR/SRR20002092_1.fastq $DATA_DIR/SRR20002092_2.fastq -o $ALN_DIR/T4SDMMPop1.sam -t 10
bwa mem $REF $DATA_DIR/SRR20002091_1.fastq $DATA_DIR/SRR20002091_2.fastq -o $ALN_DIR/T4SDMMPop2.sam -t 10
bwa mem $REF $DATA_DIR/SRR20002090_1.fastq $DATA_DIR/SRR20002090_2.fastq -o $ALN_DIR/T4SDMMPop3.sam -t 10
bwa mem $REF $DATA_DIR/SRR20002089_1.fastq $DATA_DIR/SRR20002089_2.fastq -o $ALN_DIR/T4SDMMPop4.sam -t 10

echo "aligned NDC adapted populations data"


cd $ALN_DIR
# Loop over all SAM files, convert to sorted BAM
for sam in *.sam; do
    # Skip if no SAM files are found
    [[ -e "$sam" ]] || continue
    base="${sam%.sam}"
    # Convert SAM to BAM and sort
    samtools view -bS "$sam" | samtools sort -o "${base}.sorted.bam"
done

echo "converted all SAM to sorted BAM"

# -----------------------------
# Define samples 
# -----------------------------

ADAPT_SAMPLES=(
  "T4PenD32pop1"
  "T4PenD32pop2"
  "T4PenD32pop3"
  "T4PenD32pop4"
)
WT_SAMPLES=(
  "T4SDMMPop1"
  "T4SDMMPop2"
  "T4SDMMPop3"
  "T4SDMMPop4"
)

# -----------------------------
# Function: BAM → filtered VCF
# -----------------------------

call_variants () {
  SAMPLE=$1
  BAM="${ALN_DIR}/${SAMPLE}.sorted.bam"
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
# Annotate adapted variants
# -----------------------------

REF="refs/GCF_040687945.1_ASM4068794v1_genomic.gff"
REF_GENES="refs/GCF_040687945.1_ASM4068794v1_genes.gff"
# extract only gene annotations to avoid redundant CDS annotations
grep -P '\tgene\t' ${REF} > ${REF_GENES}

# vcf -> bed convert
bcftools query \
  -f '%CHROM\t%POS0\t%POS\t%REF\t%ALT\n' \
  ${OUT_DIR}/adapted_unique/0000.vcf > ${OUT_DIR}/adapted_unique/0000.bed
# intersect with annotations
bedtools intersect \
  -a ${OUT_DIR}/adapted_unique/0000.bed \
  -b ${REF_GENES} \
  -wa -wb > ${OUT_DIR}/adapted_unique/0000_with_gene_context.tsv

