#!/usr/bin/env bash
cd demos

REF="../refs/GCF_040687945.1_ASM4068794v1_genomic.fna"
DATA_DIR="./data"
ALN_DIR="$DATA_DIR/aln/T4C"
OUT_DIR="./data/vars/T4C/"
mkdir -p $ALN_DIR
mkdir -p $OUT_DIR

# CIP adapted, NDC
fasterq-dump SRR9059534 -O ./data/
fasterq-dump SRR9059533 -O ./data/
fasterq-dump SRR9059528 -O ./data/
bwa mem $REF $DATA_DIR/SRR9059534.fastq -o $ALN_DIR/T4C-NDC120min-A.sam -t 10
bwa mem $REF $DATA_DIR/SRR9059533.fastq -o $ALN_DIR/T4C-NDC120min-B.sam -t 10
bwa mem $REF $DATA_DIR/SRR9059528.fastq -o $ALN_DIR/T4C-NDC120min-C.sam -t 10

# CIP adapted, CIP
fasterq-dump SRR9059524 -O ./data/
fasterq-dump SRR9059523 -O ./data/
fasterq-dump SRR9059352 -O ./data/
bwa mem $REF $DATA_DIR/SRR9059524.fastq -o $ALN_DIR/T4C-CIP120min-A.sam -t 10
bwa mem $REF $DATA_DIR/SRR9059523.fastq -o $ALN_DIR/T4C-CIP120min-B.sam -t 10
bwa mem $REF $DATA_DIR/SRR9059352.fastq -o $ALN_DIR/T4C-CIP120min-C.sam -t 10

echo "aligned all reads to reference, now converting SAM to sorted BAM"

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

CIP_SAMPLES=(
  "T4C-CIP120min-A"
  "T4C-CIP120min-B"
  "T4C-CIP120min-C"
)
NDC_SAMPLES=(
  "T4C-NDC120min-A"
  "T4C-NDC120min-B"
  "T4C-NDC120min-C"
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

cd ../../..
# -----------------------------
# Call variants for all samples
# -----------------------------
NDC_VCFS=()
for S in "${NDC_SAMPLES[@]}"; do
  call_variants "$S"
  NDC_VCFS+=("${OUT_DIR}/${S}_calls.filtered.vcf.gz")
done

CIP_VCFS=()
for S in "${CIP_SAMPLES[@]}"; do
  call_variants "$S"
  CIP_VCFS+=("${OUT_DIR}/${S}_calls.filtered.vcf.gz")
done



bcftools merge "${NDC_VCFS[@]}" -Ov -o "${OUT_DIR}/NDC_merged.vcf"
bcftools merge "${CIP_VCFS[@]}" -Ov -o "${OUT_DIR}/CIP_merged.vcf"

bgzip -f "${OUT_DIR}/NDC_merged.vcf"
bgzip -f "${OUT_DIR}/CIP_merged.vcf"
bcftools index -f "${OUT_DIR}/NDC_merged.vcf.gz"
bcftools index -f "${OUT_DIR}/CIP_merged.vcf.gz"


# -----------------------------
# Subtract NDC from CIP
# -----------------------------

bcftools isec \
  -C \
  "${OUT_DIR}/CIP_merged.vcf.gz" \
  "${OUT_DIR}/NDC_merged.vcf.gz" \
  -p "${OUT_DIR}/CIP_unique"

echo "Done. Final variants:"
echo "${OUT_DIR}/CIP_unique/0000.vcf"


# -----------------------------
# Annotate  variants
# -----------------------------

REF="refs/GCF_040687945.1_ASM4068794v1_genomic.gff"
REF_GENES="refs/GCF_040687945.1_ASM4068794v1_genes.gff"

# vcf -> bed convert
bcftools query \
  -f '%CHROM\t%POS0\t%POS\t%REF\t%ALT\n' \
  ${OUT_DIR}/CIP_unique/0000.vcf > ${OUT_DIR}/CIP_unique/0000.bed
# intersect with annotations
bedtools intersect \
  -a ${OUT_DIR}/CIP_unique/0000.bed \
  -b ${REF_GENES} \
  -wa -wb > ${OUT_DIR}/CIP_unique/0000_with_gene_context.tsv

rm $DATA_DIR/*.fastq
rm $ALN_DIR -r