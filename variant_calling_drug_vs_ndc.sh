#!/usr/bin/env bash
cd demos

REF="../refs/GCF_040687945.1_ASM4068794v1_genomic.fna"
DATA_DIR="./data"
ALN_DIR="$DATA_DIR/aln/T4P"
OUT_DIR="./data/vars/T4P/"
mkdir -p $ALN_DIR
mkdir -p $OUT_DIR

# PEN adapted, NDC
fasterq-dump SRR9060352 -O ./data/
fasterq-dump SRR9060353 -O ./data/
fasterq-dump SRR9059886 -O ./data/
bwa mem $REF $DATA_DIR/SRR9060352.fastq -o $ALN_DIR/T4P-NDC90min-A.sam -t 10
bwa mem $REF $DATA_DIR/SRR9060353.fastq -o $ALN_DIR/T4P-NDC90min-B.sam -t 10
bwa mem $REF $DATA_DIR/SRR9059886.fastq -o $ALN_DIR/T4P-NDC90min-C.sam -t 10

# PEN adapted, PEN
fasterq-dump SRR9060160 -O ./data/
fasterq-dump SRR9060363 -O ./data/
fasterq-dump SRR9060362 -O ./data/
bwa mem $REF $DATA_DIR/SRR9060160.fastq -o $ALN_DIR/T4P-PEN90min-A.sam -t 10
bwa mem $REF $DATA_DIR/SRR9060363.fastq -o $ALN_DIR/T4P-PEN90min-B.sam -t 10
bwa mem $REF $DATA_DIR/SRR9060362.fastq -o $ALN_DIR/T4P-PEN90min-C.sam -t 10

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

PEN_SAMPLES=(
  "T4P-PEN90min-A"
  "T4P-PEN90min-B"
  "T4P-PEN90min-C"
)
NDC_SAMPLES=(
  "T4P-NDC90min-A"
  "T4P-NDC90min-B"
  "T4P-NDC90min-C"
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

PEN_VCFS=()
for S in "${PEN_SAMPLES[@]}"; do
  call_variants "$S"
  PEN_VCFS+=("${OUT_DIR}/${S}_calls.filtered.vcf.gz")
done



bcftools merge "${NDC_VCFS[@]}" -Ov -o "${OUT_DIR}/NDC_merged.vcf"
bcftools merge "${PEN_VCFS[@]}" -Ov -o "${OUT_DIR}/PEN_merged.vcf"

bgzip -f "${OUT_DIR}/NDC_merged.vcf"
bgzip -f "${OUT_DIR}/PEN_merged.vcf"
bcftools index -f "${OUT_DIR}/NDC_merged.vcf.gz"
bcftools index -f "${OUT_DIR}/PEN_merged.vcf.gz"


# -----------------------------
# Subtract NDC from PEN
# -----------------------------

bcftools isec \
  -C \
  "${OUT_DIR}/PEN_merged.vcf.gz" \
  "${OUT_DIR}/NDC_merged.vcf.gz" \
  -p "${OUT_DIR}/PEN_unique"

echo "Done. Final variants:"
echo "${OUT_DIR}/PEN_unique/0000.vcf"


# -----------------------------
# Annotate  variants
# -----------------------------

REF="../refs/GCF_040687945.1_ASM4068794v1_genomic.gff"
REF_GENES="../refs/GCF_040687945.1_ASM4068794v1_genes.gff"

# vcf -> bed convert
bcftools query \
  -f '%CHROM\t%POS0\t%POS\t%REF\t%ALT\n' \
  ${OUT_DIR}/PEN_unique/0000.vcf > ${OUT_DIR}/PEN_unique/0000.bed
# intersect with annotations
bedtools intersect \
  -a ${OUT_DIR}/PEN_unique/0000.bed \
  -b ${REF_GENES} \
  -wa -wb > ${OUT_DIR}/PEN_unique/0000_with_gene_context.tsv

rm $DATA_DIR/*.fastq
rm $ALN_DIR -r