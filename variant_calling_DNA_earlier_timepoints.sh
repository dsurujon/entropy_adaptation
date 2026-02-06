#!/usr/bin/env bash


REF="refs/GCF_040687945.1_ASM4068794v1_genomic.fna"
DATA_DIR="demos/data"
ALN_DIR="$DATA_DIR/aln"
ABX="PEN"
OUT_DIR="demos/data/vars/variant_calling_PEN_PRJNA856279"


# Bioproject PRJNA856279
declare -A samples
samples=(
    # Day 7
    ["SRR20002100"]="T4PenD7pop1"
    ["SRR20002099"]="T4PenD7pop2"
    ["SRR20002088"]="T4PenD7pop3"
    ["SRR20002087"]="T4PenD7pop4"
    # Day 14
    ["SRR20002086"]="T4PenD14pop1"
    ["SRR20002085"]="T4PenD14pop2"
    ["SRR20002084"]="T4PenD14pop3"
    ["SRR20002083"]="T4PenD14pop4"
    # Day 21
    ["SRR20002082"]="T4PenD21pop1"
    ["SRR20002081"]="T4PenD21pop2"
    ["SRR20002098"]="T4PenD21pop3"
    ["SRR20002097"]="T4PenD21pop4"
)



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

# Download, align, call variants for each sample
for sample in "${!samples[@]}"; do
    echo "Processing: $sample, Sample name: ${samples[$sample]}"
    fasterq-dump $sample -O $DATA_DIR
    echo "Download complete, aligning reads..."
    bwa mem $REF $DATA_DIR/${sample}_1.fastq $DATA_DIR/${sample}_2.fastq -o $ALN_DIR/${samples[$sample]}.sam -t 10
    echo "Alignment complete, converting SAM to sorted BAM"
    samtools view -bS "$ALN_DIR/${samples[$sample]}.sam" | samtools sort -o "$ALN_DIR/${samples[$sample]}.sorted.bam"
    echo "Sorted BAM created, calling variants"
    call_variants "${samples[$sample]}"
    echo "Variant calling complete for sample: $sample"
done

# remove mutations that were present in the SDMM adapted populations
for sample in "${!samples[@]}"; do
    vcf_file="${OUT_DIR}/${samples[$sample]}_calls.filtered.vcf.gz"
    bcftools isec \
        -C \
        "$vcf_file" \
        "${OUT_DIR}/WT_merged.vcf.gz" \
        -p "${OUT_DIR}/adapted_unique/${samples[$sample]}"
done