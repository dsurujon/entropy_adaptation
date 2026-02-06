#!/usr/bin/env bash


REF="refs/GCF_040687945.1_ASM4068794v1_genomic.fna"
REF_GENES="refs/GCF_040687945.1_ASM4068794v1_genes.gff"

DATA_DIR="demos/data"
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
    # Day 32
    ["SRR20002096"]="T4PenD32pop1"
    ["SRR20002095"]="T4PenD32pop2"
    ["SRR20002094"]="T4PenD32pop3"
    ["SRR20002093"]="T4PenD32pop4"
)

VCFS=()
for sample in "${!samples[@]}"; do
    vcf_file="${OUT_DIR}/adapted_unique/${samples[$sample]}/0000.vcf"
    bgzip -f "$vcf_file"
    bcftools index -f "${vcf_file}.gz"
    vcf_file_final="${vcf_file}.gz"
    VCFS+=($vcf_file_final)
done

# merge vcf files from all populations
bcftools merge "${VCFS[@]}" -Ov \
    -o "${OUT_DIR}/adapted_unique/All_pops_merged.vcf"

# vcf -> bed convert
bcftools query \
  -f '%CHROM\t%POS0\t%POS\t%REF\t%ALT\n' \
  ${OUT_DIR}/adapted_unique/All_pops_merged.vcf > ${OUT_DIR}/adapted_unique/All_pops_merged.bed
  
# intersect with annotations - use left outer join to keep all mutations
bedtools intersect \
  -a ${OUT_DIR}/adapted_unique/All_pops_merged.bed \
  -b ${REF_GENES} \
  -loj > ${OUT_DIR}/adapted_unique/All_pops_merged_with_gene_context.tsv