#!/usr/bin/env bash
# this script annotates the adapted unique variants called in variant_calling_abx_adapt.sh
# the annotations come from the reference genome GFF file
ABX="VNC"
OUT_DIR="data/vars/variant_calling_${ABX}"
REF="../refs/GCF_040687945.1_ASM4068794v1_genomic.gff"
REF_GENES="../refs/GCF_040687945.1_ASM4068794v1_genes.gff"
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