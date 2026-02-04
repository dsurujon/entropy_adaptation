#!/usr/bin/env bash

REF="refs/GCF_040687945.1_ASM4068794v1_genomic.fna"
DATA_DIR="demos/data"
ALN_DIR="$DATA_DIR/aln"
ABX="PEN"
OUT_DIR="demos/data/vars/variant_calling_PEN_DNA"

mkdir -p "$OUT_DIR"

# download SRA files and align to reference


# PEN-adapted TIGR4, populations 1-4, day 32 of culture
fasterq-dump SRR20576903 -O $DATA_DIR  # 2420347_P132C
# fasterq-dump SRR20576902 -O $DATA_DIR  # 2420348_P232C  TOO LARGE
fasterq-dump SRR20576901 -O $DATA_DIR  # 2420349_P332C
fasterq-dump SRR20576900 -O $DATA_DIR  # 2420350_P432C

echo "downloaded PEN adapted populations data"

bwa mem $REF $DATA_DIR/SRR20576903_1.fastq $DATA_DIR/SRR20576903_2.fastq -o $ALN_DIR/T4-P1-D32.sam -t 10
# bwa mem $REF $DATA_DIR/SRR20576902_1.fastq $DATA_DIR/SRR20576902_2.fastq -o $ALN_DIR/T4-P2-D32.sam -t 10
bwa mem $REF $DATA_DIR/SRR20576901_1.fastq $DATA_DIR/SRR20576901_2.fastq -o $ALN_DIR/T4-P3-D32.sam -t 10
bwa mem $REF $DATA_DIR/SRR20576900_1.fastq $DATA_DIR/SRR20576900_2.fastq -o $ALN_DIR/T4-P4-D32.sam -t 10

echo "aligned PEN adapted populations data"

# NDC-adapted TIGR4, populations 1-4, day 32 of culture
fasterq-dump SRR20576899 -O $DATA_DIR  # 2420351_N132C
fasterq-dump SRR20576898 -O $DATA_DIR  # 2420352_N232C
fasterq-dump SRR20576907 -O $DATA_DIR  # 2420353_N332C
# fasterq-dump SRR20576906 -O $DATA_DIR  # 2420354_N432C   TOO LARGE

echo "downloaded NDC adapted populations data"

bwa mem $REF $DATA_DIR/SRR20576899_1.fastq $DATA_DIR/SRR20576899_2.fastq -o $ALN_DIR/T4-N1-D32.sam -t 10
bwa mem $REF $DATA_DIR/SRR20576898_1.fastq $DATA_DIR/SRR20576898_2.fastq -o $ALN_DIR/T4-N2-D32.sam -t 10
bwa mem $REF $DATA_DIR/SRR20576907_1.fastq $DATA_DIR/SRR20576907_2.fastq -o $ALN_DIR/T4-N3-D32.sam -t 10
# bwa mem $REF $DATA_DIR/SRR20576906_1.fastq $DATA_DIR/SRR20576906_2.fastq -o $ALN_DIR/T4-N4-D32.sam -t 10

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