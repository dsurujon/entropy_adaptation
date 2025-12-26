#!/usr/bin/env bash


# generate reference genome index
cd refs
bwa index GCF_000013425.1_ASM1342v1_genomic.fna.gz

cd ..



# WT TIGR4
fasterq-dump SRR9059810 -O ./data/
fasterq-dump SRR9060265 -O ./data/
fasterq-dump SRR9060266 -O ./data/
bwa mem ../refs/GCF_040687945.1_ASM4068794v1_genomic.fna ./data/SRR9059810.fastq -o ../aln/T4-NDC90MIN-A.sam -t 10
bwa mem ../refs/GCF_040687945.1_ASM4068794v1_genomic.fna ./data/SRR9060265.fastq -o ../aln/T4-NDC90MIN-B.sam -t 10
bwa mem ../refs/GCF_040687945.1_ASM4068794v1_genomic.fna ./data/SRR9060266.fastq -o ../aln/T4-NDC90MIN-C.sam -t 10

# CIP adapted
fasterq-dump SRR9059534 -O ./data/
fasterq-dump SRR9059533 -O ./data/
fasterq-dump SRR9059528 -O ./data/
bwa mem ../refs/GCF_040687945.1_ASM4068794v1_genomic.fna ./data/SRR9059534.fastq -o ../aln/T4C-NDC120min-A.sam -t 10
bwa mem ../refs/GCF_040687945.1_ASM4068794v1_genomic.fna ./data/SRR9059533.fastq -o ../aln/T4C-NDC120min-B.sam -t 10
bwa mem ../refs/GCF_040687945.1_ASM4068794v1_genomic.fna ./data/SRR9059528.fastq -o ../aln/T4C-NDC120min-C.sam -t 10

# KAN adapted
fasterq-dump SRR9060166 -O ./data/
fasterq-dump SRR9060163 -O ./data/
fasterq-dump SRR9060102 -O ./data/
bwa mem ../refs/GCF_040687945.1_ASM4068794v1_genomic.fna ./data/SRR9060166.fastq -o ../aln/T4K-NDC120min-a.sam -t 10
bwa mem ../refs/GCF_040687945.1_ASM4068794v1_genomic.fna ./data/SRR9060163.fastq -o ../aln/T4K-NDC120min-b.sam -t 10
bwa mem ../refs/GCF_040687945.1_ASM4068794v1_genomic.fna ./data/SRR9060102.fastq -o ../aln/T4K-NDC120min-c.sam -t 10

# PEN adapted
fasterq-dump SRR9060358 -O ./data/
fasterq-dump SRR9060357 -O ./data/
fasterq-dump SRR9060357 -O ./data/
bwa mem ../refs/GCF_040687945.1_ASM4068794v1_genomic.fna ./data/SRR9060358.fastq -o ../aln/T4P-NDC120min-A.sam -t 10
bwa mem ../refs/GCF_040687945.1_ASM4068794v1_genomic.fna ./data/SRR9060357.fastq -o ../aln/T4P-NDC120min-B.sam -t 10
bwa mem ../refs/GCF_040687945.1_ASM4068794v1_genomic.fna ./data/SRR9060357.fastq -o ../aln/T4P-NDC120min-C.sam -t 10

# RIF adapted
fasterq-dump SRR9060315 -O ./data/
fasterq-dump SRR9060316 -O ./data/
fasterq-dump SRR9060330 -O ./data/
bwa mem ../refs/GCF_040687945.1_ASM4068794v1_genomic.fna ./data/SRR9060315.fastq -o ../aln/T4R-NDC120min-A.sam -t 10
bwa mem ../refs/GCF_040687945.1_ASM4068794v1_genomic.fna ./data/SRR9060316.fastq -o ../aln/T4R-NDC120min-B.sam -t 10
bwa mem ../refs/GCF_040687945.1_ASM4068794v1_genomic.fna ./data/SRR9060330.fastq -o ../aln/T4R-NDC120min-C.sam -t 10

# VNC adapted
fasterq-dump SRR9059429 -O ./data/
fasterq-dump SRR9059428 -O ./data/
fasterq-dump SRR9059427 -O ./data/
bwa mem ../refs/GCF_040687945.1_ASM4068794v1_genomic.fna ./data/SRR9059429.fastq -o ../aln/T4V-024V90MIN-A.sam -t 10
bwa mem ../refs/GCF_040687945.1_ASM4068794v1_genomic.fna ./data/SRR9059428.fastq -o ../aln/T4V-024V90MIN-B.sam -t 10
bwa mem ../refs/GCF_040687945.1_ASM4068794v1_genomic.fna ./data/SRR9059427.fastq -o ../aln/T4V-024V90MIN-C.sam -t 10


cd ../aln/
# Loop over all SAM files, convert to sorted BAM
for sam in *.sam; do
    # Skip if no SAM files are found
    [[ -e "$sam" ]] || continue

    base="${sam%.sam}"

    # Convert SAM to BAM and sort
    samtools view -bS "$sam" | samtools sort -o "${base}.sorted.bam"

done