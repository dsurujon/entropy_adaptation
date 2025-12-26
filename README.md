# Adaptation analysis on the entropy dataset 

This repo contains code used in the manuscript "RNA-seq reanalysis identifies murE as a convergent target of 
antibiotic adaptation in Streptococcus pneumoniae"    

### Environment management
* Use [poetry](https://python-poetry.org/docs/basic-usage/)
* Make sure you have a copy of the `poetry.lock` file, and from the same directory run `poetry init`
* To use jupyter notebooks use `poetry run jupyter lab --allow-root`
  
### Non-python dependencies: 
Refer to `nonpython_install.sh`

## Steps to recreate the analysis
1. Download reference     
   `datasets download genome accession GCF_040687945.1`
2. Download fastqs and run alignments against the reference genome     
    `download_and_align_fastq.sh`
3. Call variants     
   `bash variant_calling_abx_adapt.sh`     
   This script will turn bam files into vcf, extract variants in each sample (3 wildtype, 3 adapted), and take the difference between adapted - wildtype to get to a list of mutations uniquely seen in the adapted strain. 
4. Annotate variants - add the names of the genes from the reference    
   `bash annotate_adapt_vcf.sh`
5. Post processing on jupyter notebooks under `demos`    
   1. `Adaptive_mutations`: Generate Table 1 (list of all mutations)
   2. `Coverage_summary`: Generate coverage plot for Figure 1, Table 2 (allele depths)