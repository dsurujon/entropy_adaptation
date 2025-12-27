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


# Results

Below is a representative subset of annotated mutations identified in coding regions:
| position | Ref | Alt   | gene | AA change | ABX | Notes              |
|----------|-----|-------|------|-----|-----|--------------------|
| 713637   | C   | A   | murE |A430S| CIP | cell wall synthesis      |
| 1006733  | C   | A,T   | gyrA |S81Y, S81F | CIP | direct target      |
| 1355418  | G   | T,A,C | parC |S79Y, S79F, S79C| CIP | direct target      |
| 54632    | C   | A     | dusB |W122C| KAN |  tRNA processing  |
| 713637   | C   | A   | murE |A430S| KAN | cell wall synthesis      |
| 1403007  | T   | G     | ciaH |T241P| KAN | antibiotic sensing |
| 1955003  | C   | T     | rplF |G174D| KAN |   ribosomal protein- target  |
| 2026741  | G   | T,A   | tsaD |A65E, A65V| KAN | tRNA modification  | 
| 713637   | C   | A   | murE |A430S| PEN | cell wall synthesis      |
| 1450247  | G   | A     | upp  |P10S| PEN |                    |
| 1832114  | G   | A     | cps4E|P80L| PEN | capsule metabolism |
| 1847706  | T   | G     | pbp2X|Y586S| PEN |  direct target   |
| 2026741  | G   | T,A   | tsaD |A65E, A65V| PEN |  tRNA modification   | 
| 41080  | G   | A     | cbpD  |G29R| RIF |  competence induced lysis    |
| 72331  | TC   | T     | dltD|Frameshift| RIF |  cell surface modification |
| 226057  | T   | G     | ulaG|W155G| RIF | vitamin C transport and metabolism |
| 295258  | T   | G   | rpoB |F474C| RIF | direct target   | 
| 295433  | C   | G   | rpoB |I532M| RIF | direct target   | 
| 295258  | T   | G   | ssrA |L116L| VNC | transfer mRNA   | 
| 1832036  | G   | A   | cps4E |S106F| VNC | capsule metabolism   | 
| 2026741  | G   | T,A   | tsaD |A65E, A65V| VNC | tRNA modification   | 


This is not a complete list. Please see the jupyter notebook for the comprehensive list of mutations.     
Across the adaptive evolution experiments, we observe that mutations in the direct targets of the antibiotics arise as expected, including:
- gyrA and parC in fluoroquinolone-adapted strains,
- pbp2X in penicillin-adapted strains,
- rpoB in rifampicin-adapted strains.

Additionally, several mutations occur in genes associated with regulation, stress response, and central metabolism, such as:
- cbpD, implicated in competence and fratricide regulation,
- dltD, involved in modification of teichoic acids,
- ulaG, linked to carbohydrate catabolism, 
- cps4E, in capsule synthesis, and
- tsaD, a tRNA modifying enzyme observed in both KAN, PEN, and VNC adapted backgrounds.

A recurrent mutation at coordinate 713637 affecting ABC806_RS03775 (annotated as UDP-N-acetylmuramoyl-L-alanyl-D-glutamate–L-lysine ligase, a MurE homolog) appears independently in CIP, PEN, and KAN adapted strains. MurE catalyzes an essential step in peptidoglycan precursor synthesis, and its repeated targeting across distinct antibiotics suggests a convergent adaptive adjustment of cell wall biosynthetic flux under heterogeneous antibiotic stress.     

### Genes with interesting mutations
#### murE
Function: MurE is the UDP-MurNAc-tripeptide–L-lysine ligase that adds L-lysine to the peptidoglycan stem peptide, a central step in cell wall biosynthesis.
Adaptation relevance: Altered flux through peptidoglycan synthesis can mitigate imbalances between cell growth and envelope assembly under stress 

#### gyrA / parC
Function: Encode subunits of DNA gyrase (GyrA) and topoisomerase IV (ParC), the canonical fluoroquinolone targets.
Adaptation relevance: Point mutations in these proteins reduce fluoroquinolone binding and confer resistance 

#### dusB
Function: Dihydrouridine synthase B catalyzes formation of dihydrouridine in tRNAs, affecting tRNA flexibility and translational accuracy.
Adaptation relevance: Mutations in tRNA modification enzymes have been associated with altered translation fidelity and antibiotic tolerance 

#### ciaH
Function: Sensor kinase of the two-component CiaRH system, regulating cell wall stress responses, autolysis, and competence.
Adaptation relevance: CiaRH mutants possess altered antibiotic susceptibility profiles and stress resilience 

#### rplF
Function: Ribosomal protein L6 is part of the 50S large subunit and contributes to ribosome stability and function.
Adaptation relevance: Mutations in ribosomal proteins can modulate aminoglycoside interaction and translational dynamics

#### tsaD
Function: Member of the t⁶A tRNA modification pathway, required for synthesis of N⁶-threonylcarbamoyladenosine on ANN-recognizing tRNAs.
Adaptation relevance: Disruption of t⁶A modification alters decoding and can confer growth and stress phenotypes linked to antibiotic responses 

#### upp
Function: Uracil phosphoribosyltransferase in pyrimidine salvage.
Adaptation relevance: Salvage pathway genes show up in diverse adaptation contexts, possibly reflecting nucleotide pool rebalancing under growth stress 

#### cps4E
Function: Capsule biosynthesis enzyme in serotype-specific capsular operons.
Adaptation relevance: Capsule alterations can affect cell surface properties and antibiotic tolerance indirectly, shows up in cell-wall synthesis inhibitors PEN, VNC 

#### pbp2X
Function: Penicillin-binding protein 2X, one of the primary transpeptidases targeted by β-lactams.
Adaptation relevance: Mutations here reduce β-lactam binding and are a well-established mechanism of penicillin resistance in pneumococci 

#### cbpD
Function: Choline-binding protein D, a murein hydrolase that participates in competence-associated fratricide.
Adaptation relevance: Mutations may alter competence dynamics and surface remodeling, potentially affecting survival under chronic stress 

#### dltD
Function: Part of the DltABCD system that adds D-alanine to teichoic acids, reducing cell surface negative charge.
Adaptation relevance: Altered teichoic acid D-alanination changes cell envelope stress tolerance and can modulate susceptibility to stressors 

#### ulaG
Function: Likely part of ascorbate uptake/metabolism (the ula operon).
Adaptation relevance: Mutations here may reflect altered metabolic routing under chronic stress rather than direct drug interaction

#### ssrA
Function: ssrA encodes tmRNA, central to the ribosome rescue system / trans-translation.    
Interpretation: Adaptive response to stalled ribosomes or translation stress.    
Mechanistically, this fits with what you saw in other antibiotics: stress on translation or protein homeostasis pathways prompts convergent selection on supporting systems (tRNA modification + ribosome rescue).

