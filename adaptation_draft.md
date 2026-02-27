murE as a convergent target of antibiotic adaptation in *Streptococcus pneumoniae*

# Abstract
Antibiotic exposure imposes strong selective pressures on bacterial populations, frequently leading to the evolution of resistance through mutations. In addition to resistance, bacteria may survive antibiotic stress through tolerance mechanisms that involve broader physiological and metabolic adjustments. Presented here is a reproducible pipeline for variant calling from transcriptomic data and an analysis of mutations acquired in *Streptococcus pneumoniae* strains independently adapted to five different antibiotics.       

Several well-established resistance-associated mutations in antibiotic target genes were recovered, validating both the approach, and the evolutionary relevance of the adapted strains. Beyond these expected changes, there were antibiotic-specific, predominantly nonsynonymous mutations in genes involved in translation, sensing, and cell envelope–associated processes. Notably, a nonsynonymous mutation in murE, encoding an essential enzyme in peptidoglycan precursor synthesis, was independently observed in strains adapted to three antibiotics with distinct mechanisms of action.    

# Introduction

Depending on their mechanism of action, antibiotics exert various primary stresses. For instance, fluoroquinolones such as ciprofloxacin (CIP) and levofloxacin (LVX) target DNA gyrase and topoisomerase IV, inhibiting DNA synthesis and inducing DNA damage–associated stress responses including the SOS system (Hooper & Jacoby 2015; Baharoglu & Mazel 2014). Beta‑lactams such as penicillin (PEN) and glycopeptides such as vancomycin (VNC) inhibit cell wall synthesis and integrity by targeting penicillin‑binding proteins or the D‑Ala–D‑Ala terminus of peptidoglycan precursors, respectively (Hakenbeck et al. 2012; Levine 2006). Other antibiotics inhibit protein synthesis (e.g., aminoglycosides like kanamycin), RNA synthesis (e.g., rifampicin), or folate metabolism (e.g., trimethoprim), all of which are central to bacterial growth and reproduction (Kohanski et al. 2010; Wright 2005). Because antibiotic stress often leads to growth arrest (bacteriostatic antibiotics) or cell death (bactericidal antibiotics), secondary consequences include metabolic imbalance and increased cellular stress (Kohanski et al. 2010; Dwyer et al. 2014).

Prolonged exposure to antibiotics over many generations selects for mutations that confer substantial fitness advantages, often through modification of the antibiotic target and resulting resistance (Andersson & Hughes 2010). For example, *Streptococcus pneumoniae* acquires mutations in penicillin‑binding protein genes (pbp1a, pbp2b, pbp2x) during the evolution of beta‑lactam resistance; these alleles are recurrently observed in clinical resistant isolates (Stanhope et al. 2008; Hakenbeck et al. 2012). While target mutations drive canonical resistance, changes in metabolism are associated with tolerance, where bacteria survive antibiotic exposure without classical resistance mutations (Balaban et al. 2019). Tolerance has been conceptualized as a physiological adjustment, often linked to reduced metabolic activity or altered stress responses, that diminishes lethal effects without changing outcomes in traditional, culture based antibiotic susceptibility testing (Balaban et al. 2019; Levin‑Reisman et al. 2017).

In 2020, Zhu et al. published a large transcriptomic dataset profiling antibiotic‑susceptible wildtype and antibiotic‑adapted S. pneumoniae strains across multiple drugs and genetic backgrounds, focusing on generalizable transcriptomic predictors of survival (Zhu et al. 2020). Although adapted strains were included to expand survival phenotypes, no systematic analysis of the genetic changes underlying adaptation in those strains was performed. The adapted strains are reported to have a modest increase in MIC: 1.5-3x the MIC of the ancestral strain after 100-150 generations, well below clinical breakpoints of resistance (Zhu et al., 2020). This study uses the published RNA‑seq data to identify mutations unique to the adapted populations. The RNA-seq based approach is first validated by comparison to whole genome sequencing in the PEN-adapted population. The approach is then expanded to strains adapted to other antibiotics. By doing so, we reveal both antibiotic‑specific and antibiotic-independent evolutionary strategies and connect these genetic changes to previously reported transcriptomic phenotypes.


# Results
## Variant calling from RNA‑seq yields concordant results to WGS

Typically, variant calling is done using whole genome sequencing (WGS) rather than RNA-seq, since RNA-seq is prone to having variable coverage across the genome with minimal coverage outside of coding sequences. However, resistance conferring mutations are often in coding regions of antibiotic-relevant genes (Alcock et al., 2023). In order to evaluate RNA-seq based variant calling as a viable method, variants detected in RNA-seq were compared to those detected in WGS. In this comparison, the WGS data reflects four independent populations of *S. pneumoniae* strain TIGR4 adapted to PEN (Nishimoto et al., 2022), and the RNA-seq data comes from a clone isolated from one of these populations (Zhu et al., 2020). Variants in the WGS data were filtered to exclude any mutations that appeared in strains that were passaged in no-antibiotic conditions (NDC) in order to eliminate mutations that confer a fitness advantage in lab culture, but are irrelevant for antibiotic adaptation. For WGS, variants were called for days 7, 14, 21 and 32 of adaptation for all 4 populations (Figure 1A). The 6 mutations identified in the final (Day 32) populations are all in coding sequences, are shared across the populations, and are concordant with what is reported in Nishimoto et al., 2022 (Figure 1B). The earlier timepoints, which were not reported in that publication, show that the target-site pbp2x Y586S mutation actually appears later (Days 21-32) in all 4 populations, while the other mutations appear earlier. Other than this however, there is remarkable reproducibility across the populations.     

Since the clone used for the RNA-seq experiment was isolated from one of these populations, and since all 6 mutations reached near 100% frequency in the final populations, the same mutations are expected to be present in the RNA-seq data. Indeed, RNA-seq based variant calling on the clone recovered the same 6 mutations with near 100% frequency in all 3 replicates (Figure 1C). There are two additional mutations present: tsaD A65E appeared in all 3 replicates, and upp P10S appeared in only one replicate. The comparator for the RNAseq-based variant calling was the wildtype strain pre-adaptation (i.e. any mutations that were also observed in the wildtype strains are filtered out). The mutation in tsaD was actually observed in both the PEN-adapted and the NDC-adapted populations in the WGS data, and therefore was filtered out. Therefore, its appearance in RNAseq is expected. The mutation in upp on the other hand appears in only one out of 3 replicates, at lower frequency. This example was used to develop a filter to rule out any artifacts and false positives. For variant calling in other adapted populations from RNA-seq data, if a mutation is not present in all 3 replicates at 75% frequency or above, it was filtered out.     

[!img](demos/data/figures/Fig1_PEN_wgs_vs_rnaseq.svg)
Figure 1: Comparison of WGS and RNA-seq derived variant calling. A. WGS-based variants shown on a timecourse of mutation frequencies across 4 independent populations adapting to PEN B. Final mutation frequencies in the adapted populations from WGS data C. Mutation frequencies from an adapted clone from RNA-seq data 

## Mutations in strains adapted to different antibiotics are relevant to antibiotic mechanism of action

Upon seeing the 6/6 concordance across traditional WGS and RNA-seq based variant calling in the PEN-adapted strains, variant calling was performed on strains that were adapted to 4 other antibiotics: ciprofloxacin (CIP), kanamycin (KAN), rifampicin (RIF), and vancomycin (VNC). While the PEN-adapted population's WGS data was published in a separate paper (Nishimoto et al. 2022), the remaining populations' WGS data are not publicly available, necessitating the re-purposing the RNA-seq data for variant calling for those strains. A total of 19 genic, nonsynonymous mutations were discovered across the different adapted strains (Table 1).   

| position | Reference | Alternative allele   | gene | AA change | ABX | Notes              |
|----------|-----|-------|------|-----|-----|--------------------|
| 713637   | C   | A   | murE* |A430S| CIP | cell wall synthesis      |
| 1006733  | C   | A,T   | gyrA |S81Y, S81F | CIP | direct target      |
| 1355418  | G   | T,A,C | parC |S79Y, S79F, S79C| CIP | direct target      |
| 54632    | C   | A     | dusB |W122C| KAN |  tRNA processing  |
| 648148    | T   | C     | ABC806_RS03415 |V75A| KAN |    |
| 713637   | C   | A   | murE |A430S| KAN | cell wall synthesis      |
| 1403007  | T   | G     | ciaH |T241P| KAN | antibiotic sensing |
| 1955003  | C   | T     | rplF |G174D| KAN |   ribosomal protein- target  |
| 38578   | C   | A   | ABC806_RS00205 |A586E| PEN |       |
| 713637   | C   | A   | murE |A430S| PEN | cell wall synthesis      |
| 876067  | T   | C     |  ABC806_RS04600 |F218L| PEN |                    |
| 949275  | CT   | C     | ABC806_RS04990  |del - Frameshift| PEN |                    |
| 1832114  | G   | A     | cps4E|P80L| PEN | capsule metabolism |
| 1847706  | T   | G     | pbp2X|Y586S| PEN |  direct target   |
| 295258  | T   | G   | rpoB |F474C| RIF | direct target   | 
| 295433  | C   | G   | rpoB |I532M| RIF | direct target   | 
| 716010  | C   | T   | ABC806_RS03780 |A327V| VNC |    | 
| 1154962  | C   | T,A   | ftsW |E39K | VNC | cell wall synthesis   |    
| 1832036  | G   | A   | cps4E |S106F| VNC | capsule metabolism   | 
Table 1: Summary of mutations in coding sequences. Antibiotic denotes which antibiotic the strain is adapted to. Mutations that were not present in all replicates, and mutations common to NDC-adapted populations in WGS data were removed. *the murE mutation in the CIP-adapted strain appears in all replicates but at low (70-80%) frequency. It is included in this summary since it was common across multiple antibiotics and may have functional significance. 


In the ciprofloxacin‑adapted strain, mutations in gyrA (S81Y, S81F) and parC (S79Y, S79F, S79C) represent canonical fluoroquinolone resistance; such mutations in the quinolone resistance‑determining region (QRDR) are widely reported to reduce drug binding (Weigel et al. 2001; Hooper & Jacoby 2015).  

In the rifampicin‑adapted strain, two rpoB mutations (F474C and I532M) were detected. Rifampicin inhibits RNA polymerase by binding the β subunit encoded by rpoB, and clinical rifampicin resistance is strongly linked to mutations in conserved regions of rpoB (Campbell et al. 2001; Ferrándiz et al. 2005; Padayachee et al. 1999). The mutations we observe are proximal to conserved domains involved in rifampicin interaction.

The penicillin‑adapted strain contained a nonsynonymous mutation in pbp2x, a principal target of beta‑lactams in S. pneumoniae (Grebe and Hakenbeck 1996, Hakenbeck et al. 2012). The Y586S amino acid change has also been observed in another publication looking at in vitro adaptive evolution of the same strain (Nishimoto 2022).

These resistance‑associated mutations, observed in clinical contexts and in these adapted lab populations, underscore convergent evolution under antibiotic selection.

## Adapted strains harbor mutations beyond direct antibiotic targets

Beyond target site mutations, several adapted strains display mutations in genes involved in broader metabolic or regulatory processes, suggesting secondary adaptations to drug‑induced stress.

In the kanamycin‑adapted strain, a missense mutation was detected in dusB (W122C), which encodes a protein involved in tRNA modification. Mutations in tRNA modification pathways have been implicated in modulating translational fidelity under stress (Yared et al. 2024). In addition to affecting translation, DusB has been implicated in oxidative stress in *Vibrio cholerae* through its intrinsic NADPH oxidase activity (Fruchard et al. 2025). Another mutation in ciaH results in a T241P amino acid change; CiaH is a sensor kinase in the CiaRH two‑component regulatory system, which influences cell wall homeostasis and antibiotic responsiveness (He et al. 2021).

In both PEN and VNC adapted strains, nonsynonymous mutations in the capsule synthesis gene cps4E were seen; P80L for PEN and S106F for VNC. The P80L mutation is also seen in a previous adaptation study under PEN selection (Nishimoto et al. 2022). Changes in capsule regulation and structure have been linked to cell envelope stress responses and antibiotic tolerance (Fernebro et al. 2004).

These mutations are largely antibiotic‑specific and may reflect secondary adaptations to drug-induced physiological stress rather than direct target modification.

## Recurrent mutation of *murE* across multiple antibiotics

A nonsynonymous mutation in *murE* (A430S) appears independently in the KAN, PEN, and CIP adapted strains, representing the only shared nonsynonymous change across multiple antibiotic conditions. *murE* encodes UDP‑MurNAc‑L‑Ala‑D‑Glu–L‑Lys ligase, an essential enzyme in peptidoglycan precursor synthesis (van Heijenoort 2001). *murE* contains 3 main domains (Gordon et al. 2001), and the A430S change lies in the C-terminal domain, exposed to the protein surface (Figure 2A). This region is implicated in the interaction between MurE and MurF, which is the enzyme that catalyzes the subsequent step in peptidoglycan synthesis (Shirakawa et al. 2023).    


      

# Discussion

The identification of canonical resistance mutations in antibiotic‑adapted populations validates the RNA‑seq variant calling approach used here and confirms that the adapted strains indeed underwent meaningful evolutionary responses consistent with the outcomes of antibiotic exposures reported by Zhu et al. (2020). However, this approach has limitations: RNA‑seq–based variant calling is constrained by expression level, and absence of evidence in unexpressed regions should not be taken as evidence of absence. Future studies with whole‑genome sequencing would be needed to comprehensively catalogue variants. That said, the concordance with WGS for the PEN-adapted population suggest that the mutations identified here are likely to be the majority of the whole set of mutations in the adapted strains rather than a small sampling. 

Most mutations identified in this work were non-synonymous, and on antibiotic target or antibiotic-relevant genes. This suggests very strong positive selection present in the adaptation experiments, even though the adapted strains didn’t fully become resistant, but instead presented with a modest increase in MIC (Zhu et al., 2020). Moreover, 3 mutations showed heterogeneity in terms of alternative alleles, including those appearing on antibiotic targets gyrA and parC in CIP. The CIP-adapted strain also showed lower frequency for the murE A430S mutation, meaning all 3 mutations detected in this clone were at <100% frequency. This suggests that there was either non-clonal isolation or potential contamination within the CIP-adapted clone used for the RNA-seq experiment.   

The recurrent murE A430S mutation identified across multiple antibiotic conditions implicates this core metabolic enzyme as a convergent adaptive node. The particular mutation we see (A430S) is on the C-terminal domain, at the interaction surface between MurE and MurF. The interaction is thought to be mediated by the presence of a number of hydrophobic residues on either protein. The change of a small, nonpolar alanine into a polar serine residue may impact the protein-protein interaction between MurE and MurF, destabilizing the protein complex. MurE’s central role in peptidoglycan precursor synthesis ties envelope biogenesis to growth, and adjusting flux through this pathway may represent a metabolic optimization under diverse stress conditions without impacting classical resistance. Moreover, sequence changes in murE in streptococcal species have been shown to impact antibiotic resistance previously (Todorova 2015).    

In summary, we observe genetic and metabolic rewiring in adapted strains involving both antibiotic‑specific and antibiotic‑agnostic strategies. The identification of murE A430S here as a reproducible, convergent adaptive mutation provides a testable hypothesis for pan-antibiotic tolerance mechanisms in S. pneumoniae and motivates future experimental validation.

# Methods
## Data sources and study design
All analyses were performed using publicly available RNA-seq data generated by Zhu et al. (2020), bioproject PRJNA542628. This dataset comprises Streptococcus pneumoniae TIGR4 wildtype and antibiotic-adapted strains exposed to multiple antibiotics with distinct mechanisms of action, including ciprofloxacin (CIP), kanamycin (KAN), penicillin (PEN), rifampicin (RIF), and vancomycin (VNC). Antibiotic-adapted strains were derived through prolonged exposure to sublethal antibiotic concentrations and were shown in the original study to exhibit increased survival and elevated minimum inhibitory concentrations (MICs). The present study focuses exclusively on identifying genetic variants unique to antibiotic-adapted populations relative to the wildtype background. No new sequencing data were generated.

## Reference genome and annotations
All reads were mapped to the Streptococcus pneumoniae TIGR4 reference genome (ASM4068794v1). Coding sequence (CDS) annotations, gene boundaries, and strand information were obtained from the corresponding RefSeq annotation files (GCF_040687945.1).

## Variant calling from WGS and RNA-seq data
Raw reads were aligned to the TIGR4 reference genome using bwa mem, with default parameters. Resulting SAM files were converted to BAM format, sorted, and indexed using SAMtools. Variants were identified from aligned reads using bcftools mpileup followed by bcftools call (Li, 2011). 

To minimize false positives, variants were filtered using the following criteria:
 - Minimum read depth ≥ 5 
 - Minimum base quality ≥ 15
 - Exclusion of variants present in the wildtype strain
 - Present in all 3 replicates at minimum 75% frequency
 - Absent in NDC-adapted strain (from WGS data)

## Variant annotation 
Variants were assigned to genes based on overlap with annotated CDS regions. For variants occurring within coding sequences, genomic coordinates were converted to CDS-relative positions using strand-aware indexing.
Synonymous, nonsynonymous, and frameshift mutations were identified by reconstructing reference and alternate codons using the reference CDS sequence. Amino acid changes were annotated using the standard bacterial codon table. Variants with multiple alternate alleles were annotated exhaustively, with amino acid changes reported for each alternate base.
Frameshift mutations were identified when the length difference between reference and alternate alleles was not divisible by three.
For each coding variant, the corresponding amino acid position and substitution were determined (e.g. S81T, A430V). Synonymous substitutions were explicitly retained and annotated to allow comparison of mutational patterns across antibiotics.

## MurE structural mapping 
The predicted MurE protein structure AF-Q97PS1-F1-v6 was obtained from AlphaFold (Jumper et al., 2021). The structure was visualized using RCSB PDB (RCSB.org). 

## Data availability and reproducibility
All analysis code, intermediate files, and processed variant tables are available in a publicly accessible GitHub repository ([link to be provided upon submission](https://github.com/dsurujon/entropy_adaptation)). Scripts were written in Python and Bash and rely on standard bioinformatics tools. Exact software versions and parameters are documented in the repository to facilitate reproducibility.


# References

Alcock BP, et al. (2023). CARD 2023: expanded curation, support for machine learning, and resistome prediction at the Comprehensive Antibiotic Resistance Database. Nucleic Acids Res. 2023 Jan 6;51(D1):D690-D699. doi: 10.1093/nar/gkac920. PMID: 36263822; PMCID: PMC9825576.

Andersson DI & Hughes D (2010). Antibiotic resistance and its cost: is it possible to reverse resistance? Nat Rev Microbiol 8: 260–271. PMID: 20208551.

Baharoglu Z & Mazel D (2014). SOS, the formidable strategy of bacteria against aggressions. FEMS Microbiol Rev 38: 1126–1145. PMID: 24923554.

Balaban NQ et al. (2019). Definitions and guidelines for research on antibiotic persistence. Nat Rev Microbiol 17: 441–448. PMID: 30980069.

Campbell EA, Korzheva N, Mustaev A, Murakami K, Nair S, Goldfarb A, Darst SA. Structural mechanism for rifampicin inhibition of bacterial rna polymerase. Cell. 2001 Mar 23;104(6):901-12. doi: 10.1016/s0092-8674(01)00286-0. PMID: 11290327.

Dwyer DJ et al. (2014). Antibiotics induce redox-related physiological alterations as part of their lethality. PNAS. 2014;111(20):E2100–E2109. doi:10.1073/pnas.1401876111. PMID:24803433.

Fernebro J, Andersson I, Sublett J, Morfeldt E, Novak R, Tuomanen E, Normark S, Normark BH. Capsular expression in Streptococcus pneumoniae negatively affects spontaneous and antibiotic-induced lysis and contributes to antibiotic tolerance. J Infect Dis. 2004 Jan 15;189(2):328-38. doi: 10.1086/380564. Epub 2003 Dec 30. PMID: 14722899.

Ferrándiz MJ, Ardanuy C, Liñares J, García-Arenzana JM, Cercenado E, Fleites A, de la Campa AG; Spanish Pneumococcal Infection Study Network. New mutations and horizontal transfer of rpoB among rifampin-resistant Streptococcus pneumoniae from four Spanish hospitals. Antimicrob Agents Chemother. 2005 Jun;49(6):2237-45. doi: 10.1128/AAC.49.6.2237-2245.2005. PMID: 15917517; PMCID: PMC1140543.

Fruchard L, Sudol C, Rouard C, Treffkorn-Maurau A, Hardy L, Bos J, Duchateau M, Giai Gianetto Q, Matondo M, Bonhomme F, Thuillier Q, Marchand V, Motorin Y, Bregeon D, Mazel D, Hamdane D, Baharoglu Z. Beyond RNA modification: a novel role for tRNA modifying enzyme in oxidative stress response and metabolism. Nucleic Acids Res. 2025 Nov 26;53(22):gkaf1276. doi: 10.1093/nar/gkaf1276. PMID: 41385321; PMCID: PMC12700106.

Gordon E, Flouret B, Chantalat L, van Heijenoort J, Mengin-Lecreulx D, Dideberg O. Crystal structure of UDP-N-acetylmuramoyl-L-alanyl-D-glutamate: meso-diaminopimelate ligase from Escherichia coli. J Biol Chem. 2001 Apr 6;276(14):10999-1006. doi: 10.1074/jbc.M009835200. Epub 2000 Dec 20. PMID: 11124264.

Grebe T, Hakenbeck R. Penicillin-binding proteins 2b and 2x of Streptococcus pneumoniae are primary resistance determinants for different classes of beta-lactam antibiotics. Antimicrob Agents Chemother. 1996 Apr;40(4):829-34. doi: 10.1128/AAC.40.4.829. PMID: 8849235; PMCID: PMC163214.

Hakenbeck R, Brückner R, Denapaite D, Maurer P. Molecular mechanisms of β-lactam resistance in Streptococcus pneumoniae. Future Microbiol. 2012 Mar;7(3):395-410. doi: 10.2217/fmb.12.2. PMID: 22393892.

He LY, Le YJ, Guo Z, Li S, Yang XY. The Role and Regulatory Network of the CiaRH Two-Component System in Streptococcal Species. Front Microbiol. 2021 Jul 14;12:693858. doi: 10.3389/fmicb.2021.693858. PMID: 34335522; PMCID: PMC8317062.

Hooper DC, Jacoby GA. Mechanisms of drug resistance: quinolone resistance. Ann N Y Acad Sci. 2015 Sep;1354(1):12-31. doi: 10.1111/nyas.12830. Epub 2015 Jul 17. PMID: 26190223; PMCID: PMC4626314.

Jumper J, Evans R, Pritzel A, Green T, Figurnov M, Ronneberger O, Tunyasuvunakool K, Bates R, Žídek A, Potapenko A, Bridgland A, Meyer C, Kohl SAA, Ballard AJ, Cowie A, Romera-Paredes B, Nikolov S, Jain R, Adler J, Back T, Petersen S, Reiman D, Clancy E, Zielinski M, Steinegger M, Pacholska M, Berghammer T, Bodenstein S, Silver D, Vinyals O, Senior AW, Kavukcuoglu K, Kohli P, Hassabis D. Highly accurate protein structure prediction with AlphaFold. Nature. 2021 Aug;596(7873):583-589. doi: 10.1038/s41586-021-03819-2. Epub 2021 Jul 15. PMID: 34265844; PMCID: PMC8371605.

Kohanski MA, Dwyer DJ, Collins JJ. How antibiotics kill bacteria: from targets to networks. Nat Rev Microbiol. 2010 Jun;8(6):423-35. doi: 10.1038/nrmicro2333. Epub 2010 May 4. PMID: 20440275; PMCID: PMC2896384.

Levine DP. Vancomycin: a history. Clin Infect Dis. 2006 Jan 1;42 Suppl 1:S5-12. doi: 10.1086/491709. PMID: 16323120.

Levin-Reisman I, Ronin I, Gefen O, Braniss I, Shoresh N, Balaban NQ. Antibiotic tolerance facilitates the evolution of resistance. Science. 2017 Feb 24;355(6327):826-830. doi: 10.1126/science.aaj2191. Epub 2017 Feb 9. PMID: 28183996.

Li H. A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. Bioinformatics. 2011 Nov 1;27(21):2987-93. doi: 10.1093/bioinformatics/btr509. Epub 2011 Sep 8. PMID: 21903627; PMCID: PMC3198575.

Nishimoto AT, Dao TH, Jia Q, Ortiz-Marquez JC, Echlin H, Vogel P, van Opijnen T, Rosch JW. Interspecies recombination, not de novo mutation, maintains virulence after β-lactam resistance acquisition in Streptococcus pneumoniae. Cell Rep. 2022 Dec 13;41(11):111835. doi: 10.1016/j.celrep.2022.111835. PMID: 36516783; PMCID: PMC9850807.

Padayachee T, Klugman KP. Molecular basis of rifampin resistance in Streptococcus pneumoniae. Antimicrob Agents Chemother. 1999 Oct;43(10):2361-5. doi: 10.1128/AAC.43.10.2361. PMID: 10508007; PMCID: PMC89483.

Shirakawa KT, Sala FA, Miyachiro MM, Job V, Trindade DM, Dessen A. Architecture and genomic arrangement of the MurE-MurF bacterial cell wall biosynthesis complex. Proc Natl Acad Sci U S A. 2023 May 23;120(21):e2219540120. doi: 10.1073/pnas.2219540120. Epub 2023 May 15. PMID: 37186837; PMCID: PMC10214165.

Stanhope MJ, Lefébure T, Walsh SL, Becker JA, Lang P, Pavinski Bitar PD, Miller LA, Italia MJ, Amrine-Madsen H. Positive selection in penicillin-binding proteins 1a, 2b, and 2x from Streptococcus pneumoniae and its correlation with amoxicillin resistance development. Infect Genet Evol. 2008 May;8(3):331-9. doi: 10.1016/j.meegid.2008.02.001. Epub 2008 Feb 14. PMID: 18394970.

Todorova K, Maurer P, Rieger M, Becker T, Bui NK, Gray J, Vollmer W, Hakenbeck R. Transfer of penicillin resistance from Streptococcus oralis to Streptococcus pneumoniae identifies murE as resistance determinant. Mol Microbiol. 2015 Sep;97(5):866-80. doi: 10.1111/mmi.13070. Epub 2015 Jun 19. PMID: 26010014.

van Heijenoort J (2001). Recent advances in the formation of the bacterial peptidoglycan monomer unit. Nat Prod Rep. 2001 Oct;18(5):503-19. doi: 10.1039/a804532a. PMID: 11699883.

Weigel LM et al. (2001). Genetic analyses of mutations contributing to fluoroquinolone resistance in clinical isolates of Streptococcus pneumoniae. Antimicrob Agents Chemother. 2001 Dec;45(12):3517-23. doi: 10.1128/AAC.45.12.3517-3523.2001. PMID: 11709333; PMCID: PMC90862.

Wright GD (2005). Bacterial resistance to antibiotics: enzymatic degradation and modification. Adv Drug Deliv Rev. 2005 Jul 29;57(10):1451-70. doi: 10.1016/j.addr.2005.04.002. PMID: 15950313.

Yared MJ, Marcelot A, Barraud P. Beyond the Anticodon: tRNA Core Modifications and Their Impact on Structure, Translation and Stress Adaptation. Genes (Basel). 2024 Mar 19;15(3):374. doi: 10.3390/genes15030374. PMID: 38540433; PMCID: PMC10969862.

Zhu Z et al. (2020). Entropy of a bacterial stress response is a generalizable predictor for fitness and antibiotic sensitivity.  Nat Commun 11, 4365 (2020). https://doi.org/10.1038/s41467-020-18134-z