# Oligo-Mouse-Microbiota (OMM12) synthetic bacterial community

The Oligo-Mouse-Microbiota (OMM12) synthetic bacterial community is a laboratory-made community of 12 bacterial strains that are representative of the mouse gut microbiota. The bacterial strains were selected based on their prevalence and abundance in the gut microbiota of laboratory mice, and their metabolic functions.The OMM12 synthetic community is used as a model system to study the role of the gut microbiota in health and disease, and to investigate the mechanisms underlying the interactions between the host and its microbiota. Researchers can use this community to study how the gut microbiota contributes to host metabolism, immune function, and susceptibility to disease. The community can also be used to evaluate the efficacy of different interventions, such as diet, antibiotics, and probiotics, in modulating the composition and function of the gut microbiota.

Accurate annotation for the OMM12 is important as it provides a framework for understanding the genetic information encoded in a genome. This information can be used to identify genes and their associated functions, as well as to predict the effect of genetic mutations on gene expression and protein function. 

We created this reposiroty to provide proper annotations in order to maintain uniformity among OMM12 researchers. We tried different annotation tools like dbCan2, kofam_scan and eggNOG tools and combined them to give an overview. Genomic and protein sequences are provided for all single organism genome assemblies that are included in NCBI's Assembly resource (www.ncbi.nlm.nih.gov/assembly/). This includes submissions to databases of the International Nucleotide Sequence Database Collaboration, which are available in NCBI's GenBank database, as well as the subset of those submissions that are included in NCBI's RefSeq Genomes project. 

## Genome and Poteome Sequences of OMM12 communities

Collection of genomic and proteomic data from Genbank and Refseq database for OMM12 communities can be downloaded from NCBI website as follows:

## Step-1:
## Download the assembly summary genbank NCBI file from ncbi
``
wget https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt
``
## Step-2:

Use the Download_NCBI.R script to download the required NCBI files for the oligoMM latest genome version
# Data Folder

It contains the genome, protein sequences for OligoMM bacterial community as well as the annotations derived from different tools. All the sequences were downloaded from the respective databases using inhouse scripts which are availbale in the script folder.

## Sub Folders:

1. **Marbouty_genbank_and_refseq_proteins:** The data in Marbouty_genbank_and_refseq_proteins is a collection of protein sequences for each genome in the OMM12 communities. These protein sequences have been taken from both GenBank and RefSeq databases. For example, the file name “muribaculum_intestinales_YL27_protein_refseq.faa” has protein sequences for muribaculum intestinales YL27 as given by refseq database.

2. **Marbouty_genbank_and_refseq_proteins_integrated:** All the unique proteins from RefSeq and Genbank are stored in this folder. For example,
“escherichia_coli_strain_Mt1B1_proteins_integrated.faa” file has all unique proteins from Ecoli, in faa format.  The same set of proteins are stored in R-readable file “escherichia_coli_strain_Mt1B1_proteins_integrated.rds”. The corresponding R script named **Genbank_RefSeq_proteins_integrate_based_on_coordinates.R** can be found in script folder. 

3. **Marbouty_annotations_for_integrated_proteins:** Contains the annotations for the integrated protein sequences from ‘Data\Marbouty_genbank_and_refseq_proteins_integrated’

## Annotation results

We used different softwares such as OperonMapper, KOFAMSCAN, dbCAN2, PULs, RGI and eggNOG to get proper annotations for integrated proteins

## Sub Folders:

**1. dbCAN2**

dbCAN2 is a database of Carbohydrate-Active Enzymes (CAZymes) that provides a comprehensive set of tools for automated CAZyme annotation in genomic and metagenomic datasets. 

Uploaded the genomes to the dbCAN2 web page (https://bcb.unl.edu/dbCAN2/blast.php). Chose the options ‘Nucleotide sequence’ for ‘Choose sequence type’ and selected all tools to run (dbCAN (E-Value < 1e-15, coverage > 0.35)  DIAMOND: CAZy (E-Value < 1e-102)  HMMER: dbCAN-sub (E-Value < 1e-15, coverage > 0.35)  CGCFinder (Distance <= 2, signature genes = CAZyme+TC))

CGC finder tool predicts CAZyme gene clusters. The gene clusters can be combination of CAZyme + Transcription Factors (type-i), CAZyme + Transporter Classification (type-II), CAZyme + Transporter Classification + Transcription Factors (type-iii).

## Output files:

* signalp.out: Output of the signalp tool
* uniInput: Prodigal predictions
* overview.txt: Contains overview of the results
* hmmer.out: Output from hmmer tool (with evalue <1e-15 and coverage > 0.35)
* h.out: Raw output from hmmer tool
* diamond.out: Output from diamond tool (with evalue <1e-102)
* dbsub.out: Output from dbCAN_sub database
* cgc_tc_tf: This folder contains outputs from CGC finder tool of type-iii
* cgc_tc: This folder contains outputs from CGC finder tool of type-ii
* cgc_tf: This folder contains outputs from CGC finder tool of type-i

The output prediction for each OMM12 communities stored under separate folders. 

**2. eggNOGDB**

eggNOG (evolutionary genealogy of genes: Non-supervised Orthologous Groups) v5 (https://eggnogdb.embl.de/) is a database of non-supervised orthologous groups and functional annotation. It is a continuation of the eggNOG v4 database, and provides a comprehensive annotation of orthologous groups across a wide range of organisms. The main purpose of eggNOG v5 is to provide a functional annotation for genes and proteins based on orthology. It classifies genes and proteins into orthologous groups based on their evolutionary relationships, and Furthermore, eggNOG v5 integrates data from other databases, such as UniProt, Ensembl, and RefSeq, to provide more comprehensive functional annotations for the genes. It also includes pre-calculated multiple sequence alignments for each orthologous group, which can be used for phylogenetic analysis and other comparative genomics studies.




