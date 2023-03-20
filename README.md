# Oligo-Mouse-Microbiota (OMM12) synthetic bacterial community

The Oligo-Mouse-Microbiota (OMM12) synthetic bacterial community is a laboratory-made community of 12 bacterial strains that are representative of the mouse gut microbiota. The bacterial strains were selected based on their prevalence and abundance in the gut microbiota of laboratory mice, and their metabolic functions.The OMM12 synthetic community is used as a model system to study the role of the gut microbiota in health and disease, and to investigate the mechanisms underlying the interactions between the host and its microbiota. Researchers can use this community to study how the gut microbiota contributes to host metabolism, immune function, and susceptibility to disease. The community can also be used to evaluate the efficacy of different interventions, such as diet, antibiotics, and probiotics, in modulating the composition and function of the gut microbiota.

Accurate annotation for the OMM12 is important as it provides a framework for understanding the genetic information encoded in a genome. This information can be used to identify genes and their associated functions, as well as to predict the effect of genetic mutations on gene expression and protein function. 

We created this reposiroty to provide proper annotations in order to maintain uniformity among OMM12 researchers. We tried different annotation tools like dbCan2, KOFAMSCAN, OperonMapper and eggNOG tools and combined them to give an overview. 

## Genome and Protein sequences of OMM12 communities

Genomic and protein sequences for all single organism genome assemblies which are used for annotation are retreived from NCBI's Assembly resource (www.ncbi.nlm.nih.gov/assembly/). 
Collection of genomic and protein sequences from Genbank and Refseq database for OMM12 communities can be downloaded from NCBI website as follows:

## Step-1:
## Download the assembly summary genbank NCBI text file from ncbi using wget
``
wget https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt
``
## Step-2:

Use the scripts/R_files/Download_NCBI.R script to download the required NCBI files for the OligoMM latest genome version (updated version: September 2022)

# Data Folder

It contains annotations derived from different tools for each member of OligoMM communities from different tools. 

## Sub Folders:

1. **Marbouty_genbank_and_refseq_proteins_integrated.zip:** All the unique proteins from RefSeq and Genbank are stored in this folder. For example,
“escherichia_coli_strain_Mt1B1_proteins_integrated.faa” file has all unique proteins from Ecoli, in faa format. The corresponding R script named **Genbank_RefSeq_proteins_integrate_based_on_coordinates.R** can be found in scripts/R_files/ folder. 

### Annotation results

We used different softwares such as dbCAN2, eggNOG, OperonMapper and KOFAMSCAN  to get proper annotations for integrated proteins.

**1. dbCAN2**

dbCAN2 is a database of Carbohydrate-Active Enzymes (CAZymes) that provides a comprehensive set of tools for automated CAZyme annotation in genomic and metagenomic datasets. 

Upload the genome sequences to the dbCAN2 web page (https://bcb.unl.edu/dbCAN2/blast.php). Choose the options ‘Nucleotide sequence’ for ‘Choose sequence type’ and selected all tools to run (dbCAN (E-Value < 1e-15, coverage > 0.35), DIAMOND: CAZy (E-Value < 1e-102),  HMMER: dbCAN-sub (E-Value < 1e-15, coverage > 0.35) and CGCFinder (Distance <= 2, signature genes = CAZyme+TC))

CGC finder tool predicts CAZyme gene clusters. The gene clusters can be combination of CAZyme + Transcription Factors (type-i), CAZyme + Transporter Classification (type-II), CAZyme + Transporter Classification + Transcription Factors (type-iii).

## Output files:

* signalp.out: Output of the signalp (signal peptide predictions) tool 
* uniInput: Prodigal predictions
* overview.txt: This file shows an overview of all the tools run. Each annotated protein is displayed along with which tools annotated it and what CAZy family they were annotated in. Each CAZy family is also a link to the CAZy web page for the appropriate family. Along with this, signal peptide predictions are displayed. The full signalp output is avaliable for download at the top of the tab. The table is also available for download along with the gene predictions (if a nucleotide sequence was uploaded
* hmmer.out: Output from hmmer tool (with evalue <1e-15 and coverage > 0.35)
* h.out: Raw output from hmmer tool
* diamond.out: Output from diamond tool (with evalue <1e-102)
* dbsub.out: Output from dbCAN_sub database. (dbCAN-sub is developed as the first comprehensive CAZyme subfamily HMM database (including CBMs) to enable substrate annotation for CAZymes)
* cgc_tc_tf: This folder contains outputs from CGC finder tool of type-iii
* cgc_tc: This folder contains outputs from CGC finder tool of type-ii
* cgc_tf: This folder contains outputs from CGC finder tool of type-i

The output prediction for each OMM12 communities stored under separate folders. Deail information about the input and output files can be found in https://bcb.unl.edu/dbCAN2/help.php. 

**2. eggNOGDB**

eggNOG (evolutionary genealogy of genes: Non-supervised Orthologous Groups) v5 (https://eggnogdb.embl.de/) is a database of non-supervised orthologous groups and functional annotation. It provides a comprehensive annotation of orthologous groups across a wide range of organisms. The main purpose of eggNOG v5 is to provide a functional annotation for genes and proteins based on orthology. It classifies genes and proteins into orthologous groups based on their evolutionary relationships, and Furthermore, eggNOG v5 integrates data from other databases, such as UniProt, Ensembl, and RefSeq, to provide more comprehensive functional annotations for the genes. It also includes pre-calculated multiple sequence alignments for each orthologous group, which can be used for phylogenetic analysis and other comparative genomics studies.

## Output files:

* out.emapper.annotations: These files provide additional information on the functional annotation for each gene in the search results. They include information on gene ontology (GO) terms, enzyme commission (EC) numbers, and protein domains, among other features.
* out.emapper.genepred.fasta: A FASTA file with the sequences of the predicted CDS.
* out.emapper.genepred.gff: A GFF file with the sequences of the predicted CDS. 
* out.emapper.hits: A file with the results from the search phase, from HMMER, Diamond or MMseqs2.
* out.emapper.orthologs: A file with the list of orthologs found for each query. This file is created only if using the --report_orthologs option.
* out.emapper.seed_orthologs: A file with the results from parsing the hits. Each row links a query with a seed ortholog. This file has the same format independently of which searcher was used, except that it can be in short format (4 fields), or full.

**3. KOFAMSCAN**

KOFAMSCAN generates several output files after performing functional annotation of protein-coding genes in a metagenome or metatranscriptome dataset.

Uploaded the integrated proteins from  folder to https://www.genome.jp/tools/kofamkoala/. Downloaded the result file and renamed as ‘_proteins_kofam_results.txt’ and saved the files in the present folder. Since the new version of kofamscan code (exec_annotation ) didn’t finish within the allotted time of 48 hours, the online version of analysis is preferred.

KEGG Version: 103.0
Version release date: 01-08-2022
The data behind the analysis is downloaded and kept in the folder from https://www.genome.jp/ftp/tools/kofam_scan/ and https://www.genome.jp/ftp/db/kofam/ (As on August 29,2022)

Some of the commonly generated output files include:

* kegg_annotation.txt: This file contains the KEGG Orthology (KO) annotations for each predicted protein-coding gene in the input sequence dataset.

**4. OperonMapper**

OperonMapper is a bioinformatics tool that is used to predict operons in prokaryotic genomes. An operon is a group of genes that are transcribed together as a single unit and are involved in a common biological pathway or function. The prediction of operons is important for understanding the functional organization of prokaryotic genomes, as well as for identifying potential gene regulatory mechanisms.

OperonMapper uses a comparative genomics approach to predict operons in prokaryotic genomes. It compares the gene order and orientation of the input genome with those of reference genomes that have experimentally verified operons. Based on this comparison, OperonMapper predicts operons in the input genome and assigns them to specific biological functions or pathways.

OperonMapper provides several output files that can be used to analyze and visualize the predicted operons. These output files include:
* Operonic gene pairs: Predicted operonic gene pairs with their corresponding confidence values of being part of the same operon
* ORFs coordinates: Predicted ORFs coordinates
* Predicted protein sequences: Protein sequences of the translated predicted ORFs
* DNA sequences of the predicted ORFs
* Orthology assignment of proteins to their corresponding COGs groups
* functional_descriptions: Proteins functional descriptions
* List of operons: 
* Predicted COGs:
* predicted orfs: List of operons with their conforming genes

### newJointAnnotation

This folder contains excel sheets for each memeber of OligoMM community with all the annotation derived from different annotation tools.

### Geneius analysis input

This folder contains suitable gff formats derived from different annotation tools such as dbCAN2, KOFAMSCAN, eggNOG and Operon mapper for Geneius visulaization purposes.



