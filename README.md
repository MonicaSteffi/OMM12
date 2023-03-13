# Oligo-Mouse-Microbiota(OMM12)synthetic bacterial community

The Oligo-Mouse-Microbiota (OMM12) synthetic bacterial community is a laboratory-made community of 12 bacterial strains that are representative of the mouse gut microbiota. The bacterial strains were selected based on their prevalence and abundance in the gut microbiota of laboratory mice, and their metabolic functions.The OMM12 synthetic community is used as a model system to study the role of the gut microbiota in health and disease, and to investigate the mechanisms underlying the interactions between the host and its microbiota. Researchers can use this community to study how the gut microbiota contributes to host metabolism, immune function, and susceptibility to disease. The community can also be used to evaluate the efficacy of different interventions, such as diet, antibiotics, and probiotics, in modulating the composition and function of the gut microbiota.

Accurate annotation for the OMM12 is important as it provides a framework for understanding the genetic information encoded in a genome. This information can be used to identify genes and their associated functions, as well as to predict the effect of genetic mutations on gene expression and protein function. 

We created this reposiroty to provide proper annotations in order to maintain uniformity among OMM12 researchers. We tried different annotation tools like dbCan2, kofam_scan and eggNOG tools and combined them to give an overview. Genomic and protein sequences are provided for all single organism genome assemblies that are included in NCBI's Assembly resource (www.ncbi.nlm.nih.gov/assembly/). This includes submissions to databases of the International Nucleotide Sequence Database Collaboration, which are available in NCBI's GenBank database, as well as the subset of those submissions that are included in NCBI's RefSeq Genomes project. 

# Data Folder

It contains the genome, protein sequences for OligoMM bacterial community as well as the annotations derived from different tools. All the sequences were downloaded from the respective databases using inhouse scripts which are availbale in the script folder.

## Sub Folders:

Collection of genomic and proteomic data from Genbank and Refseq database for OMM12 communities can be downloaded from NCBI website as follows:

## Step-1:
## Download the assembly summary genbank file from ncbi
<wget https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt
wget https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt>

## Step-2:
## Download the genbank NCBI file

### Function to download the required NCBI files for the oligoMM latest genome version
##### Source the R function and then call the function.

The data is organized into 12 different folders, one for each member of the OMM12 communities. Each folder contains the genome sequences for that particular community member. The README files in each directory likely provide information on the fields and formats used in the data, as well as any additional information on the source and processing of the data.
Marbouty_refseq: It also contain similar information as the previous folder. However, the data in Marbouty_refseq appears to have been downloaded from the RefSeq database instead of GenBank 
Files in each folder: 
File formats and content:


Marbouty_genbank_and_refseq_proteins: The data in Marbouty_genbank_and_refseq_proteins is a collection of protein sequences for each genome in the OMM12 communities. These protein sequences have been taken from both GenBank and RefSeq databases. For example, the file name “muribaculum_intestinales_YL27_protein_refseq.faa” has protein sequences for muribaculum intestinales YL27 as given by refseq database.
3. Marbouty_genbank_and_refseq_proteins_integrated: Contains the integrated protein sequences for each genome, taken from genbank and refseq.
4. Marbouty_annotations_for_integrated_proteins: Contains the results for the annotations done for the integrated protein sequences from ‘Data\Marbouty_genbank_and_refseq_proteins_integrated’

# Annotation results

We used different softwares such as OperonMapper, KOFAMSCAN, dbCAN2, PULs, RGI and eggNOG to get proper annotations for integrated proteins





