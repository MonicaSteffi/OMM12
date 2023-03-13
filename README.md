# Oligo-Mouse-Microbiota(OMM12)synthetic bacterial community

The Oligo-Mouse-Microbiota (OMM12) synthetic bacterial community is a laboratory-made community of 12 bacterial strains that are representative of the mouse gut microbiota. The bacterial strains were selected based on their prevalence and abundance in the gut microbiota of laboratory mice, and their metabolic functions.The OMM12 synthetic community is used as a model system to study the role of the gut microbiota in health and disease, and to investigate the mechanisms underlying the interactions between the host and its microbiota. Researchers can use this community to study how the gut microbiota contributes to host metabolism, immune function, and susceptibility to disease. The community can also be used to evaluate the efficacy of different interventions, such as diet, antibiotics, and probiotics, in modulating the composition and function of the gut microbiota.

Accurate annotation for the OMM12 is important as it provides a framework for understanding the genetic information encoded in a genome. This information can be used to identify genes and their associated functions, as well as to predict the effect of genetic mutations on gene expression and protein function. 

We created this reposiroty to provide proper annotations in order to main uniformity among OMM12 researchers. Genomic and protein sequences are provided for all single organism genome assemblies that are included in NCBI's Assembly resource (www.ncbi.nlm.nih.gov/assembly/). This includes submissions to databases of the International Nucleotide Sequence Database Collaboration, which are available in NCBI's GenBank database, as well 
as the subset of those submissions that are included in NCBI's RefSeq Genomes project. 

# Data Folder

It contains the genome and protein sequences for OligoMM bacterial community. All the sequences were downloaded from the respective databases using inhouse scripts which are availbale in the script folder.


Sub Folders:
1. Downloaded_from_genbank_and_refseq: Genome sequences of OligoMM microbes were downloaded from GenBank and RefSeq databases using inhouser scripts. 
2. Marbouty_genbank_and_refseq_proteins: Contains the protein sequences for each genome, taken from genbank and refseq.
3. Marbouty_genbank_and_refseq_proteins_integrated: Contains the integrated protein sequences for each genome, taken from genbank and refseq.
4. Marbouty_annotations_for_integrated_proteins: Contains the results for the annotations done for the integrated protein sequences from ‘Data\Marbouty_genbank_and_refseq_proteins_integrated’

# Annotation results

We used different softwares such as OperonMapper, KOFAMSCAN, dbCAN2, PULs, RGI and eggNOG to get proper annotations for integrated proteins





