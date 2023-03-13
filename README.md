# Oligo-Mouse-Microbiota(OMM12)synthetic bacterial community

The Oligo-Mouse-Microbiota (OMM12) synthetic bacterial community is a laboratory-made community of 12 bacterial strains that are representative of the mouse gut microbiota. The bacterial strains were selected based on their prevalence and abundance in the gut microbiota of laboratory mice, and their metabolic functions.The OMM12 synthetic community is used as a model system to study the role of the gut microbiota in health and disease, and to investigate the mechanisms underlying the interactions between the host and its microbiota. Researchers can use this community to study how the gut microbiota contributes to host metabolism, immune function, and susceptibility to disease. The community can also be used to evaluate the efficacy of different interventions, such as diet, antibiotics, and probiotics, in modulating the composition and function of the gut microbiota.

Accurate annotation for the OMM12 is important as it provides a framework for understanding the genetic information encoded in a genome. This information can be used to identify genes and their associated functions, as well as to predict the effect of genetic mutations on gene expression and protein function. 

We created this reposiroty to provide proper annotations in order to main uniformity among OMM12 researchers. We tried different annotation tools like dbCan2, kofam_scan and eggNOG tools. Genomic and protein sequences are provided for all single organism genome assemblies that are included in NCBI's Assembly resource (www.ncbi.nlm.nih.gov/assembly/). This includes submissions to databases of the International Nucleotide Sequence Database Collaboration, which are available in NCBI's GenBank database, as well as the subset of those submissions that are included in NCBI's RefSeq Genomes project. 

# Data Folder

It contains the genome, protein sequences for OligoMM bacterial community as well as the annotations derived from different tools. All the sequences were downloaded from the respective databases using inhouse scripts which are availbale in the script folder.

## Sub Folders:
1. Downloaded_from_genbank_and_refseq: 
Marbouty_Genbank: is a collection of genomic data for OMM12 communities that have been downloaded from GenBank. The data is organized into 12 different folders, one for each member of the OMM12 communities. Each folder contains the genome sequences for that particular community member. The README files in each directory likely provide information on the fields and formats used in the data, as well as any additional information on the source and processing of the data.
Marbouty_refseq: It also contain similar information as the previous folder. However, the data in Marbouty_refseq appears to have been downloaded from the RefSeq database instead of GenBank 
Files in each folder: 
File formats and content:
   assembly_status.txt
       A text file reporting the current status of the version of the assembly
       for which data is provided. Any assembly anomalies are also reported.
       
   ***_assembly_report.txt** : Tab-delimited text file reporting the name, role and sequence accession.version for objects in the assembly. The file header contains 
       meta-data for the assembly including: assembly name, assembly accession.version, scientific name of the organism and its taxonomy ID,  assembly submitter, and          sequence release date.
       
   ***_assembly_stats.txt**: Tab-delimited text file reporting statistics for the assembly including: 
       total length, ungapped length, contig & scaffold counts, contig-N50, 
       scaffold-L50, scaffold-N50, scaffold-N75, and scaffold-N90
   *_assembly_regions.txt
       Provided for assemblies that include alternate or patch assembly units. 
       Tab-delimited text file reporting the location of genomic regions and 
       listing the alt/patch scaffolds placed within those regions.
   *_assembly_structure directory
       This directory will only be present if the assembly has internal 
       structure. When present, it will contain AGP files that define how 
       component sequences are organized into scaffolds and/or chromosomes. 
       Other files define how scaffolds and chromosomes are organized into 
       non-nuclear and other assembly-units, and how any alternate or patch 
       scaffolds are placed relative to the chromosomes. Refer to the README.txt
       file in the assembly_structure directory for additional information.
   *_cds_from_genomic.fna.gz
       FASTA format of the nucleotide sequences corresponding to all CDS 
       features annotated on the assembly, based on the genome sequence. See 
       the "Description of files" section below for details of the file format.
   *_feature_count.txt.gz
       Tab-delimited text file reporting counts of gene, RNA, CDS, and similar
       features, based on data reported in the *_feature_table.txt.gz file.
       See the "Description of files" section below for details of the file 
       format.
   *_feature_table.txt.gz
       Tab-delimited text file reporting locations and attributes for a subset 
       of annotated features. Included feature types are: gene, CDS, RNA (all 
       types), operon, C/V/N/S_region, and V/D/J_segment. Replaces the .ptt & 
       .rnt format files that were provided in the old genomes FTP directories.
       See the "Description of files" section below for details of the file 
       format.
   *_genomic.fna.gz file
       FASTA format of the genomic sequence(s) in the assembly. Repetitive 
       sequences in eukaryotes are masked to lower-case (see below).
       The FASTA title is formatted as sequence accession.version plus 
       description. The genomic.fna.gz file includes all top-level sequences in
       the assembly (chromosomes, plasmids, organelles, unlocalized scaffolds,
       unplaced scaffolds, and any alternate loci or patch scaffolds). Scaffolds
       that are part of the chromosomes are not included because they are
       redundant with the chromosome sequences; sequences for these placed 
       scaffolds are provided under the assembly_structure directory.
   *_genomic.gbff.gz file
       GenBank flat file format of the genomic sequence(s) in the assembly. This
       file includes both the genomic sequence and the CONTIG description (for 
       CON records), hence, it replaces both the .gbk & .gbs format files that 
       were provided in the old genomes FTP directories.
   *_genomic.gff.gz file
       Annotation of the genomic sequence(s) in Generic Feature Format Version 3
       (GFF3). Sequence identifiers are provided as accession.version.
       Additional information about NCBI's GFF files is available at 
       ftp://ftp.ncbi.nlm.nih.gov/genomes/README_GFF3.txt.
   *_genomic.gtf.gz file
       Annotation of the genomic sequence(s) in Gene Transfer Format Version 2.2
       (GTF2.2). Sequence identifiers are provided as accession.version.
   *_genomic_gaps.txt.gz
       Tab-delimited text file reporting the coordinates of all gaps in the 
       top-level genomic sequences. The gaps reported include gaps specified in
       the AGP files, gaps annotated on the component sequences, and any other 
       run of 10 or more Ns in the sequences. See the "Description of files" 
       section below for details of the file format.
   *_protein.faa.gz file
       FASTA format sequences of the accessioned protein products annotated on
       the genome assembly. The FASTA title is formatted as sequence 
       accession.version plus description.
   *_protein.gpff.gz file
       GenPept format of the accessioned protein products annotated on the 
       genome assembly
   *_rm.out.gz file
       RepeatMasker output; 
       Provided for Eukaryotes 
   *_rm.run file
       Documentation of the RepeatMasker version, parameters, and library; 
       Provided for Eukaryotes 
   *_rna.fna.gz file
       FASTA format of accessioned RNA products annotated on the genome 
       assembly; Provided for RefSeq assemblies as relevant (Note, RNA and mRNA 
       products are not instantiated as a separate accessioned record in GenBank
       but are provided for some RefSeq genomes, most notably the eukaryotes.)
       The FASTA title is provided as sequence accession.version plus 
       description.
   *_rna.gbff.gz file
       GenBank flat file format of RNA products annotated on the genome 
       assembly; Provided for RefSeq assemblies as relevant
   *_rna_from_genomic.fna.gz
       FASTA format of the nucleotide sequences corresponding to all RNA 
       features annotated on the assembly, based on the genome sequence. See 
       the "Description of files" section below for details of the file format.
   ***_translated_cds.faa**: FASTA sequences of individual CDS features annotated on the genomic 
 records, conceptually translated into protein sequence. The sequence 
       corresponds to the translation of the nucleotide sequence provided in the
       *_cds_from_genomic.fna.gz file. 
  **md5checksums.txt file**: file checksums are provided for all data files in the directory

Marbouty_genbank_and_refseq_proteins: The data in Marbouty_genbank_and_refseq_proteins is a collection of protein sequences for each genome in the OMM12 communities. These protein sequences have been taken from both GenBank and RefSeq databases. For example, the file name “muribaculum_intestinales_YL27_protein_refseq.faa” has protein sequences for muribaculum intestinales YL27 as given by refseq database.
3. Marbouty_genbank_and_refseq_proteins_integrated: Contains the integrated protein sequences for each genome, taken from genbank and refseq.
4. Marbouty_annotations_for_integrated_proteins: Contains the results for the annotations done for the integrated protein sequences from ‘Data\Marbouty_genbank_and_refseq_proteins_integrated’

# Annotation results

We used different softwares such as OperonMapper, KOFAMSCAN, dbCAN2, PULs, RGI and eggNOG to get proper annotations for integrated proteins





