joint_Marbouty_AllAnnotations_integrate=function()
{
  
#### Rename "Escherichia_coli_Mt1B1" to "escherichia_coli_strain_Mt1B1".
######################
library("stringr")
library("dplyr")
library("data.table")
  
list_files=list.files(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/",pattern = "Escherichia_coli_Mt1B1_",full.names = T,recursive = T)
fl2=gsub("Escherichia_coli_Mt1B1_", "escherichia_coli_strain_Mt1B1_", list_files)
file.rename(list_files,fl2)

list_files=list.files(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/",pattern = "Escherichia_coli_strain_Mt1B1_",full.names = T,recursive = T)
fl2=gsub("Escherichia_coli_strain_Mt1B1_", "escherichia_coli_strain_Mt1B1_", list_files)
file.rename(list_files,fl2)

list_files=list.files(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/",pattern = "muribaculum_intestinale_YL27",full.names = T,recursive = T)
fl2=gsub("muribaculum_intestinale_YL27", "muribaculum_intestinales_YL27", list_files)
file.rename(list_files,fl2)

list_dirs=list.dirs(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/", full.names = T,recursive = T)
fl2=gsub("Escherichia_coli_strain_Mt1B1", "escherichia_coli_strain_Mt1B1", list_dirs)
file.rename(list_dirs,fl2)

list_files=list.files(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/",pattern = "Akkermansia_muciniphila_YL44",full.names = T,recursive = T)
fl2=gsub("Akkermansia_muciniphila_YL44", "akkermansia_muciniphila_YL44", list_files)
file.rename(list_files,fl2)

list_files=list.files(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/",pattern = "Escherichia_coli_Mt1B1_",full.names = T,recursive = T)
fl2=gsub("Escherichia_coli_Mt1B1_", "escherichia_coli_strain_Mt1B1_", list_files)
file.rename(list_files,fl2)

list_files=list.files(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/",pattern = "Escherichia_coli_strain_Mt1B1_",full.names = T,recursive = T)
fl2=gsub("Escherichia_coli_strain_Mt1B1_", "escherichia_coli_strain_Mt1B1_", list_files)
file.rename(list_files,fl2)

list_files=list.files(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/",pattern = "Akkermansia_muciniphila_YL44",full.names = T,recursive = T)
fl2=gsub("Akkermansia_muciniphila_YL44", "akkermansia_muciniphila_YL44", list_files)
file.rename(list_files,fl2)


list_files=list.files(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/",pattern = "Escherichia_coli_strain_Mt1B1_",full.names = T,recursive = T)
fl2=gsub("Escherichia_coli_strain_Mt1B1_", "escherichia_coli_strain_Mt1B1_", list_files)
file.rename(list_files,fl2)


list_dirs=list.dirs(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/", full.names = T,recursive = T)
fl2=gsub("Escherichia_coli_strain_Mt1B1", "escherichia_coli_strain_Mt1B1", list_dirs)
file.rename(list_dirs,fl2)

list_dirs=list.dirs(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/", full.names = T,recursive = T)
fl2=gsub("Escherichia_coli_strain_Mt1B1", "escherichia_coli_strain_Mt1B1", list_dirs)
file.rename(list_dirs,fl2)

list_dirs=list.dirs(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/", full.names = T,recursive = T)
fl2=gsub("Akkermansia", "akkermansia", list_dirs)
file.rename(list_dirs,fl2)

list_dirs=list.dirs(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/", full.names = T,recursive = T)
fl2=gsub("muribaculum_intestinale_YL27", "muribaculum_intestinales_YL27", list_dirs)
file.rename(list_dirs,fl2)


##### Cleaned files.
dbCAN2_files=list.files(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/",pattern = "_overall_cgc.rds",full.names = T)
eggNOGDB_files=list.files(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/eggNOG_2022_new/eggNOGDBRes/",pattern = "_eggNOGDB_Cleaned.rds$",full.names = T)
NovelFamDB_files=list.files(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/eggNOG_2022_new/NovelFamDBRes/",pattern = "_novelFamCleaned.rds",full.names = T)
gnmFile=read.csv(file = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/GenomeID_to_Organism_accession.csv",sep = ";",header = T)

gff_files=list.files(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/",pattern = "_gff_results.txt.rds",full.names = T)

kofamscan_files=list.files(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/KOFAMSCAN_results_new/",pattern = "_kofam_results_eval1e03.rds",full.names = T)


operonMapper_files=list.files(pattern = "_operonList_cleaned.rds",path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/",full.names = T)
PULs_files=list.files(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/PULs_new/",pattern = "pul.tsv$",full.names = T)
RGI_files=list.files(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/RGI_5.1.1/RGI_for_genome/",pattern = "_genome_ann.txt",full.names = T)

kofam_for_eggnog_files=list.files("S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/eggNOG_2022_new/kofamscan_for_eggNOG_pred_proteins/",pattern = "_kofam_results_eval1e03.rds",full.names = T)  
kofam_for_operon_files=list.files("S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/OperonMapper_proteins_annotations/kofamscan/",pattern = "_kofam_results_eval1e03.rds",full.names = T)
eggNOG_for_operon_files=list.files("S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/OperonMapper_proteins_annotations/eggNOGDB/",pattern = "annotations$",full.names = T,recursive = T)
novFam_for_operon_files=list.files("S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/OperonMapper_proteins_annotations/NovelFamDB/",pattern = "annotations$",full.names = T,recursive = T)

# faa_files=list.files(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_genbank_and_refseq_proteins_integrated/",pattern = ".faa",full.names = T)
#### 
eggNOG_for_integrated_proteins=read.csv(file = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/eggNOG_integrated_proteins/out.emapper (2).annotations",header = F,sep = "\t",comment.char = "#")
eggNOG_for_integrated_proteins=eggNOG_for_integrated_proteins[,c(1,5)]
colnames(eggNOG_for_integrated_proteins)=c("query","COGgene")
eggNOG_for_integrated_proteins$COGgene=str_replace_all(string = eggNOG_for_integrated_proteins$COGgene,pattern = "@.*",replacement = "")
tmpSplit=str_split(string = eggNOG_for_integrated_proteins$query,pattern = "_")
tmpGenes=unlist(lapply(tmpSplit,function(x){str_c(x[(length(x)-1):length(x)],collapse = "_")}))
tmpGenome=unlist(lapply(tmpSplit,function(x){str_c(x[1:(length(x)-2)],collapse = "_")}))

eggNOG_for_integrated_proteins$locus_tag=tmpGenes
eggNOG_for_integrated_proteins$Genome=tmpGenome

#######################################


eggNOG_NovFam_gff_files=list.files(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/eggNOG_2022_new/",pattern = "genepred.gff$",full.names = T,recursive = T)
eggNOG_gff_files=eggNOG_NovFam_gff_files[grep(pattern = "eggNOGDB",x=eggNOG_NovFam_gff_files)]
eggNOG_genbank=str_replace_all(string = eggNOG_gff_files,pattern = ".*GC",replacement = "GC")
eggNOG_genbank=str_replace_all(string = eggNOG_genbank,pattern = "_eggNOGDB/.*",replacement = "")

eggNOG_genome=as.character(gnmFile$Genome_compatible)[match(eggNOG_genbank,as.character(gnmFile$Genbank_accession))]

eggNOG_gff=lapply(eggNOG_gff_files,function(x){
  
  rr=read.csv(file = x,header = F,sep = "\t",comment.char = "#");
  
  })


for(e1 in 1: length(eggNOG_gff))
{
  rr=eggNOG_gff[[e1]];
  rr$str_Crd=apply(rr[,c(4,5,7)],1,function(x){
    str_c(str_trim(x,side = "both"),collapse = "_");
  })
  
  rr$eggNOG_GeneID=apply(rr[,c(1,9)],1,function(x){
    
    x1=x[1];
    x2=str_replace_all(string = x[2],pattern = ";.*",replacement = "");
    x2=str_replace_all(string = x2,pattern = ".*_",replacement = "");
    str_c(x1,"_",x2)
  })
  
  eggNOG_gff[[e1]]=rr
}

names(eggNOG_gff)=eggNOG_genome


NovelFam_gff_files=eggNOG_NovFam_gff_files[grep(pattern = "NovelFamDB",x=eggNOG_NovFam_gff_files)]
NovelFam_genbank=str_replace_all(string = NovelFam_gff_files,pattern = ".*GC",replacement = "GC")
NovelFam_genbank=str_replace_all(string = NovelFam_genbank,pattern = "_NovelFamDB/.*",replacement = "")

NovelFam_genome=as.character(gnmFile$Genome_compatible)[match(NovelFam_genbank,as.character(gnmFile$Genbank_accession))]

NovelFam_gff=lapply(NovelFam_gff_files,function(x){
  
  rr=read.csv(file = x,header = F,sep = "\t",comment.char = "#");
  
})


for(e1 in 1: length(NovelFam_gff))
{
  rr=NovelFam_gff[[e1]];
  rr$str_Crd=apply(rr[,c(4,5,7)],1,function(x){
    str_c(str_trim(x,side = "both"),collapse = "_");
  })
  
  rr$NovelFam_GeneID=apply(rr[,c(1,9)],1,function(x){
    
    x1=x[1];
    x2=str_replace_all(string = x[2],pattern = ";.*",replacement = "");
    x2=str_replace_all(string = x2,pattern = ".*_",replacement = "");
    str_c(x1,"_",x2)
  })
  
  NovelFam_gff[[e1]]=rr
}

names(NovelFam_gff)=NovelFam_genome

eggNOG_NovFam_gff_files=list.files(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/eggNOG_2022_new/",pattern = "genepred.gff$",full.names = T,recursive = T)
eggNOG_gff_files=eggNOG_NovFam_gff_files[grep(pattern = "eggNOGDB",x=eggNOG_NovFam_gff_files)]
eggNOG_location="S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/eggNOG_2022_new/eggNOGDBRes/"
eggNOGDB_gnm=str_replace_all(string = eggNOGDB_files,pattern = ".*/",replacement = "")  
eggNOGDB_gnm=str_replace_all(string = eggNOGDB_gnm,pattern = "_eggNOG.*",replacement = "")  


for(e1 in 1: length(eggNOG_gff))
{
  
  eggNOG_cleaned=readRDS(file = eggNOGDB_files[e1])
  colnames(eggNOG_cleaned)[which(colnames(eggNOG_cleaned)=="evalue")]="eggNOG_evalue"
  colnames(eggNOG_cleaned)[which(colnames(eggNOG_cleaned)=="score")]="eggNOG_score"
  colnames(eggNOG_cleaned)[which(colnames(eggNOG_cleaned)=="query")]="eggNOG_query"
  
  
  kofam4eggNOG=readRDS(file = kofam_for_eggnog_files[grep(eggNOGDB_gnm[e1],x = kofam_for_eggnog_files)])
  colnames(kofam4eggNOG)=str_c("ko4eggNOG_",colnames(kofam4eggNOG))
  colnames(kofam4eggNOG)[2]="eggNOG_query"
  
  
  
  eggNOGDBRes=merge(eggNOG_cleaned,kofam4eggNOG,by="eggNOG_query",all=TRUE)
  
  kegg_brite_cols=grep(pattern = "KEGG_|BRITE",x=colnames(eggNOGDBRes))
  if(length(kegg_brite_cols) > 0)
  {
    eggNOGDBRes=eggNOGDBRes[,-(kegg_brite_cols)]
    
  }
  
  rr=eggNOG_gff[[eggNOGDB_gnm[e1]]];
  eggNOGDBRes$str_Crd=rr$str_Crd[match(as.character(eggNOGDBRes$eggNOG_query),as.character(rr$eggNOG_GeneID))]
  
  
  
  saveRDS(file = str_c(eggNOG_location,eggNOGDB_gnm[e1],"_eggNOGDB_Cleaned_with_coordinates.rds"),object = eggNOG_cleaned)
  saveRDS(file = str_c(eggNOG_location,eggNOGDB_gnm[e1],"_eggNOGDB_Cleaned_with_coordinates_with_kofam.rds"),object = eggNOGDBRes)
  
  
  }

eggNOG_NovFam_gff_files=list.files(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/eggNOG_2022_new/",pattern = "genepred.gff$",full.names = T,recursive = T)
NovelFam_gff_files=eggNOG_NovFam_gff_files[grep(pattern = "NovelFamDB",x=eggNOG_NovFam_gff_files)]
NovelFam_location="S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/eggNOG_2022_new/NovelFamDBRes/"
NovelFamDB_gnm=str_replace_all(string = NovelFamDB_files,pattern = ".*/",replacement = "")  
NovelFamDB_gnm=str_replace_all(string = NovelFamDB_gnm,pattern = "_novelFamCleaned.*",replacement = "")  

for(e1 in 1: length(NovelFamDB_gnm))
{
  rr=NovelFam_gff[[NovelFamDB_gnm[e1]]];
  NovelFam_cleaned=readRDS(file = NovelFamDB_files[e1])
  NovelFam_cleaned$str_Crd=rr$str_Crd[match(as.character(NovelFam_cleaned$query),as.character(rr$NovelFam_GeneID))]
  saveRDS(file = str_c(NovelFam_location,NovelFamDB_gnm[e1],"_NovelFamDB_Cleaned_with_coordinates.rds"),object = NovelFam_cleaned)
  
}


#######################################



GenomeNames=str_replace_all(string = kofamscan_files,pattern = ".*/",replacement = "")
GenomeNames=str_replace_all(string = GenomeNames,pattern = "_kofam.*",replacement = "")
eggNOGDB_crd_ko_files=list.files(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/eggNOG_2022_new/eggNOGDBRes/",full.names = T,
                              pattern="_eggNOGDB_Cleaned_with_coordinates_with_kofam.rds"  )

NovelFam_crd_files=list.files(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/eggNOG_2022_new/NovelFamDBRes/",full.names = T,
                              pattern="_NovelFamDB_Cleaned_with_coordinates.rds"  )


source("S:/Abilash/Marbouty_annotation_refseq_genbank_2022/scripts/R_files/annotation_main.R")

for(i in 1: length(GenomeNames))
{
  print(str_c(GenomeNames[i]," started!"))
  annotation_main(GenomeNamesHere = GenomeNames[i],dbCAN2_files,eggNOGDB_crd_ko_files,NovelFam_crd_files,gff_files,kofamscan_files,operonMapper_files,PULs_files,RGI_files,kofam_for_eggnog_files,kofam_for_operon_files,eggNOG_for_operon_files,novFam_for_operon_files,eggNOG_for_integrated_proteins,gnmFile)  
  print(str_c(GenomeNames[i]," finished!"))
  }


}
joint_Marbouty_AllAnnotations_integrate()
