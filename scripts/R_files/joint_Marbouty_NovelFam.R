##### Read the nov family database outputs from https://novelfams.cgmlab.org/static/nmgfams_app/downloads/NOV-gene_family_composition.tsv.gz and https://github.com/AlvaroRodriguezDelRio/NovFamilies/raw/main/Supplementary%20tables.xlsx
library("tidyverse")
library("data.table")
library("readxl")

nov_gene_family=fread(file = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/eggNOG_2022_new/NovelFam_database/NOV-gene_family_composition.tsv/microbial_genomes-v1.clustering.folded.parsed.NoISO.3sp.codes.tsv",sep = "\t",header = F)

koList=lapply(novFam$`Table S3`$`Novel famiy codes`,function(x){unlist(str_split(string = x,pattern = ","))})
names(koList)=as.character(novFam$`Table S3`$KO)

novelFamResFolder="S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/eggNOG_2022_new/NovelFamDBRes/"
if(!dir.exists(novelFamResFolder))
{
  dir.create(novelFamResFolder)
  
}
path <- "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/eggNOG_2022_new/NovelFam_database/Supplementary tables (2).xlsx"
novFam=path %>% 
  excel_sheets() %>% 
  set_names() %>% 
  map(read_excel, path = path)

################## Create GeneID -> Start -> Stop -> Strand -> ID:
library("seqinr")

novelfamFolders=list.dirs(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/eggNOG_2022_new/NovelFamDB/")
novelfamFolders=novelfamFolders[-1]
novelfamFolders=str_replace_all(string = novelfamFolders,pattern = "//",replacement = "/")

eggNOGFolders=list.dirs(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/eggNOG_2022_new/eggNOGDB/")
eggNOGFolders=eggNOGFolders[-1]
eggNOGFolders=str_replace_all(string = eggNOGFolders,pattern = "//",replacement = "/")

gbk2gnm=read.csv(file = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/GenbankID2GenomeName.txt",header = T,sep = "\t")
GenomeNames=str_replace_all(string = novelfamFolders,pattern = ".*NovelFamDB/",replacement = "")
GenomeNames=str_replace_all(string = GenomeNames,pattern = "_NovelFamDB",replacement = "")
GenomeNames=as.character(gbk2gnm$GenomeName)[match(GenomeNames,as.character(gbk2gnm$GenomeID))]
  
## For NovelFAM:

for(i in 1: length(novelfamFolders))
{
  
  fastaFile=read.fasta(file = str_c(novelfamFolders[i],"/out.emapper.genepred.fasta"),seqtype = "AA",whole.header = T)
  
  list_info=lapply(names(fastaFile),function(x){str_trim(unlist(str_split(string = x,pattern = "#")),side = "both")})
  gnmID=unlist(lapply(list_info,function(x){x[1]}))
  strt_info=unlist(lapply(list_info,function(x){x[2]}))
  stp_info=unlist(lapply(list_info,function(x){x[3]}))
  strnd_info=unlist(lapply(list_info,function(x){x[4]}))
  gnmID2_info=unlist(lapply(list_info,function(x){x[5]}))
  geneDf=data.frame(GeneID=gnmID,Start=strt_info,Stop=stp_info,Strand=strnd_info,ID=gnmID2_info)
  geneDf$ID=str_replace_all(string = geneDf$ID,pattern = ";.*",replacement = "")
  geneDf$ID=str_replace_all(string = geneDf$ID,pattern = "ID=",replacement = "")
  geneDf$Strandv2="+"
  geneDf$Strandv2[which(geneDf$Strand=="-1")]="-"
  geneDf$str_Crd=apply(geneDf[,c("Start","Stop","Strandv2")],1,function(x){str_c(x,collapse = "_")})
  
  readFile=list.files(path = novelfamFolders[i],pattern = ".annotations$",full.names = T)
  emapper_annotations=tryCatch(read.csv( readFile,header = F,sep = "\t",comment.char = "#"),error=function(x){NA})
  
  if(!is.na(emapper_annotations))
  {
   colnames(emapper_annotations)=c("query",	"target",	"evalue",	"score",	"novel_fam")  
   emapper_annotations$str_Crd=geneDf$str_Crd[match(emapper_annotations$query,as.character(geneDf$GeneID))]
   emapper_annotations$Novel_familyID=str_replace_all(string = emapper_annotations$novel_fam,pattern = "_.*",replacement = "")
   
   koIdx=lapply(emapper_annotations$Novel_familyID,function(x){grep(pattern = x,x=koList)})
   names(koIdx)=as.character(emapper_annotations$Novel_familyID)
   
   if(any(unlist(lapply(koIdx,function(x){length(x)}))) >0) {
     koCol=lapply(koIdx,function(x){str_c(novFam$`Table S3`$KO[x],collapse = ";")})
     koDesc=lapply(koIdx,function(x){str_c(novFam$`Table S3`$`KO description`[x],collapse = ";")})
     novelFamID=names(koIdx)
     koDf=data.frame(KO=koCol,KO_Desc=koDesc,Novel_familyID=novelFamID)
     
   } else {
     
     koDf=data.frame(KO=rep("",length(koIdx)),KO_Desc=rep("",length(koIdx)),Novel_familyID=names(koIdx))
     
     
   }
   
   emapper_annotations=merge(emapper_annotations,koDf,by="Novel_familyID",all=TRUE)
   saveRDS(emapper_annotations,str_c(novelFamResFolder,GenomeNames[i],"_novelFamCleaned.rds"))
   } 
  
  
}
