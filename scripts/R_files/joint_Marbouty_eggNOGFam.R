
################## Create GeneID -> Start -> Stop -> Strand -> ID:
library("seqinr")
library("stringr")

eggNOGFolders=list.dirs(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/eggNOG_2022_new/eggNOGDB/")
eggNOGFolders=eggNOGFolders[-1]
eggNOGFolders=str_replace_all(string = eggNOGFolders,pattern = "//",replacement = "/")

gbk2gnm=read.csv(file = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/GenbankID2GenomeName.txt",header = T,sep = "\t")
GenomeNames=str_replace_all(string = eggNOGFolders,pattern = ".*eggNOGDB/",replacement = "")
GenomeNames=str_replace_all(string = GenomeNames,pattern = "_eggNOGDB",replacement = "")
GenomeNames=as.character(gbk2gnm$GenomeName)[match(GenomeNames,as.character(gbk2gnm$GenomeID))]

eggNOGResFolder="S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/eggNOG_2022_new/eggNOGDBRes/"
if(!dir.exists(eggNOGResFolder))
{
  dir.create(eggNOGResFolder)
  
}

## For NovelFAM:

for(i in 1: length(eggNOGFolders))
{
  
  fastaFile=read.fasta(file = str_c(eggNOGFolders[i],"/out.emapper.genepred.fasta"),seqtype = "AA",whole.header = T)
  
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
  
  readFile=list.files(path = eggNOGFolders[i],pattern = ".annotations$",full.names = T)
  emapper_annotations=tryCatch(read.csv( readFile,header = F,sep = "\t",comment.char = "#"),error=function(x){NA})
  
  if(!is.na(emapper_annotations))
  {
    colnames(emapper_annotations)=c("query",	"seed_ortholog",	"evalue",	"score",	"eggNOG_OGs",	"max_annot_lvl",	"COG_category",	"Description",	"Preferred_name",	"GOs",	"EC",	"KEGG_ko",	"KEGG_Pathway",	"KEGG_Module",	"KEGG_Reaction",	"KEGG_rclass",	"BRITE",	"KEGG_TC",	"CAZy",	"BiGG_Reaction",	"PFAMs")  
    emapper_annotations$str_Crd=geneDf$str_Crd[match(emapper_annotations$query,as.character(geneDf$GeneID))]
    
    
    
    saveRDS(emapper_annotations,str_c(eggNOGResFolder,GenomeNames[i],"_eggNOGDB_Cleaned.rds"))
  } 
  
  
}

##########################

