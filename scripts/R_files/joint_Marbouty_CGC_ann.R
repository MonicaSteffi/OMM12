joint_Marbouty_CGC_ann=function()
{
  #----------------- dbCAN2 result for genome sequences ----------------------------------------------------------------
  library("stringr")
  library("RCurl")
  library("data.table")
  
  cgcPath="S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/dbCAN2_new//"
  
  list_cgc_files=list.files(path=cgcPath,recursive=T,full.names=T,pattern="cgc_standard.out.txt|CGC[0-9]+.gff.txt")
  gnmNms=str_replace_all(string = list_cgc_files,pattern = ".*//",replacement = "")
  gnmNms=str_replace_all(string = gnmNms,pattern = "/.*",replacement = "")
  unq_genomes=unique(gnmNms)
  
  ecoli_dbsub_files=list.files(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/dbCAN2_new/escherichia_coli_strain_Mt1B1/",pattern = "dbsub.out.txt",recursive = T,full.names = T)
  ecoli_dbsub_list=lapply(ecoli_dbsub_files,function(x){read.csv(file = x,header = T,sep = "\t")})
  ecoli_dbsub=do.call("rbind",ecoli_dbsub_list)
  ecoli_dbsub=ecoli_dbsub[!duplicated(ecoli_dbsub),]
  write.table(ecoli_dbsub,file ="S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/dbCAN2_new/escherichia_coli_strain_Mt1B1/dbsub.out.txt",quote=F,sep = "\t",row.names = F )
  
  for(i2 in 1:length(unq_genomes))
  {
   
    list_cgc_filesThis=list_cgc_files[which(gnmNms==unq_genomes[i2])]
    message(str_c("dbCAN is starting!!"))
    source("S:/Abilash/Marbouty_annotation_refseq_genbank_2022/scripts/R_files/joint_Marbouty_dbCANSht.R")
    overall_cgc=joint_Marbouty_dbCANSht(list_cgc_filesThis = list_cgc_filesThis,gnmNmsThis = unq_genomes[i2])
    message(str_c("dbCAN is done!!"))
    
    message(str_c("dbCAN substrate adding started!!"))
    
    dbsub.out.file=tryCatch(read.csv(file = str_c(cgcPath,unq_genomes[i2],"/dbsub.out.txt"),header = T,sep = "\t"),error=function(x){NA})
    
    if(!is.na(dbsub.out.file))
    {
      overall_cgc$dbCAN2_substrate=dbsub.out.file$Substrate[match(overall_cgc$dbCAN2_GeneID,as.character(dbsub.out.file$Gene.ID))]
      overall_cgc$dbCAN2_EC=dbsub.out.file$Subfam.EC[match(overall_cgc$dbCAN2_GeneID,as.character(dbsub.out.file$Gene.ID))]
      overall_cgc$dbCAN2_EC=unlist(lapply(overall_cgc$dbCAN2_EC,function(x){
        xx=unlist(str_split(string = x,pattern = "\\|"));
        xx=str_replace_all(string = xx,pattern = ":.*",replacement = "");
        return(str_c(xx,collapse = "|"))
      }))
      
    }
    
    saveRDS(overall_cgc,str_c("S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/",unq_genomes[i2],"_overall_cgc.rds"))
    
    
    
  }
}

joint_Marbouty_CGC_ann()


