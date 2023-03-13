joint_Marbouty_Operon=function()
{
  operonPath="S:/Abilash/Marbouty_annotation_refseq_genbank 11_2021/Data/Marbouty_annotations_for_integrated_proteins/OperonMapper_new/"
  list_files_operon=list.files(path = operonPath,pattern = "list_of_operons",all.files = T,full.names = T,recursive = T)

  for(j in 1: length(list_files_operon))
  {
    
    operon_df=read.csv(file = list_files_operon[j],header = T,sep = "\t")
    operon_id_list=which(operon_df$Operon!="")
    
    for(i in 1:length(operon_id_list))
    {
      operon_df$Operon[operon_id_list[i]:nrow(operon_df)]=i
      
    }
    
    to_be_removed=intersect(which(is.na(operon_df$PosLeft)),which(is.na(operon_df$postRight)))
    operon_df=operon_df[-(to_be_removed),]
    operon_df$ID=str_c(gnmNms[i2],"_",operon_df$PosLeft,"_",operon_df$postRight,"_",operon_df$Strand)
    operon_df$Operon=str_c(gnmNms[i2],"_",operon_df$Operon)
    colnames(operon_df)[grep(pattern = "IdGene",x=colnames(operon_df))]="GeneID"
    operon_df$GeneID=str_replace_all(string = operon_df$GeneID,pattern = "cds-",replacement = "")
    
    OperonfileName=str_c(operonPath,gnmNms[i2],"_operonList_cleaned")
    write.table(file = OperonfileName,x = operon_df,quote = F,sep = "\t",row.names = F)
    saveRDS(operon_df,file = str_c(OperonfileName,".rds"))
    
    
  }
  
  
  
  }

joint_Marbouty_Operon()