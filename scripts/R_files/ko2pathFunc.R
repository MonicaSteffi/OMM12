ko2pathFunc=function()
{
library("stringr") 
  library("KEGGREST")
rxnDB="S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/rxn_from_KEGG/"
mod2ko2gene_files=list.files(path = rxnDB,pattern = "_Mod2KO2Gene.rds",full.names = T)
GenomeNames=str_replace_all(string = mod2ko2gene_files,pattern = ".*/",replacement = "")
GenomeNames=str_replace_all(string = GenomeNames,pattern = "_Mod2.*",replacement = "")


for(m1 in 1: length(mod2ko2gene_files))
{
  
  Mod2KO2Gene=readRDS(mod2ko2gene_files[m1])
  
  unqKO=unique(Mod2KO2Gene$KO)
  pathID=list()
  pathName=list()
  k5=1
  for(m2 in 1: length(unqKO))
  {
    print(str_c(m2," started!"))
    
    koInfo=keggGet(unqKO[m2])
    if("PATHWAY" %in% names(koInfo[[1]]))
    {
      pathID[[m2]]=unlist(names(koInfo[[1]]$PATHWAY))  
      pathName[[m2]]=unlist(unname(koInfo[[1]]$PATHWAY))  
      
      
      
    } else {
      pathID[[m2]]=""
      pathName[[m2]]=""
    }
    print(str_c(m2," finished!"))
    
  }
  names(pathID)=unqKO
  names(pathName)=unqKO
  
  for(m3 in 1: length(pathID))
  {
    if(all(pathID[[m3]]!=""))
    {
      tmpIdx=which(Mod2KO2Gene$KO==names(pathID)[m3])
      pathIDs=unname(pathID[[m3]])
      pathNms=unname(pathName[[m3]])
      
      for(m4 in 1: length(pathIDs))
      {
        tmpDf=Mod2KO2Gene[tmpIdx,]
        tmpDf$Pathway_ID=pathIDs[m4]
        tmpDf$Pathway_Name=pathNms[m4]
        
        if(k5==1)
        {
          
          OverallPathDf=tmpDf
          k5=k5+1
        } else {
          
          OverallPathDf=rbind(OverallPathDf,tmpDf)
        }
        
        
      }
    }
    
  }
  
  saveRDS(OverallPathDf,str_c(rxnDB,GenomeNames[m1],"_KO2Pathway2Mod2Gene.rds"))
  
}
}

ko2pathFunc()