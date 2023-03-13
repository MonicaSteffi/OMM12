
gnmFile=read.csv(file = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/GenomeID_to_Organism_accession.csv",sep = ";",header = T)
eggNOG_GOList=lapply(list.files(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/eggNOG_2022_new/eggNOGDBRes/",full.names = T),
       function(x){
                    x_GenomeNm=str_replace_all(string = x,pattern = ".*eggNOGDBRes/",replacement = "");
                    x_GenomeNm=str_replace_all(string = x_GenomeNm,pattern = "_eggNOG.*",replacement = "");
                    
                    eggNOG_Cleaned=readRDS(x);
                    eggNOG_GO=as.character(eggNOG_Cleaned$GO);
                    names(eggNOG_GO)=eggNOG_Cleaned$query;
                    return(eggNOG_GO)
                    
                  })
x_GenomeNm=str_replace_all(string = list.files(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/eggNOG_2022_new/eggNOGDBRes/",full.names = T),pattern = ".*eggNOGDBRes/",replacement = "");
x_GenomeNm=str_replace_all(string = x_GenomeNm,pattern = "_eggNOG.*",replacement = "");

names(eggNOG_GOList)=x_GenomeNm

goList=unique(unlist(str_split(unname(unlist(eggNOG_GOList)),pattern = ",")))


go2chebi=function(goList)
{
  # install.packages("rjson")
  library("rjson")
  library("stringr")
  
  rltn=vector()
  idtn=vector()
  termtn=vector()
  
  
  for(j2 in 1: length(goList))
  {
    message(str_c(j2," is started!"))
    
    urll=str_c("https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/",goList[j2],"/complete")
    
    # hmDes="/dss/dssfs02/lwp-dss-0001/u7b03/u7b03-dss-0000/ra72seh/DB/"
    # DstnFl=str_c(hmDes,goList[j2])
    # json_data <- tryCatch(fromJSON(file=DstnFl),error=function(x){"###"})
    # file.remove(DstnFl)
    
    json_data=tryCatch(fromJSON(file = str_c("https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/",goList[j2],"/complete")),error=function(x){"###"})
    
    if(json_data!="###")
    {
      
      res <- unlist(json_data)
      
      rltn[j2]=str_c(unique(unname(unlist(res[ grepl("xRelations.*relation", names(res)) ]))),collapse = ";")
      idtn[j2]=str_c(unique(unname(unlist(res[ grepl("xRelations.*id", names(res)) ]))),collapse = ";")
      termtn[j2]=str_c(unique(unname(unlist(res[ grepl("xRelations.*term", names(res)) ]))),collapse = ";")
      
    }
    
    
    message(str_c(j2," is ended!"))
  }
  
  
  # rltnVec=str_c(rltn,collapse = ";")
  # idtnVec=str_c(idtn,collapse = ";")
  # termtnVec=str_c(termtn,collapse = ";")
  GOdf=data.frame(GOID=goList,term=termtn,Relation=rltn,CheBI=idtn)
  return(GOdf)
}
OverallGoDf=go2chebi(goList = goList)
saveRDS(file = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/godb2chebi.rds",object = OverallGoDf)
