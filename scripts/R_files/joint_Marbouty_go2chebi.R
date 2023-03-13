CheBIUniRDS=list.files(path = "S:/Abilash/Marbouty_annotation_refseq_genbank 11_2021/Data/Marbouty_annotations_for_integrated_proteins/",pattern = "_FinalAnnWithCheBIUniprotFin.rds",full.names = T)
goListOverall=vector()

for(h2 in 1: length(CheBIUniRDS))
{
  FinalAnnv2=readRDS(CheBIUniRDS[h2])
  for(ii in 1:nrow(FinalAnnv2))
  {
    goList=unlist(str_split(FinalAnnv2$GOID[ii],pattern = ";"))
    goList=goList[goList!=""]
    
    goListOverall=unique(c(goListOverall,goList))
  }
  
}



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
OverallGoDf=go2chebi(goList = goListOverall)

View(FinalAnnv2)
flNm=str_replace_all(string = CheBIUniRDS,pattern = ".*/",replacement = "")
flNm=str_replace_all(string = flNm,pattern = "_FinalAnn.*",replacement = "")
baseDir="S:/Abilash/Marbouty_annotation_refseq_genbank 11_2021/Data/Marbouty_annotations_for_integrated_proteins/"

for(i in 1: length(CheBIUniRDS))
{
  FinalAnnv2=readRDS(CheBIUniRDS[i])
  GODf=OverallGoDf
  GID=vector()
  for(j in 1: nrow(GODf))
  {
    tmpIdx=grep(pattern = GODf$GOID[j],x = FinalAnnv2$GOID)
    GID[j]=str_c(unique(FinalAnnv2$GeneID.x[tmpIdx]),collapse = ";")
  }
  GODf$GeneID=GID
  flName=str_c(baseDir,flNm[i],"UniprotGO2CheBI.rds")
  saveRDS(object = GODf,file = flName)
}

rdsfile=list.files(path = baseDir,pattern = "UniprotGO2CheBI.rds",full.names = T)
flNm=str_replace_all(string = rdsfile,pattern = ".*/",replacement = "")
flNm=str_replace_all(string = flNm,pattern = "UniprotGO2CheBI.rds",replacement = "")

for(i in 1: length(rdsfile))
{
  
  rds1=readRDS(file = rdsfile[i])
  rds1=rds1[(rds1$GeneID!=""),]  
  
  write.table(x = rds1,file = str_c(baseDir,flNm[i],"UniprotGO2CheBI.csv"),quote=F,sep = "\t",row.names = F)  
}
