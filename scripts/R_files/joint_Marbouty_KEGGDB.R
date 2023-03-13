joint_Marbouty_KEGGDB=function()
{
  
library("stringr")
library("KEGGREST")
library(data.table)


rxnDB="S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/rxn_from_KEGG/"

if(!dir.exists(rxnDB))
{
  dir.create(rxnDB)
  
  
}
ko4eggNOG_unqKO_Genome=lapply(list.files("S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/eggNOG_2022_new/kofamscan_for_eggNOG_pred_proteins/",pattern = "_kofam_results_eval1e03.rds",full.names = T),function(x){
  
  xx=readRDS(x);
  unqKO=unique(xx$KO)
  return(unqKO)
  
})

NovelFamDB_files=list.files(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/eggNOG_2022_new/NovelFamDBRes/",pattern = ".rds",full.names = T)
ko4novFam_unqKO_Genome=lapply(NovelFamDB_files,function(x){
  
  xx=readRDS(x);
  unqKO=unique(xx$KO)
  return(unqKO)
  
})

kofam_for_operon_files=list.files("S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/OperonMapper_proteins_annotations/kofamscan/",pattern = "_kofam_results_eval1e03.rds",full.names = T)

ko4Operon_unqKO_Genome=lapply(kofam_for_operon_files,function(x){
  
  xx=readRDS(x);
  unqKO=unique(xx$KO)
  return(unqKO)
  
})

kofamscan_files=list.files(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/KOFAMSCAN_results_new/",pattern = "_kofam_results_eval1e03.rds",full.names = T)

ko4gbk_rfsq_unqKO_Genome=lapply(kofamscan_files,function(x){
  
  xx=readRDS(x);
  unqKO=unique(xx$KO)
  return(unqKO)
  
})


overallKO=unique(unlist(c(ko4gbk_rfsq_unqKO_Genome,
                          ko4Operon_unqKO_Genome,
                          ko4novFam_unqKO_Genome,
                          ko4eggNOG_unqKO_Genome)))


library("KEGGREST")


overallKO=overallKO[!is.na(overallKO)]
overallKO=overallKO[overallKO!=""]
brite_vec=rep("",length(overallKO))
rxn_vec=rep("",length(overallKO))
module_vec=rep("",length(overallKO))
module_id_vec=rep("",length(overallKO))


for(i in 1: length(overallKO))
{
  print(str_c(i, " started!"))
  keggGtOP=tryCatch(keggGet(dbentries = overallKO[i]),error=function(x){NA})
  
  print(str_c(i, " BRITE started!"))
  
  if("BRITE" %in% names(keggGtOP[[1]]))  
  {
brite_vec[i]    = str_c(keggGtOP[[1]]$BRITE,collapse = "////")
    
    
  }
  print(str_c(i, " BRITE finished!"))
  
  print(str_c(i, " MODULE started!"))
  
  if("DBLINKS" %in% names(keggGtOP[[1]]))
  {
    
   rxn1=keggGtOP[[1]]$DBLINKS[grep(keggGtOP[[1]]$DBLINKS,pattern = "RN:")]
   rxn1=unlist(str_split(str_replace_all(string = rxn1,pattern = "RN: ",replacement = ""),pattern = " "))
   
   if(length(rxn1)>0)
   {
     lll=lapply(rxn1,function(x){
       
       xx=tryCatch(keggGet(x),error=function(x){NA});
       rxnFname=str_c(rxnDB,x,"_annotation.rds")
       if(!file.exists(rxnFname))
       {
         saveRDS(file = rxnFname,object = xx)
         
       }
       return(xx)
       
       
     })  
     names(lll)=rxn1
     rxn_vec[i]=str_c(rxn1,collapse = "////")
      }
   
   
   }
  print(str_c(i, " RXNS finished!"))
  
  if("MODULE" %in% names(keggGtOP[[1]]))  
  {
    module_vec[i]    = str_c(keggGtOP[[1]]$MODULE,collapse = "////")
    module_id_vec[i]    = str_c(names(keggGtOP[[1]]$MODULE),collapse = "////")
    
  }
  print(str_c(i, " MODULE finished!"))
  
  
}

OverallKODf=data.frame(KO=overallKO,reactions=rxn_vec,brite=brite_vec,module=module_vec,module_id=module_id_vec)
saveRDS(OverallKODf,str_c(rxnDB,"OverallKODf.rds"))


## rxnID to annotations data.frame:
rxn_list_files=list.files(path = rxnDB,pattern = "_annotation.rds",full.names = T)
enzVec=rep("",length(rxn_list_files))
rclassVec=rep("",length(rxn_list_files))
commentVec=rep("",length(rxn_list_files))
eqVec=rep("",length(rxn_list_files))
defVec=rep("",length(rxn_list_files))
NmVec=rep("",length(rxn_list_files))


for(i in 1: length(rxn_list_files))
{
  rxnList=readRDS(file = rxn_list_files[i])[[1]]
  
  enzVec[i]=tryCatch(str_c(rxnList$ENZYME,collapse = "////"), error=function(x){""})
  rclassVec[i]=tryCatch(str_c(rxnList$RCLASS,collapse = "////"), error=function(x){""})
  commentVec[i]=tryCatch(str_c(rxnList$COMMENT,collapse = "////"), error=function(x){""})
  eqVec[i]=tryCatch(str_c(rxnList$EQUATION,collapse = "////"), error=function(x){""})
  defVec[i]=tryCatch(str_c(rxnList$DEFINITION,collapse = "////"), error=function(x){""})
  NmVec[i]=tryCatch(str_c(rxnList$NAME,collapse = "////"), error=function(x){""})
  
  
}

rxnID=str_replace_all(string = rxn_list_files,pattern = ".*/",replacement = "")
rxnID=str_replace_all(string = rxnID,pattern = "_annotation.rds",replacement = "")

rxnDf=data.frame(ReactionID=rxnID,Enzyme=enzVec,RClass=rclassVec,Comment=commentVec,Equation=eqVec,Definition=defVec,Name=NmVec)
defList=lapply(rxnDf$Definition,function(x){unlist(str_split(string = x,pattern = "<=>"))})
EqnList=lapply(rxnDf$Equation,function(x){unlist(str_split(string = x,pattern = "<=>"))})

SubIDEdit=function(x)
{
  xx1=str_trim(unlist(str_split(string = unlist(x)[1],pattern = "\\+")),side = "both")
  xx2=str_replace_all(string = xx1,pattern = "[0-9]+ ",replacement = "")
  return(str_c(xx2,collapse = "////"))  
  
}

PrdIDEdit=function(x)
{
  xx1=str_trim(unlist(str_split(string = unlist(x)[2],pattern = "\\+")),side = "both")
  xx2=str_replace_all(string = xx1,pattern = "[0-9]+ ",replacement = "")
  return(str_c(xx2,collapse = "////"))  
  
}


substrateNameList=lapply(defList,SubIDEdit)

ProductNameList=lapply(defList,PrdIDEdit)


substrateIDList=lapply(EqnList,SubIDEdit)

ProductIDList=lapply(EqnList,PrdIDEdit)



rxnDf$SubstrateName=unlist(substrateNameList)
rxnDf$SubstrateID=unlist(substrateIDList)
rxnDf$ProductName=unlist(ProductNameList)
rxnDf$ProductID=unlist(ProductIDList)
saveRDS(rxnDf,file = str_c(rxnDB,"rxnDf.rds" ))

##  Compound ID hierarchy
library("KEGGREST")
unqCompounds= unique(str_replace_all(str_trim(unlist(str_split(unlist(str_split(string = rxnDf$Equation,pattern = "<=>")),pattern = "\\+")),side = "both"),pattern = "[0-9]+ ",replacement = ""))


CompoundBrite=rep("",length(unqCompounds))

for(i in 1: length(unqCompounds))
{
  print(str_c(i, " started!"))
  cmpdList=tryCatch(keggGet(unqCompounds[i]),error=function(x){NA})
  if(!is.na(cmpdList))
  {
    tmpBrite=tryCatch(cmpdList[[1]]$BRITE,error=function(x){NA})
    if(!is.null(tmpBrite))
    {
      spVec=unlist(lapply(tmpBrite,function(x)
      {
        
        uu=unlist(str_split(x,pattern = ""))
        uu1=rep(0,length(uu))
        uu1[grep(pattern = " ",x=uu)]=1
        length(unlist(str_split(unlist(str_split(string = as.character(str_c(uu1,collapse = "")),pattern = "0"))[1],pattern = "")))
      }))
      
      for(j in 1: length(tmpBrite))
      {
        
        if(spVec[j]==0)
        {
          tmpBrite[j]=str_c("::>",tmpBrite[j])
        } else {
          
          prfx1=str_c(c(rep("=",spVec[j]),">"),collapse = "")
          tmpBrite[j]=str_c(prfx1,str_trim(tmpBrite[j],side = "both"))
          
        }
        
      }
      
      
      CompoundBrite[i]=str_c(unlist(tmpBrite),collapse = "////")
      
      
    }
    
  }
  
  print(str_c(i, " finished!"))
  
}

names(CompoundBrite)=unqCompounds
saveRDS(CompoundBrite,str_c(rxnDB,"CompoundBrite.rds"))

mainVec=vector()
fstVec=vector()
secondVec=vector()
th3Vec=vector()
f4Vec=vector()
k3=1
compoundnm=vector()
for(i3 in 1: length(CompoundBrite))
{
  print(str_c(i3, " started!"))
  tmpList=unlist(str_split(string = CompoundBrite[i3],pattern = "::>"))
  tmpList=tmpList[tmpList!=""]
  
  if(length(tmpList) >0)
  {
    
  
  for(j3 in 1: length(tmpList))  
  {
  catList=unlist(str_split(string = tmpList[j3],pattern = "////")  )
  mainVec[k3]=catList[1]
  fstVec[k3]=catList[2]
  secondVec[k3]=catList[3]
  th3Vec[k3]=catList[4]
  f4Vec[k3]=catList[5]
  compoundnm[k3]=names(CompoundBrite)[i3]
  k3=k3+1
  }
  print(str_c(i3, " finished!"))
  }
}
CompoundDf=data.frame(CompoundID=compoundnm,Primary_category=mainVec,Sec_category=fstVec,Third_category=secondVec,Fourth_category=th3Vec,Fifth_category=f4Vec)
CompoundDf$Sec_category=str_replace_all(string = CompoundDf$Sec_category,pattern = "=>",replacement = "")
CompoundDf$Third_category=str_replace_all(string = CompoundDf$Third_category,pattern = "==>",replacement = "")
CompoundDf$Fourth_category=str_replace_all(string = CompoundDf$Fourth_category,pattern = "===>",replacement = "")
CompoundDf$Fifth_category=str_replace_all(string = CompoundDf$Fifth_category,pattern = "====>",replacement = "")
saveRDS(CompoundDf,str_c(rxnDB,"CompoundDf.rds"))
}

joint_Marbouty_KEGGDB()