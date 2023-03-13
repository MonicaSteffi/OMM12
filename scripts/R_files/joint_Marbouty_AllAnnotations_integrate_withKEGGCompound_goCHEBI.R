baseDir="S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/"
new_FinMat_files=list.files(path="S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/",pattern = "_new_FinMatShtv2.rds",full.names = T)
OverallGoDf=readRDS(file = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/godb2chebi.rds")

rxnDB="S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/rxn_from_KEGG/"
OverallKODf=readRDS(str_c(rxnDB,"OverallKODf.rds"))
CompoundBrite=readRDS(str_c(rxnDB,"CompoundBrite.rds"))
CompoundDf=readRDS(str_c(rxnDB,"CompoundDf.rds"))
rxnDf=readRDS(str_c(rxnDB,"rxnDf.rds"))

shetty_mod2koMelt=readRDS(file = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Extra/Shetty_et_al_MDbMM16-master/shetty_mod2koMelt.rds")
shettyDf=readRDS(file = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Extra/Shetty_et_al_MDbMM16-master/shetty_annotation_integrated.rds")
shetty_mod2kolist=readRDS(file = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Extra/Shetty_et_al_MDbMM16-master/shetty_mod2kolist.rds")
shetty_mod2koMat=readRDS(file = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Extra/Shetty_et_al_MDbMM16-master/shetty_mod2koMat.rds")

gnmNames=str_replace_all(string = new_FinMat_files,pattern = ".*/",replacement = "")
gnmNames=str_replace_all(string = gnmNames,pattern = "_new_FinMatShtv2.rds",replacement = "")


for(i in 1: length(gnmNames))
{
  
  allAnnotate=readRDS(str_c(baseDir,gnmNames[i],"_new_FinMatSht.rds"))
  allAnnotateNew=allAnnotate
  ### Final cleaning of the newFinMatSht object:
  for(jj in 1: ncol(allAnnotate))
  {
    xx=allAnnotate[,jj] 
    xxList=str_split(string = xx,pattern = "//")

    l4=unlist(lapply(xxList,function(x)
    {
      l1=str_replace_all(string = x,pattern = "/",replacement = "");
      l2=l1[l1!=""];
      l3=l2[!is.na(l2)];
      return(str_c(unique(l3),collapse = "///"))
    }))
    
    allAnnotateNew[,jj]=l4
      }
  allAnnotateNew$GeneIDFin=allAnnotateNew$locus_tag
  allAnnotateNew$GeneIDFin[which((allAnnotateNew$GeneIDFin==""))]=allAnnotateNew$str_Crd[which((allAnnotateNew$GeneIDFin==""))]
  
  saveRDS(file = str_c(baseDir,gnmNames[i],"_new_FinMatShtv2.rds"),object = allAnnotateNew)
  
}


for(i in 1: length(gnmNames))
{
  allAnnotateNew=readRDS(file = str_c(baseDir,gnmNames[i],"_new_FinMatShtv2.rds"))
  allAnnotateNew$KO_overall=apply(allAnnotateNew[,colnames(allAnnotateNew)[grep(pattern = "KO$",x=colnames(allAnnotateNew))]],1,function(xx){
    xx=xx[xx!=""];
    xx=xx[!is.na(xx)];
    xx=unique(xx);
    
    if(length(xx) > 0)
    {
      return(str_c(xx,collapse = "///"))
    } else {
      return("")
    }
    
    
  })
  
  saveRDS(file = str_c(baseDir,gnmNames[i],"_new_FinMatShtv2.rds"),object = allAnnotateNew)
  
  
}

briteList=list()
ReactionList=list()
ModuleIDList=list()
ModuleNmList=list()
ModuleList=list()

for(i in 1: nrow(OverallKODf))
{
 tmpBrite=unlist(str_split(as.character(OverallKODf$brite[i]),pattern = "////"))
 firstLines=grep(pattern = "^ ",x = tmpBrite,invert = TRUE)
 
 firstLineVec=rep("", length(firstLines))
 for(j in 1: length(firstLines)) 
 {
   if(j == length(firstLines))
   {
     firstLineVec[j]=str_c(tmpBrite[firstLines[j]: length(tmpBrite)],collapse = "///")
     
   } else {
 
     firstLineVec[j]=str_c(tmpBrite[firstLines[j]: (firstLines[j+1]-1)],collapse = "///")
         
   }
 }
 
 briteList[[i]]=firstLineVec
 ReactionList[[i]]=unlist(str_split(string = OverallKODf$reactions[i],pattern = "////"))
 ModuleIDList[[i]]=unlist(str_split(string = OverallKODf$module_id[i],pattern = "////"))
 ModuleNmList[[i]]=unlist(str_split(string = OverallKODf$module[i],pattern = "////"))
 
 
 
}

names(briteList)=OverallKODf$KO
names(ReactionList)=OverallKODf$KO
names(ModuleIDList)=OverallKODf$KO
names(ModuleNmList)=OverallKODf$KO
ModuleList=lapply(1:length(ModuleNmList),function(x){
  
  idList=ModuleIDList[[x]];
  NmList=ModuleNmList[[x]];

  if(length(idList) > 0 && any(idList !=""))
  {
    return(str_c(idList,": ",NmList))
    
  } else {
    
    return("")
  }
  
  
  })

names(ModuleList)=OverallKODf$KO
saveRDS(briteList,str_c(rxnDB,"briteList.rds"))
saveRDS(ReactionList,str_c(rxnDB,"ReactionList.rds"))
saveRDS(ModuleIDList,str_c(rxnDB,"ModuleIDList.rds"))
saveRDS(ModuleNmList,str_c(rxnDB,"ModuleNmList.rds"))
saveRDS(ModuleList,str_c(rxnDB,"ModuleList.rds"))

mod2koList=list()
unqMod=unique(unname(unlist(ModuleList)))
unqMod=unqMod[unqMod!=""]

for(y in 1: length(unqMod))
{
  modGrp=grep(pattern = unqMod[y],x = ModuleList)  
  
  if(length(modGrp) > 0)
  {
    mod2koList[[y]]=names(ModuleList)[modGrp]
    
  }

  }

names(mod2koList)=unqMod
saveRDS(mod2koList,str_c(rxnDB,"mod2koList.rds"))

brite2koList=list()
unqBrite=unique(unname(unlist(briteList)))
unqBrite=unqBrite[unqBrite!=""]

for(y in 1: length(unqBrite))
{
  print(str_c(y, " started!"))
  modGrp=grep(pattern = unqBrite[y],x = briteList,fixed = T)  
  
  if(length(modGrp) > 0)
  {
    brite2koList[[y]]=names(briteList)[modGrp]
    
  }
  print(str_c(y, " finished!"))
  
}

names(brite2koList)=unqBrite
saveRDS(brite2koList,str_c(rxnDB,"brite2koList.rds"))



for(j in 1: length(brite2koList))
{
  
  print(str_c(j, " started!"))
  tmpKO=unique(unlist(brite2koList[[j]]))
  tmpList=OverallKODf[(OverallKODf$KO %in% tmpKO),]
  
  if(dim(tmpList)[1] > 0)
  {
    tmpDf=tmpList
    tmpDf$Brite=names(brite2koList)[j]
  }
  
  if(j ==1)
  {
    Brite2KO2Gene=tmpDf
    
  } else {
    
    Brite2KO2Gene=rbind(Brite2KO2Gene,tmpDf)
  }
  print(str_c(j, " finished!"))
  
}


saveRDS(Brite2KO2Gene,str_c(rxnDB,"Brite2KO2Gene.rds"))


allAnnotateNewFiles=list.files(path = baseDir,pattern = "_new_FinMatShtv2.rds",full.names = T)
gnmNames=str_replace_all(string = allAnnotateNewFiles,pattern = ".*/",replacement = "")
gnmNames=str_replace_all(string = gnmNames,pattern = "_new.*",replacement = "")


for(i in 1: length(allAnnotateNewFiles))
{
  allAnnotateNew=readRDS(allAnnotateNewFiles[i])
  ko_overallList=lapply(allAnnotateNew$KO_overall,function(x){
    unlist(str_split(string = x,pattern = "///"))
    
  })
  
  names(ko_overallList)=allAnnotateNew$GeneIDFin
  
  ko2GeneList=list()
  unqKO=unique(unname(unlist(ko_overallList)))
  unqKO=unqKO[unqKO!=""]
  k1=1
  for(y in 1: length(unqKO))
  {
    koGrp=grep(pattern = unqKO[y],x = ko_overallList)  
    
    if(length(koGrp) > 0)
    {
      ko2GeneList[[y]]=names(ko_overallList)[koGrp]
      tmpKO2GDf=data.frame(GeneIDFin=ko2GeneList[[y]])
      tmpKO2GDf$KO=unqKO[y]
      
      if(k1==1)
      {
        OverallKo2GnDf=tmpKO2GDf
        k1=k1+1
      } else {
        
        OverallKo2GnDf=rbind(OverallKo2GnDf,tmpKO2GDf)
      }
      }
    
  }
  
  names(ko2GeneList)=unqKO
  
  
  for(j in 1: length(mod2koList))
  {
    tmpKO=unique(unlist(mod2koList[[j]]))
    tmpList=OverallKo2GnDf[(OverallKo2GnDf$KO %in% tmpKO),]
  
    if(dim(tmpList)[1] > 0)
    {
      tmpDf=tmpList
      tmpDf$Module=names(mod2koList)[j]
    }
    
    if(j ==1)
    {
      Mod2KO2Gene=tmpDf
      
    } else {
      
      Mod2KO2Gene=rbind(Mod2KO2Gene,tmpDf)
    }
    
  }
  saveRDS(Mod2KO2Gene,str_c(rxnDB,gnmNames[i],"_Mod2KO2Gene.rds"))
  k2=1
  for(j2 in 1: length(unqKO))
  {
     
    tmprxn=ReactionList[[unqKO[j2]]]
    tmprxn=tmprxn[tmprxn!=""]
    tmprxn=tmprxn[!is.na(tmprxn)]
  
    if(length(tmprxn) > 0)
    {
      
      tmprxnDf=rxnDf[match(tmprxn,rxnDf$ReactionID),]
      tmpGene=as.character(Mod2KO2Gene$GeneIDFin)[which(Mod2KO2Gene$KO==unqKO[j2])]
      
      if(length(tmpGene) > 0)
      {
        xx1=do.call("rbind", replicate(length(tmpGene), tmprxnDf, simplify = FALSE))
        tmpKO2Gene2rxn=cbind(unlist(lapply(tmpGene,function(x){rep(x,nrow(tmprxnDf))})),xx1)
        tmpKO2Gene2rxn$KO=unqKO[j2]
        
        if(k2 ==1 )
        {
          Overall_KO2Gene2rxn=tmpKO2Gene2rxn
          k2=k2+1
          
        }  else {
          Overall_KO2Gene2rxn=rbind(Overall_KO2Gene2rxn,tmpKO2Gene2rxn)
          
        }
      }
      
      
      }
    }
  colnames(Overall_KO2Gene2rxn)[1]="GeneIDFin"
  saveRDS(Overall_KO2Gene2rxn,str_c(rxnDB,gnmNames[i],"_Overall_KO2Gene2rxn.rds"))
  
  for(y5 in 1: length(brite2koList))
  {
    tmpKO=unique(unlist(brite2koList[[y5]]))
    tmpList=Brite2KO2Gene[(Brite2KO2Gene$KO %in% tmpKO),]
    
    if(dim(tmpList)[1] > 0)
    {
      tmpDf=tmpList
      tmpDf=tmpDf[,-which(colnames(tmpDf)=="brite")]
    }
    tmpGenes=as.character(OverallKo2GnDf$GeneIDFin)[(OverallKo2GnDf$KO %in% tmpKO)]
    
    if(length(tmpGenes) > 0)
    {
      
      xx1=do.call("rbind", replicate(length(tmpGenes), tmpDf, simplify = FALSE))
      tmpBrite2Mod2KO2Gene=cbind(unlist(lapply(tmpGenes,function(x){rep(x,nrow(tmpDf))})),xx1)
      
      
    }
    if(y5 ==1)
    {
      Brite2Mod2KO2Gene=tmpBrite2Mod2KO2Gene
      
    } else {
      
      Brite2Mod2KO2Gene=rbind(Brite2Mod2KO2Gene,tmpBrite2Mod2KO2Gene)
    }
  }
  
  colnames(Brite2Mod2KO2Gene)[1]="GeneIDFin"
  Brite2Mod2KO2Gene[c('First_category','Second_category','Third_category','Fourth_category','Fifth_category')]=str_split_fixed(string = Brite2Mod2KO2Gene$Brite,pattern = "///",n = 5)
  Brite2Mod2KO2Gene$First_category=str_trim(string = Brite2Mod2KO2Gene$First_category,side = "both")  
  Brite2Mod2KO2Gene$Second_category=str_trim(string = Brite2Mod2KO2Gene$Second_category,side = "both")  
  Brite2Mod2KO2Gene$Third_category=str_trim(string = Brite2Mod2KO2Gene$Third_category,side = "both")  
  Brite2Mod2KO2Gene$Fourth_category=str_trim(string = Brite2Mod2KO2Gene$Fourth_category,side = "both")  
  Brite2Mod2KO2Gene$Fifth_category=str_trim(string = Brite2Mod2KO2Gene$Fifth_category,side = "both")  
  
  saveRDS(Brite2Mod2KO2Gene,str_c(rxnDB,gnmNames[i],"_Brite2Mod2KO2Gene.rds"))
  
  
  
}

### Shetty module connection to Genes:
shettyMod2KOList=lapply(shettyDf$KO,function(x){
  
  unlist(str_split(string = x,pattern=","))
})
names(shettyMod2KOList)=shettyDf$module_id

OverallKODf=readRDS(str_c(rxnDB,"OverallKODf.rds"))
Fin_genomes=str_replace_all(string = allAnnotateNewFiles,pattern = ".*/",replacement = "")
Fin_genomes=str_replace_all(string = Fin_genomes,pattern = "_new.*",replacement = "")

for(j3 in 1: length(allAnnotateNewFiles))
{
  print(str_c(j3, " started!"))
  allAnnotateNew=readRDS(allAnnotateNewFiles[j3])
  allAnnotateNew_KOList=lapply(allAnnotateNew$KO_overall,function(x){
    
    unique(unlist(str_split(string = x,pattern = "///")))
  })
  names(allAnnotateNew_KOList)=allAnnotateNew$GeneIDFin
  
  for(i3 in 1: length(shettyMod2KOList))
  {
    tmpKO=unique(intersect(shettyMod2KOList[[i3]],allAnnotateNew$KO_overall))
    tmpKO=tmpKO[tmpKO!=""]
    if(length(tmpKO) > 0)
    {
      
      for(k3 in 1: length(tmpKO))
      {
        tmpIdx=which(allAnnotateNew_KOList %in% tmpKO[k3]==TRUE)
        
        tmpAnn=allAnnotateNew[tmpIdx,]
        ss=which(tmpAnn==tmpKO[k3],arr.ind=T)
        ko_col=setdiff(colnames(tmpAnn)[unique(ss[,2])],"KO_overall")[1]
        
        ko_def=unique(tmpAnn[,str_c(ko_col,"_definition")])
        ko_def=ko_def[ko_def!=""]
        
        tmpAnn=tmpAnn[,c("GeneIDFin","str_Crd","locus_tag","locus_tag_refseq","locus_tag_genbank","ProteinID","Product")]
        tmpAnn$module_id=names(shettyMod2KOList)[i3]
        tmpAnn$KO=tmpKO[k3]
        tmpAnn$KO_definition=ko_def
        
        if(k3==1)
        {
          Ov_MF=tmpAnn
          
        } else {
          
          Ov_MF=rbind(Ov_MF,tmpAnn)
        }
        
      }
      tmpshetty_Ov_MF=merge(Ov_MF,shettyDf,by="module_id")
      colnames(tmpshetty_Ov_MF)[which(colnames(tmpshetty_Ov_MF)=="KO.y")]="Overall_KO_for_Module"
    } else {
      
      col_names=c("module_id","GeneIDFin","str_Crd","locus_tag","locus_tag_refseq","locus_tag_genbank" ,"ProteinID","Product",          
                                             "KO.x","KO_definition","module_desc","Overall_KO_for_Module","module_name"  ,"category_1",        "category_2"  ,"category_3",       
                                             "note","trophic_included","trophic_level","trophic_class","trophic_guild" ,"Degradation"  )
    tmpshetty_Ov_MF <- data.frame(matrix(ncol = length(col_names), nrow = 0))
      colnames(tmpshetty_Ov_MF)=col_names
      }
    
    if(i3 == 1)
    {
      
      shetty_Ov_MF=tmpshetty_Ov_MF
    } else {
      
      shetty_Ov_MF=rbind(shetty_Ov_MF,tmpshetty_Ov_MF)
    }
    
  
  }
  
  saveRDS(shetty_Ov_MF,str_c(baseDir,Fin_genomes[j3],"_shetty_modules_with_genes_integration.rds"))
  print(str_c(j3, " finished!"))
  
}

### Reaction from KEGG connection to Genes:

brite2mod2ko2geneFiles=list.files(path = rxnDB,pattern = "_Brite2Mod2KO2Gene.rds",full.names = T)
genomeNames=str_replace_all(string = brite2mod2ko2geneFiles,pattern = ".*/",replacement = "")
genomeNames=str_replace_all(string = genomeNames,pattern = "_Brite2.*",replacement = "")


for(j4 in 1: length(genomeNames))
{
  Brite2Mod2KO2Gene = readRDS(brite2mod2ko2geneFiles[j4])
  Mod2KO2Gene=readRDS(str_c(rxnDB,genomeNames[j4],"_Mod2KO2Gene.rds"))
  Mod2KO2Gene=Mod2KO2Gene[,c("GeneIDFin","KO")]
  Mod2KO2Gene=Mod2KO2Gene[!duplicated(Mod2KO2Gene),]
  
  Brite2Mod2KO2Gene =Brite2Mod2KO2Gene[,c("KO","reactions")]
  colnames(Brite2Mod2KO2Gene)[2]="ReactionID"
  Mod2Brite=merge(Brite2Mod2KO2Gene,Mod2KO2Gene,by="KO")
  Mod2Brite=Mod2Brite[!duplicated(Mod2Brite),]
  tbr_rows=which(Mod2Brite$ReactionID=="")
  if(length(tbr_rows) > 0)
  {
    Mod2Brite=Mod2Brite[-tbr_rows,]
    
  }
  Moddf=Mod2Brite[,c(1,2)]
  Moddf=Moddf[!duplicated(Moddf),]
  
  
  rxnList=lapply(Moddf$ReactionID,function(x){
    
    unlist(str_split(string = x,pattern = "////"))
  })
  names(rxnList)=Moddf$KO
  
  for(k4 in 1: length(rxnList))
  {
    tmpRxn=rxnList[[k4]]
    tmpRxndf=rxndf[(rxndf$ReactionID %in% tmpRxn),]
    tmpRxndf$KO=names(rxnList)[k4]
    tmpDf=merge(tmpRxndf,Mod2Brite,by="KO")
    colnames(tmpDf)[grep(pattern = "ReactionID.y",x=colnames(tmpDf))]="Overall_ReactionID_for_KO"
   
    
    if(k4==1)
    {
      
      OverallRxndf=tmpDf
    } else {
      
      OverallRxndf=rbind(OverallRxndf,tmpDf)
    }
    }
  
  saveRDS(OverallRxndf,str_c(baseDir,genomeNames[j4],"_Gene2Rxn2Substrate_Product.rds"))
  }


##### GO-chebi connection to Genes:

for(j5 in 1: length(genomeNames))
{
  allAnnotateNew=readRDS(allAnnotateNewFiles[grep(pattern = genomeNames[j5],x=allAnnotateNewFiles)])
  GO_cols=grep(pattern = "GOs",x=colnames(allAnnotateNew))
  
  for(k5 in 1: length(GO_cols))
  {
    GOList=lapply(allAnnotateNew[,GO_cols[k5]],function(x){
      
      unique(unlist(str_split(string = x,pattern = ",")))
    })
  names(GOList)=allAnnotateNew$GeneIDFin  
    
  for(l5 in 1: length(GOList))
  {
    tmpGO=GOList[[l5]]
    tmpGO=tmpGO[tmpGO!=""]  
    tmpGO=tmpGO[tmpGO!="-"]
    
    if(length(tmpGO) > 0)
    {
      
    tmpGODf=OverallGoDf[(OverallGoDf$GOID %in% tmpGO),]  
    tmpGODf$GeneIDFin=names(GOList)[l5]  
    } else {
      col_names=c("GOID","term" ,"Relation", "CheBI","GeneIDFin")
      tmpGODf <- data.frame(matrix(ncol = length(col_names), nrow = 0))
      colnames(tmpGODf)=col_names
      
      
    }
    
    if(l5==1)
    {
      
      OverallGO2GeneDf=tmpGODf
    } else {
      
      OverallGO2GeneDf=rbind(OverallGO2GeneDf,tmpGODf)
      
    }
   
   
   }
    
  if(k5 ==1)
  {
    
    OV_GO_Df=OverallGO2GeneDf
  } else {
    OV_GO_Df=rbind(OV_GO_Df,OverallGO2GeneDf)
  }
    
  }
  
  OV_GO_Df=OV_GO_Df[!duplicated(OV_GO_Df),]
  saveRDS(OV_GO_Df,str_c(baseDir,genomeNames[j5],"_GO2substrate2Gene.rds"))
}


##### Final edit to allAnnotateNew:

for(j5 in 1: length(genomeNames))
{
  print(str_c(j5, " started!"))
  allAnnotateNew=readRDS(allAnnotateNewFiles[grep(pattern = genomeNames[j5],x=allAnnotateNewFiles)])
  empty_cluster_rows=which(allAnnotateNew$dbCAN2_Cluster=="")
  allAnnotateNewv2=allAnnotateNew[empty_cluster_rows,]
  to_change_allAnnotateNew=allAnnotateNew[-empty_cluster_rows,]
  
  clusterList=lapply(to_change_allAnnotateNew$dbCAN2_Cluster,function(x){
    
    unique(unlist(str_split(string = x,pattern = "///")))
  })
  names(clusterList)=as.character(to_change_allAnnotateNew$GeneIDFin)
  
  for(k5 in 1: length(clusterList))
  {
    tmpDf=to_change_allAnnotateNew[which(to_change_allAnnotateNew$GeneIDFin==names(clusterList)[k5]),]
    
    xx1=do.call("rbind", replicate(length(unlist(clusterList[[k5]])), tmpDf, simplify = FALSE))
    xx1$dbCAN2_Cluster=unlist(clusterList[[k5]])
    if(k5==1)
    {
      
      to_changeDf=xx1
    } else {
      
      to_changeDf=rbind(to_changeDf,xx1)
    }
    }
  
  allAnnotateNewFin=rbind(allAnnotateNewv2,to_changeDf)
  saveRDS(allAnnotateNewFin,str_c(baseDir,genomeNames[j5],"_allAnnotationFin.rds"))
  print(str_c(j5, " finished!"))
  
  }

##### Download Module definition and Module Hierarchy
library("KEGGREST")
library("stringr")

modList <- keggList("module")
names(modList)=str_replace_all(string = names(modList),pattern = "md:",replacement = "")
mod_def=vector()
mod_hierarchy=vector()
for(i in 1: length(names(modList)))
{
  print(str_c(i, ": ",names(modList)[i]," started!"))
  kk=keggGet(names(modList)[i])
  mod_def[i]=kk[[1]]$DEFINITION
  mod_hierarchy[i]=kk[[1]]$CLASS
  print(str_c(i, ": ",names(modList)[i]," finished!"))
  
  
}
names(mod_def)=names(modList)
names(mod_hierarchy)=names(modList)

mod2ko2gene_files=list.files(path = rxnDB,pattern = "_Mod2KO2Gene.rds",full.names = T)
modGenomes=str_replace_all(string = mod2ko2gene_files,pattern = ".*/",replacement = "")
modGenomes=str_replace_all(string = modGenomes,pattern = "_Mod2K.*",replacement = "")


for(m1 in 1: length(mod2ko2gene_files))
{
  
  Mod2KO2Gene=readRDS(mod2ko2gene_files[m1])
  Mod2KO2Gene$ModuleID=str_replace_all(string = Mod2KO2Gene$Module,pattern = ":.*",replacement = "")
  Mod2KO2Gene$Module_definition=mod_def[match(Mod2KO2Gene$ModuleID,as.character(names(mod_def)))]
  Mod2KO2Gene$Module_hierarchy=mod_hierarchy[match(Mod2KO2Gene$ModuleID,as.character(names(mod_hierarchy)))]
  Mod2KO2Gene[c('Module_category_I','Module_category_II','Module_category_III')]=str_split_fixed(string = Mod2KO2Gene$Module_hierarchy,pattern = ";",n = 3)
  
  saveRDS(file = str_c(rxnDB,modGenomes[m1],"_Mod_hierarchy2KO2Gene.rds"),object = Mod2KO2Gene )
  }


############################################### Convert rds files to csv files

list_files_baseDir=list.files(path = baseDir,pattern = "_allAnnotationFin.rds",full.names = T,recursive = F)
fileNames=str_replace_all(string = list_files_baseDir,pattern = ".rds",replacement = ".csv")

for(i in 1: length(list_files_baseDir))
{
  rr=readRDS(list_files_baseDir[i])
  write.table(rr,file = fileNames[i],sep = "\t",row.names = F,append = F,quote = F)
  
}


list_files_baseDir=list.files(path = baseDir,pattern = "_GO2substrate2Gene.rds",full.names = T,recursive = F)
fileNames=str_replace_all(string = list_files_baseDir,pattern = ".rds",replacement = ".csv")

for(i in 1: length(list_files_baseDir))
{
  rr=readRDS(list_files_baseDir[i])
  write.table(rr,file = fileNames[i],sep = "\t",row.names = F,append = F,quote = F)
  
}


list_files_baseDir=list.files(path = rxnDB,pattern = "_Brite2Mod2KO2Gene.rds",full.names = T,recursive = F)
fileNames=str_replace_all(string = list_files_baseDir,pattern = ".rds",replacement = ".csv")

for(i in 1: length(list_files_baseDir))
{
  rr=readRDS(list_files_baseDir[i])
  write.table(rr,file = fileNames[i],sep = "\t",row.names = F,append = F,quote = F)
  
}

# 306
list_files_baseDir=list.files(path = rxnDB,pattern = "_Overall_KO2Gene2rxn.rds",full.names = T,recursive = F)
fileNames=str_replace_all(string = list_files_baseDir,pattern = ".rds",replacement = ".csv")

for(i in 1: length(list_files_baseDir))
{
  rr=readRDS(list_files_baseDir[i])
  write.table(rr,file = fileNames[i],sep = "\t",row.names = F,append = F,quote = F)
  
}


# 
list_files_baseDir=list.files(path = rxnDB,pattern = "_Mod2KO2Gene.rds",full.names = T,recursive = F)
fileNames=str_replace_all(string = list_files_baseDir,pattern = ".rds",replacement = ".csv")

for(i in 1: length(list_files_baseDir))
{
  rr=readRDS(list_files_baseDir[i])
  write.table(rr,file = fileNames[i],sep = "\t",row.names = F,append = F,quote = F)
  
}


list_files_baseDir=list.files(path = baseDir,pattern = "_shetty_modules_with_genes_integration.rds",full.names = T,recursive = F)
fileNames=str_replace_all(string = list_files_baseDir,pattern = ".rds",replacement = ".csv")

for(i in 1: length(list_files_baseDir))
{
  rr=readRDS(list_files_baseDir[i])
  write.table(rr,file = fileNames[i],sep = "\t",row.names = F,append = F,quote = F)
  
}

list_files_baseDir=list.files(path = rxnDB,pattern = "_Mod_hierarchy2KO2Gene.rds",full.names = T,recursive = F)
fileNames=str_replace_all(string = list_files_baseDir,pattern = ".rds",replacement = ".csv")

for(i in 1: length(list_files_baseDir))
{
  rr=readRDS(list_files_baseDir[i])
  write.table(rr,file = fileNames[i],sep = "\t",row.names = F,append = F,quote = F)
  
}



list_files_baseDir=list.files(path = rxnDB,pattern = "_KO2Pathway2Mod2Gene.rds",full.names = T,recursive = F)
fileNames=str_replace_all(string = list_files_baseDir,pattern = ".rds",replacement = ".csv")

for(i in 1: length(list_files_baseDir))
{
  rr=readRDS(list_files_baseDir[i])
  write.table(rr,file = fileNames[i],sep = "\t",row.names = F,append = F,quote = F)
  
}

list_files_baseDir=list.files(path = baseDir,pattern = "_shetty_module_completeness.rds",full.names = T,recursive = F)
fileNames=str_replace_all(string = list_files_baseDir,pattern = ".rds",replacement = ".csv")

for(i in 1: length(list_files_baseDir))
{
  rr=readRDS(list_files_baseDir[i])
  write.table(rr,file = fileNames[i],sep = "\t",row.names = F,append = F,quote = F)
  
}

list_files_baseDir=list.files(path = baseDir,pattern = "_shetty_module_completeness.rds",full.names = T,recursive = F)
fileNames=str_replace_all(string = list_files_baseDir,pattern = ".rds",replacement = ".csv")

for(i in 1: length(list_files_baseDir))
{
  rr=readRDS(list_files_baseDir[i])
  write.table(rr,file = fileNames[i],sep = "\t",row.names = F,append = F,quote = F)
  
}


list_files=list.files(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/",pattern = "bacteroides_caecemuris_I48",full.names = T,recursive = T)
fl2=gsub("bacteroides_caecemuris_I48", "bacteroides_caecimuris_I48", list_files)
file.rename(list_files,fl2)


require(openxlsx)
FinalExcelFolder="S:/Abilash/Marbouty_annotation_refseq_genbank_2022/newJointAnnotationExcel/"
for(jj in 1: length(genomeNames))
{
  allAnnotationFin=readRDS(str_c(baseDir,genomeNames[jj],"_allAnnotationFin.rds"))
  GO2substrate2Gene=readRDS(str_c(baseDir,genomeNames[jj],"_GO2substrate2Gene.rds"))
  Brite2Mod2KO2Gene=readRDS(str_c(rxnDB,genomeNames[jj],"_Brite2Mod2KO2Gene.rds"))
  Overall_KO2Gene2rxn=readRDS(str_c(rxnDB,genomeNames[jj],"_Overall_KO2Gene2rxn.rds"))
  Mod2KO2Gene=readRDS(str_c(rxnDB,genomeNames[jj],"_Mod_hierarchy2KO2Gene.rds"))
  shetty_integrated=readRDS(str_c(baseDir,genomeNames[jj],"_shetty_modules_with_genes_integration.rds"))
  KO2Pathway2Mod2Gene=readRDS(str_c(path = rxnDB,genomeNames[jj],"_KO2Pathway2Mod2Gene.rds"))
  
  KEGG_module_reconstruct=read.csv(file = str_c("S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Extra/KEGG_and_Shetty_module_reconstruct/KEGG_module_reconstruct_output/",genomeNames[jj],"_reconstructed_modules_df.txt"),header = F,sep = "\t")
  colnames(KEGG_module_reconstruct)=c("KEGG_Module_ID","Number_of_blocks_present_in_genome","Total_number_of_blocks_required_for_module_completeness","Fraction_of_blocks_present","Module_description")
  
  
  Shetty_module_reconstruct=readRDS(file = str_c("S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/",genomeNames[jj],"_shetty_module_completeness.rds"))
  dbCAN2_data=allAnnotationFin[,c("str_Crd","dbCAN2_Cluster","dbCAN2_GeneID","dbCAN2_substrate")]
  dbCAN2_data=dbCAN2_data[dbCAN2_data$dbCAN2_Cluster!="",]
  
  list_of_datasets <- list("Overall_annotations" = allAnnotationFin, "Gene2subtrate_GO" = GO2substrate2Gene, "dbCAN2_Gene2substrate" = dbCAN2_data,
                           "Gene2KEGG_annotations" = Brite2Mod2KO2Gene,
                           "Gene2Reactions_KEGG"= Overall_KO2Gene2rxn,
                           "KEGGModule2Genes"=Mod2KO2Gene,
                           "KEGGPathways2Genes"=KO2Pathway2Mod2Gene,
                           "KEGG_modules_completeness"=KEGG_module_reconstruct,
                           "Metabolic_modules_shetty2Genes" = shetty_integrated,
                           "Shetty_module_completeness"=Shetty_module_reconstruct)
  
  write.xlsx(list_of_datasets, file = str_c(FinalExcelFolder,genomeNames[jj],"_allAnnJoint.xlsx"))
  
  
  
  
  }

