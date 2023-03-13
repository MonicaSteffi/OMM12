
library(readtext)
library("data.table")
library("stringr")

shetty_res=fread("S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Extra/Shetty_et_al_MDbMM16-master/rna_seq/data_files/ModulesCustom/curatedGMMs.v1.07_new.txt",sep = "@",quote = "",header = F)
shetty_hierarchy=fread("S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Extra/Shetty_et_al_MDbMM16-master/data/curated_gmm_trophic_levels.txt",sep = "\t",header = T)
shetty_cat=fread("S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Extra/Shetty_et_al_MDbMM16-master/rna_seq/data_files/ModulesCustom/CuratedGMM_names.txt",sep = "\t",header = F)

end_points=which(shetty_res$V1=="///")
start_points=grep(pattern = "^M",x = shetty_res$V1)

gcIDVec=rep("",length(start_points))
gcInfoVec=rep("",length(start_points))
KOsVec=rep("",length(start_points))

for(i in 1: length(start_points))
{
  
  GC_info=as.character(shetty_res$V1[start_points[i]])
  GC_ID=str_replace_all(string = GC_info,pattern = "\t.*",replacement = "")
  GC_Desc=unlist(str_split(string = GC_info,pattern = "\t"))
  GC_Desc=GC_Desc[-1]
  GC_Desc=str_c(GC_Desc,collapse = " ")
  
  KOs=shetty_res$V1[(start_points[i]+1):(end_points[i]-1)]
  
  KOs=unlist(str_split(string = KOs,pattern = "\t"))
  KOs=unlist(str_split(string = KOs,pattern = ","))
  KOs=str_replace_all(string = KOs,pattern = "[[:punct:]]",replacement = "")
  
  KOs=unlist(str_split(string = KOs,pattern = ","))
  
  KOsVec[i]=str_c(KOs,collapse = ",")
  gcIDVec[i]=GC_ID
  gcInfoVec[i]=GC_Desc
  
  
}

shettyDf=data.frame(module_id=gcIDVec,module_desc=gcInfoVec,KO=KOsVec)
shettyDf=merge(x = shettyDf,y = shetty_hierarchy,by="module_id",all=TRUE)
shettyDf$Degradation=shetty_cat$V3[match(shettyDf$module_id,as.character(shetty_cat$V1))]

saveRDS(file = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Extra/Shetty_et_al_MDbMM16-master/shetty_annotation_integrated.rds",object = shettyDf)
write.table(x = shettyDf,"S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Extra/Shetty_et_al_MDbMM16-master/shetty_annotation_integrated.csv",quote=F,sep = "\t",row.names = F)


shetty_mod2kolist=lapply(shettyDf$KO,function(x){
  
  xx=unlist(str_split(string = x,pattern = ","));
  xx=xx[xx!=""];
  xx=xx[!is.na(xx)];
  return(unique(xx))
})
names(shetty_mod2kolist)=shettyDf$module_id

saveRDS(file = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Extra/Shetty_et_al_MDbMM16-master/shetty_mod2kolist.rds",object = shetty_mod2kolist)
write.table(x = shetty_mod2kolist,"S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Extra/Shetty_et_al_MDbMM16-master/shetty_mod2kolist.csv",quote=F,sep = "\t",row.names = F)



shetty_unqKO=unique(unname(unlist(shetty_mod2kolist)))
shetty_unqMod=names(shetty_mod2kolist)

shetty_mod2koMat=matrix(0,nrow=length(shetty_unqKO),ncol=length(shetty_unqMod))
rownames(shetty_mod2koMat)=shetty_unqKO
colnames(shetty_mod2koMat)=shetty_unqMod

for(i in 1: length(shetty_mod2kolist))
{
  
  tmpKO=shetty_mod2kolist[[i]]
  tmpKO=tmpKO[tmpKO!=""]
  tmpKO=unique(tmpKO[!is.na(tmpKO)])
  
  tmpMod=names(shetty_mod2kolist)[i]
  shetty_mod2koMat[tmpKO,tmpMod]=rep(1,length(tmpKO))
  
}

saveRDS(file = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Extra/Shetty_et_al_MDbMM16-master/shetty_mod2koMat.rds",object = shetty_mod2koMat)
write.table(x = shetty_mod2koMat,"S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Extra/Shetty_et_al_MDbMM16-master/shetty_mod2koMat.csv",quote=F,sep = "\t",row.names = F)


for(i in 1: length(shetty_mod2kolist))
{
  
  tmpKO2=unname(unlist(shetty_mod2kolist[[i]]))
  tmpKO2=tmpKO2[tmpKO2!=""]
  tmpKO2=unique(tmpKO2[!is.na(tmpKO2)])
  
  if(length(tmpKO2)>0)
  {
    
    tmpDf=data.frame(KO=tmpKO2)
    tmpDf$Module=names(shetty_mod2kolist)[[i]]  
    
    if(i==1)
    {
      shetty_mod2koMelt=tmpDf
      
    } else {
      
      shetty_mod2koMelt=rbind(shetty_mod2koMelt,tmpDf)
    }
    
  }
  
  

  }
saveRDS(file = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Extra/Shetty_et_al_MDbMM16-master/shetty_mod2koMelt.rds",object = shetty_mod2koMelt)
write.table(x = shetty_mod2koMelt,"S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Extra/Shetty_et_al_MDbMM16-master/shetty_mod2koMelt.csv",quote=F,sep = "\t",row.names = F)
