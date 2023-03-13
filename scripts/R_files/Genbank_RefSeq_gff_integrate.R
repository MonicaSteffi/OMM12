library("ape")
library("stringr")
MarboutyFolder="S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/"


list_files=list.files(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/",pattern = "Escherichia_coli_Mt1B1_",full.names = T,recursive = T)
fl2=gsub("Escherichia_coli_Mt1B1_", "escherichia_coli_strain_Mt1B1_", list_files)
file.rename(list_files,fl2)

list_files=list.files(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/",pattern = "Escherichia_coli_strain_Mt1B1_",full.names = T,recursive = T)
fl2=gsub("Escherichia_coli_strain_Mt1B1_", "escherichia_coli_strain_Mt1B1_", list_files)
file.rename(list_files,fl2)

list_files=list.files(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/",pattern = "Akkermansia_muciniphila_YL44",full.names = T,recursive = T)
fl2=gsub("Akkermansia_muciniphila_YL44", "akkermansia_muciniphila_YL44", list_files)
file.rename(list_files,fl2)


list_files=list.files(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/",pattern = "Escherichia_coli_strain_Mt1B1_",full.names = T,recursive = T)
fl2=gsub("Escherichia_coli_strain_Mt1B1_", "escherichia_coli_strain_Mt1B1_", list_files)
file.rename(list_files,fl2)


list_dirs=list.dirs(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/", full.names = T,recursive = T)
fl2=gsub("Escherichia_coli_strain_Mt1B1", "escherichia_coli_strain_Mt1B1", list_dirs)
file.rename(list_dirs,fl2)

list_dirs=list.dirs(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/", full.names = T,recursive = T)
fl2=gsub("Escherichia_coli_strain_Mt1B1", "escherichia_coli_strain_Mt1B1", list_dirs)
file.rename(list_dirs,fl2)

list_dirs=list.dirs(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/", full.names = T,recursive = T)
fl2=gsub("Akkermansia", "akkermansia", list_dirs)
file.rename(list_dirs,fl2)

list_dirs=list.dirs(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/", full.names = T,recursive = T)
fl2=gsub("muribaculum_intestinale_YL27", "muribaculum_intestinales_YL27", list_dirs)
file.rename(list_dirs,fl2)


list_files_genbank=list.files(path="S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_genbank_and_refseq_gffFiles/Marbouty_Genbank/",full.names = T)
gnkNm=str_replace_all(string = list_files_genbank,pattern = ".*Marbouty_Genbank/",replacement = "")
gnkNm=str_replace_all(string = gnkNm,pattern = "_genomic.gff",replacement = "")

setClass(Class = "gffAnnotate",representation("FinalGff"="data.frame","FinalGff2"="data.frame",
                                              "FinalGff3"="data.frame","FinalGff4"="data.frame",
                                              overlapGn="list",sbstrr="list",
                                              RefSeq_coord="data.frame",Genbank_coord="data.frame"))

gnmID2Org=read.csv(file = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/GenomeID_to_Organism.csv",header = T,sep = ";")
gnmID2Org$GenbankID_edited=str_replace_all(string = gnmID2Org$GenbankID,pattern = "\\..*",replacement = "")

# list_files_refseq=list.files(path="S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_genbank_and_refseq_gffFiles/Marbouty_RefSeq/",full.names = T)
list_files_refseq=list.files(path="S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_genbank_and_refseq_gffFiles/Marbouty_RefSeq/",full.names = T,recursive = T)
RfqNm=str_replace_all(string = list_files_refseq,pattern = ".*Marbouty_RefSeq/",replacement = "")
RfqNm=str_replace_all(string = RfqNm,pattern = "_genomic.gff",replacement = "")

ecol_refseq="S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_genbank_and_refseq_gffFiles/Marbouty_RefSeq/Escherichia_coli_Mt1B1_genomic.gff"
ecol_rfqnm="Escherichia_coli_Mt1B1"

ecol_gnk="S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_genbank_and_refseq_gffFiles/Marbouty_Genbank/Escherichia_coli_Mt1B1_genomic.gff"
ecol_gnknm="Escherichia_coli_Mt1B1"


ovLnk=c(list_files_genbank,list_files_refseq)
ovLnkDf=data.frame(ovLnk=ovLnk)
ovLnkDf$ovGff=ovLnkDf$ovLnk
ovLnkDf$ovGff=str_replace_all(string = ovLnkDf$ovGff,pattern = ".*Marbouty_RefSeq/",replacement = "")
ovLnkDf$ovGff=str_replace_all(string = ovLnkDf$ovGff,pattern = ".*Marbouty_Genbank/",replacement = "")
ovLnkDf$ovGff=str_replace_all(string = ovLnkDf$ovGff,pattern = "_genomic.gff",replacement = "")

ovLnkDf$Source=str_replace_all(string = ovLnkDf$ovLnk,pattern = ".*Marbouty_",replacement = "")
ovLnkDf$Source=str_replace_all(string = ovLnkDf$Source,pattern = "/.*",replacement = "")

ovGff=ovLnkDf$ovGff
ovLnk=ovLnkDf$ovLnk
unqGff=unique(ovGff)
srcTyp=as.character(ovLnkDf$Source)


lst2df=function(attrList)
{
  
unqNamesAttr=unique(names(attrList))
attrVec=vector()
for(unqNmCnt in 1: length(unqNamesAttr))
{
  tmpList=unique(unlist(attrList[grep(pattern = unqNamesAttr[unqNmCnt],x=names(attrList))]))
  tmpList=str_split(string = tmpList,pattern = "=")
  for(t1 in 1: length(tmpList))
  {
    tmpList[[t1]][1]=str_c(tmpList[[t1]][1],"_",unqNamesAttr[unqNmCnt],"=")
    
  }
  attrVec[unqNmCnt]=str_c(unlist(lapply(tmpList, function(x){
    str_c(x,collapse = "")
  })),collapse = ";")
  
  }

return(str_c(attrVec,collapse = ";"))
}
# lst2df=function(attrList)
# {
#   attrRws=unique(unlist(lapply(attrList,function(x){str_replace_all(string = x,pattern = "=.*",replacement = "")})))
#   attrDf=data.frame(Attributes=attrRws)
#   
#  for(i in 1:length(attrList))
#  {
#    Ai=unlist(attrList[i])
#    tmpVec=unlist(lapply(as.character(attrDf$Attributes),function(x){grep(pattern = str_c("^",x,"="),x=unlist(Ai))}))
#    
#    matchingRws=match(str_replace_all(unlist(Ai),pattern = "=.*",replacement = ""),attrDf$Attributes)
#    
#    attrDf [matchingRws,(1+i)]=Ai[tmpVec]
#    attrDf[,(1+i)]=str_replace_all(string = attrDf[,(1+i)],pattern = ".*=",replacement = "")
#  }
#   
#   colnames(attrDf)[2:ncol(attrDf)]=names(attrList)
#   attrDf[is.na(attrDf)]="Not_available"
#   
#   for(j3 in 2: ncol(attrDf))
#   {
#     
#     attrDf[,j3]=str_c(colnames(attrDf)[j3],"=",attrDf[,j3])
#     
#   }
#   
#   for(j4 in 1:nrow(attrDf))
#   {
#     
#     attrDf[j4,-1]=str_c(attrDf[j4,1],"_",attrDf[j4,-1])
#     
#   }
#   
#   # attrVec=vector()
#   # for(j in 1: nrow(attrDf))
#   # {
#   #   unqAttr=unique(as.character(attrDf[j,2:ncol(attrDf)]))
#   #   
#   #   if(length(unqAttr)==1)
#   #   {
#   #     attrVec[j]=str_c(attrDf[j,1],"=",unqAttr)
#   #     
#   #   } else {
#   #     
#   #     attrVec[j]=str_c(str_c(attrDf[j,1],"_",(str_c(colnames(attrDf)[-1],"=",unqAttr))),collapse = ";")
#   #   }
#   # 
#   #   
#   #   }
#   # 
#   return(str_c(unique(unlist(attrDf[,-1])),collapse = ";"))
#   }

## Old function:
## Had RS deleted for old_locus_tag;
## Not working with Blautia, where gene with Locus_tag I5Q86_RS20695 had an old locus tag I5Q86_20700.
## That is solved with the new function unique_gffv2.
unique_gff=function(rr1,rr2,flType)
{
  colnames(rr1)=str_c("RefSeq_",colnames(rr1))
  colnames(rr2)=str_c("Genbank_",colnames(rr2))
  
  stdff1=setdiff(unique(rr1$RefSeq_ID),unique(rr2$Genbank_ID))
  intsct=intersect(unique(rr1$RefSeq_ID),unique(rr2$Genbank_ID))
  stdff2=setdiff(unique(rr2$Genbank_ID),unique(rr1$RefSeq_ID))
 
  
  #### Adding Locus_tag column
   
  rr1$Locus_tag=rr1$RefSeq_attributes
  rr1$Locus_tag=str_replace_all(string = rr1$Locus_tag,pattern = ".*old_locus_tag=",replacement = "")
  rr1$Locus_tag=str_replace_all(string = rr1$Locus_tag,pattern = ";.*",replacement = "")
  rr1$Locus_tag=str_replace_all(string = rr1$Locus_tag,pattern = "ID=.*",replacement = "")
  rr1$Locus_tag=str_replace_all(string = rr1$Locus_tag,pattern = "%2.*",replacement = "")
  
  unq_ID=unique(rr1$RefSeq_ID)
  emptyLoc=which(rr1$Locus_tag=="")
  
  for(i in 1 : length(unq_ID))
  {
  tmpID=which(rr1$RefSeq_ID==unq_ID[i])
  tmpEmp=intersect(tmpID,emptyLoc)
    
  if(length(tmpID) ==2 && length(tmpEmp) ==2)
  {
    tmpLc=str_replace_all(string = rr1$RefSeq_attributes[tmpID],pattern = ".*;locus_tag=",replacement = "")
    tmpLc=unique(str_replace_all(string = tmpLc,pattern = ";.*",replacement = ""))
    tmpLc=tmpLc[tmpLc!=""]
    tmpLc=str_replace_all(string = tmpLc,pattern = "RS",replacement = "")
    rr1$Locus_tag[tmpID]=tmpLc
    
  } else {
    if(length(tmpID)>1 && length(tmpEmp) == 1)
    {
      rr1$Locus_tag[tmpEmp]=rr1$Locus_tag[setdiff(tmpID,tmpEmp)]
      
    } 
    
    if(length(tmpID) ==1 && length(tmpEmp) ==1)
    {
      
      tmpLc=str_replace_all(string = rr1$RefSeq_attributes[tmpID],pattern = ".*;locus_tag=",replacement = "")
      tmpLc=unique(str_replace_all(string = tmpLc,pattern = ";.*",replacement = ""))
      tmpLc=tmpLc[tmpLc!=""]
      tmpLc=str_replace_all(string = tmpLc,pattern = "RS",replacement = "")
      rr1$Locus_tag[tmpID]=tmpLc
      
    }
    
  }
  
  
  
  
  }
  
  # biotype=rRNA
  rr2$Locus_tag=rr2$Genbank_attributes
  rr2$Locus_tag=str_replace_all(string = rr2$Locus_tag,pattern = ".*locus_tag=",replacement = "")
  rr2$Locus_tag=str_replace_all(string = rr2$Locus_tag,pattern = ";.*",replacement = "")
  rr2$Locus_tag=str_replace_all(string = rr2$Locus_tag,pattern = "%2.*",replacement = "")
  ###################################################################################################
  ##### Category-1: Joining all the rows with same ID:
  ###################################################################################################
  unq_RefSeq=rr1[(rr1$RefSeq_ID %in% intsct),]
  unq_GenBank=rr2[(rr2$Genbank_ID %in% intsct),]
  clnms=str_replace_all(string = colnames(rr1),pattern = "RefSeq_",replacement = "")
  colnames(unq_RefSeq)=clnms
  colnames(unq_GenBank)=clnms
  
  unq_RefSeq$typ_st_end_str=apply(unq_RefSeq[,c(3:5,7)],1,function(x){str_c(str_trim(x,side = "both"),collapse="_")})
  unq_GenBank$typ_st_end_str=apply(unq_GenBank[,c(3:5,7)],1,function(x){str_c(str_trim(x,side="both"),collapse="_")})
  unq_RefSeq$DB="RefSeq"
  unq_GenBank$DB="GenBank"
  
  OvGff=rbind(unq_RefSeq,unq_GenBank)
  unique_typID=unique(OvGff$typ_st_end_str)
  
  for(i in 1:length(unique_typID))
  {
    tmptypID=which(OvGff$typ_st_end_str==unique_typID[i])
    
    IntsctGff=data.frame(seqid=str_c(unique(OvGff$seqid[tmptypID]),collapse ="&"),source=str_c(unique(OvGff$source[tmptypID]),collapse="&"),type=str_c(unique(OvGff$type[tmptypID]),collapse="&"))
    IntsctGff$start=OvGff$start[tmptypID[1]]
    IntsctGff$end=OvGff$end[tmptypID[1]]
    IntsctGff$score=str_c(unique(OvGff$score[tmptypID]),collapse="&")
    IntsctGff$strand=OvGff$strand[tmptypID[1]]
    IntsctGff$phase=str_c(unique(OvGff$phase[tmptypID]),collapse="&")
    
    attrList=(lapply(OvGff$attributes[tmptypID],function(x){unlist(str_split(string = x,pattern = ";"))}))
    names(attrList)=OvGff$DB[tmptypID]
   
    attrVec=lst2df(attrList)
    attrVec2=str_c("Start_End_strand=",str_trim(unique_typID[i],side= "both"))
    attrVec3=str_replace_all(string = OvGff$Locus_tag[tmptypID],pattern = "^ID=",replacement = "")
    attrVec3=str_c("Locus_tag=",str_c(str_trim(unique(attrVec3),side = "both"),collapse="&"))
    ovAttr=str_c(attrVec,";",attrVec2,";",attrVec3)
    IntsctGff$attributes=ovAttr
    if(i==1)
    {
     FinalGff=IntsctGff
      
    } else {
      
      FinalGff=rbind(FinalGff,IntsctGff)
    }
    
    
  }
  
  
  
  ###################################################################################################
  ##### Category-2: Those with different ID, Joining all the rows with same Locus_tag
  ###################################################################################################
  
  stdfrr1=rr1[(rr1$RefSeq_ID %in% stdff1),]
  stdfrr2=rr2[(rr2$Genbank_ID %in% stdff2),]
  
  merge_std=(merge(stdfrr1,stdfrr2,by="Locus_tag"))
  merge_std=merge_std[!duplicated(merge_std),]
  
  Cat2unq_RefSeq=merge_std[,grep(pattern = "RefSeq",x=colnames(merge_std))]
  Cat2unq_RefSeq$Locus_tag=merge_std$Locus_tag
  Cat2unq_RefSeq$typ_st_end_str=apply(Cat2unq_RefSeq[,c("RefSeq_type","Locus_tag")],1,function(x){str_c(str_trim(x,side = "both"),collapse="_")})
  
  Cat2unq_GenBank=merge_std[,grep(pattern = "Genbank",x=colnames(merge_std))]
  Cat2unq_GenBank$Locus_tag=merge_std$Locus_tag
  Cat2unq_GenBank$typ_st_end_str=apply(Cat2unq_GenBank[,c("Genbank_type","Locus_tag")],1,function(x){str_c(str_trim(x,side="both"),collapse="_")})
  
  Cat2unq_RefSeq$DB="RefSeq"
  Cat2unq_GenBank$DB="GenBank"
  colnames(Cat2unq_GenBank)=str_replace_all(string = colnames(Cat2unq_GenBank),pattern = "Genbank_",replacement = "")
  colnames(Cat2unq_RefSeq)=str_replace_all(string = colnames(Cat2unq_RefSeq),pattern = "RefSeq_",replacement = "")
  
  Cat2OvGff=rbind(Cat2unq_RefSeq,Cat2unq_GenBank)
  unique_typID=unique(Cat2OvGff$typ_st_end_str)
  
  if(length(unique_typID)>0)
    {
  for(i in 1:length(unique_typID))
  {
    tmptypID=which(Cat2OvGff$typ_st_end_str==unique_typID[i])
    
    Cat2IntsctGff=data.frame(seqid=str_c(unique(Cat2OvGff$seqid[tmptypID]),collapse ="&"),source=str_c(unique(Cat2OvGff$source[tmptypID]),collapse="&"),type=str_c(unique(Cat2OvGff$type[tmptypID]),collapse="&"))
    Lngs=sort(apply(Cat2OvGff[tmptypID,c(4,5)],1,function(x){abs(x[1]-x[2])+1}),decreasing = T)
    
    Cat2IntsctGff$start=Cat2OvGff$start[as.numeric(names(Lngs)[1])]
    Cat2IntsctGff$end=Cat2OvGff$end[as.numeric(names(Lngs)[1])]
    Cat2IntsctGff$score=str_c(unique(Cat2OvGff$score[tmptypID]),collapse="&")
    Cat2IntsctGff$strand=Cat2OvGff$strand[tmptypID[1]]
    Cat2IntsctGff$phase=str_c(unique(Cat2OvGff$phase[tmptypID]),collapse="&")
    
    attrList=(lapply(Cat2OvGff$attributes[tmptypID],function(x){unlist(str_split(string = x,pattern = ";"))}))
    names(attrList)=Cat2OvGff$DB[tmptypID]
    
    attrVec=lst2df(attrList)
    attrVec2=str_c("Start_End_strand=",str_trim(unique_typID[i],side= "both"))
    attrVec3=str_replace_all(string = Cat2OvGff$Locus_tag[tmptypID],pattern = "^ID=",replacement = "")
    attrVec3=str_c("Locus_tag=",str_c(str_trim(unique(attrVec3),side = "both"),collapse="&"))
    
    attrTmp=unique(apply(Cat2OvGff[tmptypID,c(2,4,5)],1,function(x){str_c(x[1],"_start=",x[2],";",x[1],"_end=",x[3])}))
    attrVec4=paste(attrTmp,collapse=";")
    
    ovAttr=str_c(attrVec,";",attrVec2,";",attrVec3,";",attrVec4)
    Cat2IntsctGff$attributes=ovAttr
    if(i==1)
    {
      FinalGff2=Cat2IntsctGff
      
    } else {
      
      FinalGff2=rbind(FinalGff2,Cat2IntsctGff)
    }
    
    
  }
  } else {
    
    FinalGff2=data.frame()
  }
  
  ###############################################################################################################
  ##### Category-3: Those with different ID and Locus_tags, Joining all the rows subset of another gene : sbstrr
  ##### Category-4: Those with different ID and Locus_tags, Joining all the rows where genes overlap : overlapGn
  ###############################################################################################################
  
  unq_std=setdiff(unique(c(stdff1,stdff2)),unique(c(merge_std$RefSeq_ID, merge_std$Genbank_ID)))
  unq_std_locus=c(stdfrr1$Locus_tag[match(unq_std,stdfrr1$RefSeq_ID)],
                  stdfrr2$Locus_tag[match(unq_std,stdfrr2$Genbank_ID)])
  names(unq_std_locus)=c(stdfrr1$RefSeq_ID[match(unq_std,stdfrr1$RefSeq_ID)],
                         stdfrr2$Genbank_ID[match(unq_std,stdfrr2$Genbank_ID)])
  
  unq_std_locus=unq_std_locus[!is.na(unq_std_locus)]
  
  
  rr1_unq_std=rr1[(rr1$RefSeq_ID %in% unq_std),]
  rr2_unq_std=rr2[(rr2$Genbank_ID %in% unq_std),]
  
  rr1_unq_std$GnCoord=apply(rr1_unq_std[,c("RefSeq_start", "RefSeq_end")],1,function(x){str_c(as.character(seq(x[1],x[2],by=1)),collapse = ",")})
  rr2_unq_std$GnCoord=apply(rr2_unq_std[,c("Genbank_start", "Genbank_end")],1,function(x){str_c(as.character(seq(x[1],x[2],by=1)),collapse = ",")})
  
  
  sbstrr1= vector("list", length = nrow(rr1_unq_std))
  sbstrr2=vector("list", length = nrow(rr2_unq_std))
  names(sbstrr1)=rr1_unq_std$Locus_tag
  names(sbstrr2)=rr2_unq_std$Locus_tag
  
  overlapGn=list()
  k=1
 if(nrow(rr1_unq_std)==0 && nrow(rr2_unq_std)!=0)
 {
   colnames(rr2_unq_std)=str_replace_all(string = colnames(rr2_unq_std),pattern = "RefSeq_|Genbank_",replacement = "")
   rr2_unq_std$attributes=str_c(rr2_unq_std$attributes,";Locus_tag=",rr2_unq_std$Locus_tag)
   FinalGff2=rbind(FinalGff2,rr2_unq_std[,1:ncol(FinalGff2)])
   FinalGff3=data.frame()
   FinalGff4=data.frame()
 } else {
   
   if(nrow(rr2_unq_std)==0 && nrow(rr1_unq_std) !=0)
   {
     colnames(rr1_unq_std)=str_replace_all(string = colnames(rr1_unq_std),pattern = "RefSeq_|Genbank_",replacement = "")
     rr1_unq_std$attributes=str_c(rr1_unq_std$attributes,";Locus_tag=",rr1_unq_std$Locus_tag)
     FinalGff2=rbind(FinalGff2,rr1_unq_std[,1:ncol(FinalGff2)])
     FinalGff3=data.frame()  
     FinalGff4=data.frame()
   } else {
     if(nrow(rr1_unq_std) !=0 && nrow(rr2_unq_std) != 0)
     {
     
     for(i in 1:nrow(rr1_unq_std))
     {
       gnCrds=as.numeric(unlist(str_split(rr1_unq_std$GnCoord[i],pattern = ",")  ))
       
       
       for(j in 1:nrow(rr2_unq_std))
       {
         
         gnCrds2=as.numeric(unlist(str_split(rr2_unq_std$GnCoord[j],pattern = ",")  ))
         instCrds=intersect(gnCrds,gnCrds2)  
         
         if(length(instCrds)==length(gnCrds) && rr1_unq_std$RefSeq_strand[i] == rr2_unq_std$Genbank_strand[j])
         {
           
           sbstrr1[rr1_unq_std$Locus_tag[i]]=c(unname(unlist(sbstrr1[rr1_unq_std$Locus_tag[i]])),rr2_unq_std$Locus_tag[j])
         }
         
         if(length(instCrds)==length(gnCrds2) && rr1_unq_std$RefSeq_strand[i] == rr2_unq_std$Genbank_strand[j] )
         {
           
           sbstrr2[rr2_unq_std$Locus_tag[j]]=c(unname(unlist(sbstrr2[rr2_unq_std$Locus_tag[j]])),rr1_unq_std$Locus_tag[i])
         }
         
         if(length(instCrds)>1 && rr1_unq_std$RefSeq_strand[i] == rr2_unq_std$Genbank_strand[j])
         {
           if(length(instCrds) != length(gnCrds) && length(instCrds)!=length(gnCrds2) )
           {
             overlapGn[[k]]=c(rr1_unq_std$Locus_tag[i],rr2_unq_std$Locus_tag[j])
             k=k+1  
           }
           
         }
       }
     }
     
     
     ## Remove NULL lists, then remove redundant lists.
     overlapGn=Filter(Negate(is.null), overlapGn)
     
     res.list <- overlapGn
     # take set difference between contents of list elements and accumulated elements
     res.list[-1] <- mapply("setdiff", res.list[-1],
                            head(Reduce(c, overlapGn, accumulate=TRUE), -1))
     
     overlapGn=res.list[lapply(res.list,length)>0]
     
     
     sbstrr1=Filter(Negate(is.null),sbstrr1)
     res.list <- sbstrr1
     
     # take set difference between contents of list elements and accumulated elements
     res.list[-1] <- mapply("setdiff", res.list[-1],
                            head(Reduce(c, sbstrr1, accumulate=TRUE), -1))
     
     sbstrr1=res.list[lapply(res.list,length)>0]
     
     sbstrr2=Filter(Negate(is.null),sbstrr2)
     res.list <- sbstrr2
     # take set difference between contents of list elements and accumulated elements
     res.list[-1] <- mapply("setdiff", res.list[-1],
                            head(Reduce(c, sbstrr2, accumulate=TRUE), -1))
     sbstrr2=res.list[lapply(res.list,length)>0]
     sbstrr=c(sbstrr1,sbstrr2)
     
     ### 1. Create gff file for these overlapping and subset genes:
     ### 2. Add Rows where the overlapping region and non-overlapping region have separate rows.
     
     if(length(sbstrr)>0)
     {
       for(i in 1: length(sbstrr))
       {
         
         tmpG1rr1=rr1_unq_std[( rr1_unq_std$Locus_tag %in% unname(sbstrr)[i]),]
         tmpG1rr2=rr2_unq_std[( rr2_unq_std$Locus_tag %in% unname(sbstrr)[i]),]
         
         tmpG2rr1=rr1_unq_std[( rr1_unq_std$Locus_tag %in% names(sbstrr)[i]),]
         tmpG2rr2=rr2_unq_std[( rr2_unq_std$Locus_tag %in% names(sbstrr)[i]),]
         
         colnames(tmpG1rr1)=str_replace_all(string = colnames(tmpG1rr1),pattern = "RefSeq_|Genbank_",replacement = "")
         colnames(tmpG1rr2)=str_replace_all(string = colnames(tmpG1rr2),pattern = "RefSeq_|Genbank_",replacement = "")
         colnames(tmpG2rr1)=str_replace_all(string = colnames(tmpG2rr1),pattern = "RefSeq_|Genbank_",replacement = "")
         colnames(tmpG2rr2)=str_replace_all(string = colnames(tmpG2rr2),pattern = "RefSeq_|Genbank_",replacement = "")
         
         tmpG1rr1$Length=apply(tmpG1rr1[,c(4,5)],1,function(x){abs(x[1]-x[2])+1})
         tmpG1rr2$Length=apply(tmpG1rr2[,c(4,5)],1,function(x){abs(x[1]-x[2])+1})
         
         tmpG2rr1$Length=apply(tmpG2rr1[,c(4,5)],1,function(x){abs(x[1]-x[2])+1})
         tmpG2rr2$Length=apply(tmpG2rr2[,c(4,5)],1,function(x){abs(x[1]-x[2])+1})
         
         gh=rbind(tmpG1rr1,tmpG1rr2,tmpG2rr1,tmpG2rr2)
         gh=gh[order(gh$Length),]
         gh$Set=""
        if(length(unique(gh$Locus_tag))==1)
        {
          ghGff3=gh 
        } else { 
         bgLoc=unique(gh$Locus_tag[which(gh$Length==max(gh$Length))])
         smLoc=unique(gh$Locus_tag[which(gh$Length!=max(gh$Length))])
         
         for(jj in 1:length(smLoc))
         {
           smGnCrd=gh$GnCoord[which(gh$Locus_tag==smLoc[jj])]
           bgGnCrd=gh$GnCoord[which(gh$Locus_tag==bgLoc)]
           
           Newbg1=gh[gh$Locus_tag==bgLoc,]
           Newbg1$Set="Overlapping"
           
           Newbg2=gh[gh$Locus_tag==bgLoc,]
           Newbg2$Set="Non_overlapping"
           
           Newbg=rbind(Newbg1,Newbg2)
           
           intsctGnCrd=intersect(as.numeric(unlist(str_split(string = bgGnCrd,pattern = ","))),
                                 as.numeric(unlist(str_split(string = smGnCrd,pattern = ","))))
           
           stdffsctGnCrd=setdiff(as.numeric(unlist(str_split(string = bgGnCrd,pattern = ","))),
                                 as.numeric(unlist(str_split(string = smGnCrd,pattern = ","))))
           
           if(all(Newbg$strand=="+"))
           {
             
             Newbg1$start=min(intsctGnCrd)
             Newbg1$end=max(intsctGnCrd)
             
             Newbg2$start=min(stdffsctGnCrd)
             Newbg2$end=max(stdffsctGnCrd)
             
             Newbg=rbind(Newbg1,Newbg2)
             
             
           } 
           
           if(all(Newbg$strand=="-"))
           {
             
             
             Newbg1$start=max(intsctGnCrd)
             Newbg1$end=min(intsctGnCrd)
             
             Newbg2$start=max(stdffsctGnCrd)
             Newbg2$end=min(stdffsctGnCrd)
             
             Newbg=rbind(Newbg1,Newbg2)
             
             
             
             
             
           }
           
           if(jj==1)
           {
             ExtraRws=Newbg
             
           } else {
             
             ExtraRws=rbind(ExtraRws,Newbg)
           }
           
         }
         
         ghGff3=rbind(gh,ExtraRws)
       }
         if(i==1)
         {
           FinGff3=ghGff3
         } else {
           
           FinGff3=rbind(FinGff3,ghGff3)
         }
       }
       
       FinalGff3=FinGff3[,1:9]
       
       for(i in 1:nrow(FinGff3))
       {
         tmpLocus=ifelse(FinGff3$Set[i]!="",str_c(FinGff3$Locus_tag[i],"_",FinGff3$Set[i]  ),FinGff3$Locus_tag[i])
         FinalGff3$attributes[i]=str_c(FinalGff3$attributes[i],";Locus_tag=",tmpLocus)
         
       }
       
       
     } else {
       
       FinalGff3=data.frame()
     }
     
   if(length(overlapGn)>0)
     {
       for(i in 1: length(overlapGn))
       {
         tmpGenes=unlist(overlapGn[[i]])
         
         tmpG1rr1=rr1_unq_std[( rr1_unq_std$Locus_tag %in% tmpGenes[1]),]
         tmpG1rr2=rr2_unq_std[( rr2_unq_std$Locus_tag %in% tmpGenes[2]),]
         
         tmpG2rr1=rr1_unq_std[( rr1_unq_std$Locus_tag %in% tmpGenes[1]),]
         tmpG2rr2=rr2_unq_std[( rr2_unq_std$Locus_tag %in% tmpGenes[2]),]
         
         colnames(tmpG1rr1)=str_replace_all(string = colnames(tmpG1rr1),pattern = "RefSeq_|Genbank_",replacement = "")
         colnames(tmpG1rr2)=str_replace_all(string = colnames(tmpG1rr2),pattern = "RefSeq_|Genbank_",replacement = "")
         colnames(tmpG2rr1)=str_replace_all(string = colnames(tmpG2rr1),pattern = "RefSeq_|Genbank_",replacement = "")
         colnames(tmpG2rr2)=str_replace_all(string = colnames(tmpG2rr2),pattern = "RefSeq_|Genbank_",replacement = "")
         
         
         gh=rbind(tmpG1rr1,tmpG1rr2,tmpG2rr1,tmpG2rr2)
         gh=gh[!duplicated(gh),]
         gh$Set=""
        
         if(length(unique(gh$Locus_tag))==1)
         {
           
           ghGff4=gh
         } else { 
         
         smLoc=unique(gh$Locus_tag)[1]
         bgLoc=unique(gh$Locus_tag)[2]
         
         smGnCrd=unique(gh$GnCoord[which(gh$Locus_tag==smLoc)])
         bgGnCrd=unique(gh$GnCoord[which(gh$Locus_tag==bgLoc)])
         
         Newbg1=gh[gh$Locus_tag==bgLoc,]
         Newbg1$Set="Overlapping"
         
         Newbg2=gh[gh$Locus_tag==bgLoc,]
         Newbg2$Set="Non_overlapping"
         
         
         Newbg3=gh[gh$Locus_tag==smLoc,]
         Newbg3$Set="Overlapping"
         
         Newbg4=gh[gh$Locus_tag==smLoc,]
         Newbg4$Set="Non_overlapping"
         
         Newbg=rbind(Newbg1,Newbg2,Newbg3,Newbg4)
         
         intsctGnCrd=intersect(as.numeric(unlist(str_split(string = bgGnCrd,pattern = ","))),
                               as.numeric(unlist(str_split(string = smGnCrd,pattern = ","))))
         
         stdffsctGnCrdbg=setdiff(as.numeric(unlist(str_split(string = bgGnCrd,pattern = ","))),
                                 as.numeric(unlist(str_split(string = smGnCrd,pattern = ","))))
         
         stdffsctGnCrdsm=setdiff(as.numeric(unlist(str_split(string = smGnCrd,pattern = ","))),
                                 as.numeric(unlist(str_split(string = bgGnCrd,pattern = ","))))
         
         if(all(Newbg$strand=="+"))
         {
           
           Newbg1$start=min(intsctGnCrd)
           Newbg1$end=max(intsctGnCrd)
           
           Newbg2$start=min(stdffsctGnCrdbg)
           Newbg2$end=max(stdffsctGnCrdbg)
           
           
           Newbg3$start=min(intsctGnCrd)
           Newbg3$end=max(intsctGnCrd)
           
           Newbg4$start=min(stdffsctGnCrdsm)
           Newbg4$end=max(stdffsctGnCrdsm)
           
           
           Newbg=rbind(Newbg1,Newbg2,Newbg3,Newbg4)
           
           
         } 
         
         if(all(Newbg$strand=="-"))
         {
           
           Newbg1$start=min(intsctGnCrd)
           Newbg1$end=max(intsctGnCrd)
           
           Newbg2$start=min(stdffsctGnCrdbg)
           Newbg2$end=max(stdffsctGnCrdbg)
           
           Newbg3$start=min(intsctGnCrd)
           Newbg3$end=max(intsctGnCrd)
           
           Newbg4$start=min(stdffsctGnCrdsm)
           Newbg4$end=max(stdffsctGnCrdsm)
           
           Newbg=rbind(Newbg1,Newbg2,Newbg3,Newbg4)
           
         }
         
         ghGff4=rbind(gh,Newbg)
       }  
         if(i==1)
         {
           FinGff4=ghGff4
         } else {
           
           FinGff4=rbind(FinGff4,ghGff4)
         }
       }
       
       FinalGff4=FinGff4[,1:9]
       
       for(i in 1:nrow(FinGff4))
       {
         tmpLocus=ifelse(FinGff4$Set[i]!="",str_c(FinGff4$Locus_tag[i],"_",FinGff4$Set[i]  ),FinGff4$Locus_tag[i])
         FinalGff4$attributes[i]=str_c(FinalGff4$attributes[i],";Locus_tag=",tmpLocus)
         
       }
       
     } else {
       FinalGff4=data.frame()
       
     }
     
     
     
     } else {
       
       FinalGff3=data.frame()
       FinalGff4=data.frame()
     }
     }
 }
  
  #### For intersect, all the attributes in V9 for each unique protein
  # View(rbind(rr1[(rr1$RefSeq_ID %in% intsct),],rr2[(rr2$Genbank_ID %in% intsct),]))
  
  #### For setdiff, check whether they are subset of another protein:
  
  ###### If yes, check whether they both are with same locus_tag in Genbank and old_locus_tag in RefSeq.
  
  
  ######### If yes, then include only the start and end part of the largest region.
  ######### If no, then they both are different. Retain both of them as separate proteins.
  ###### If no, then they both are different. Retain both of them as separate proteins.
  
return(new(Class = "gffAnnotate",FinalGff=FinalGff,FinalGff2=FinalGff2,
           FinalGff3=FinalGff3,FinalGff4=FinalGff4,overlapGn=overlapGn,sbstrr=sbstrr))  
  
  
  }

unique_gffv2=function(rr1,rr2,flType)
{
  colnames(rr1)=str_c("RefSeq_",colnames(rr1))
  colnames(rr2)=str_c("Genbank_",colnames(rr2))
  
  merge_std_cols=c("RefSeq_seqid","RefSeq_source", "RefSeq_type","RefSeq_start","RefSeq_end","RefSeq_score","RefSeq_strand",
  "RefSeq_phase","RefSeq_attributes","RefSeq_ID","Locus_tag","RefSeq_GeneID","Refseq_length")    
  
  stdff1=setdiff(unique(rr1$RefSeq_ID),unique(rr2$Genbank_ID))
  intsct=intersect(unique(rr1$RefSeq_ID),unique(rr2$Genbank_ID))
  stdff2=setdiff(unique(rr2$Genbank_ID),unique(rr1$RefSeq_ID))
  
  ##################################################################################################################################################
  ### Each gene has most of the times multiple columns of type: CDS, gene etc. But the locus tag could be there in only one of the columns. To fill in the locus tags in such 
  #### such empty columns, we do the following step in Refseq object.
  #### Adding Locus_tag column
  ##################################################################################################################################################  
  rr1$Locus_tag=rr1$RefSeq_attributes
  rr1$Locus_tag=str_replace_all(string = rr1$Locus_tag,pattern = ".*;old_locus_tag=",replacement = "")
  rr1$Locus_tag=str_replace_all(string = rr1$Locus_tag,pattern = ";.*",replacement = "")
  rr1$Locus_tag=str_replace_all(string = rr1$Locus_tag,pattern = "ID=.*",replacement = "")
  rr1$Locus_tag=str_replace_all(string = rr1$Locus_tag,pattern = "%2.*",replacement = "")
  
  unq_ID=unique(rr1$RefSeq_ID)
  emptyLoc=which(rr1$Locus_tag=="")
  
  for(i in 1 : length(unq_ID))
  {
    tmpID=which(rr1$RefSeq_ID==unq_ID[i])
    tmpEmp=intersect(tmpID,emptyLoc)
    
    if(length(tmpID) ==2 && length(tmpEmp) ==2)
    {
      tmpLc=str_replace_all(string = rr1$RefSeq_attributes[tmpID],pattern = ".*;locus_tag=",replacement = "")
      tmpLc=unique(str_replace_all(string = tmpLc,pattern = ";.*",replacement = ""))
      tmpLc=tmpLc[tmpLc!=""]
      # tmpLc=str_replace_all(string = tmpLc,pattern = "RS",replacement = "")
      rr1$Locus_tag[tmpID]=tmpLc
      
    } else {
      if(length(tmpID)>1 && length(tmpEmp) == 1)
      {
        rr1$Locus_tag[tmpEmp]=rr1$Locus_tag[setdiff(tmpID,tmpEmp)]
        
      } 
      
      if(length(tmpID) ==1 && length(tmpEmp) ==1)
      {
        
        tmpLc=str_replace_all(string = rr1$RefSeq_attributes[tmpID],pattern = ".*;locus_tag=",replacement = "")
        tmpLc=unique(str_replace_all(string = tmpLc,pattern = ";.*",replacement = ""))
        tmpLc=tmpLc[tmpLc!=""]
        # tmpLc=str_replace_all(string = tmpLc,pattern = "RS",replacement = "")
        rr1$Locus_tag[tmpID]=tmpLc
        
      }
      
    }
    
    
    
    
  }
  
  # biotype=rRNA
  rr2$Locus_tag=rr2$Genbank_attributes
  rr2$Locus_tag=str_replace_all(string = rr2$Locus_tag,pattern = ".*;locus_tag=",replacement = "")
  rr2$Locus_tag=str_replace_all(string = rr2$Locus_tag,pattern = ";.*",replacement = "")
  rr2$Locus_tag=str_replace_all(string = rr2$Locus_tag,pattern = "%2.*",replacement = "")
  ###################################################################################################
  ##### Category-1: Joining all the rows with same ID:
  ###################################################################################################
  unq_RefSeq=rr1[(rr1$RefSeq_ID %in% intsct),]
  unq_GenBank=rr2[(rr2$Genbank_ID %in% intsct),]
  clnms=str_replace_all(string = colnames(rr1),pattern = "RefSeq_",replacement = "")
  colnames(unq_RefSeq)=clnms
  colnames(unq_GenBank)=clnms
  
  unq_RefSeq$typ_st_end_str=apply(unq_RefSeq[,c(3:5,7)],1,function(x){str_c(str_trim(x,side = "both"),collapse="_")})
  unq_GenBank$typ_st_end_str=apply(unq_GenBank[,c(3:5,7)],1,function(x){str_c(str_trim(x,side="both"),collapse="_")})
  unq_RefSeq$DB="RefSeq"
  unq_GenBank$DB="GenBank"
  
  OvGff=rbind(unq_RefSeq,unq_GenBank)
  OvGff$seqid=str_replace_all(string = OvGff$seqid,pattern = "NZ_",replacement = "")
  unique_typID=unique(OvGff$typ_st_end_str)
  
  for(i in 1:length(unique_typID))
  {
    tmptypID=which(OvGff$typ_st_end_str==unique_typID[i])
    
    IntsctGff=data.frame(seqid=str_c(unique(OvGff$seqid[tmptypID]),collapse ="&"),source=str_c(unique(OvGff$source[tmptypID]),collapse="&"),type=str_c(unique(OvGff$type[tmptypID]),collapse="&"))
    IntsctGff$start=min(unlist(unique(OvGff[tmptypID,c("start","end")])))
    IntsctGff$end=max(unlist(unique(OvGff[tmptypID,c("start","end")])))
    IntsctGff$score=str_c(unique(OvGff$score[tmptypID]),collapse="&")
    IntsctGff$strand=OvGff$strand[tmptypID[1]]
    IntsctGff$phase=str_c(unique(OvGff$phase[tmptypID]),collapse="&")
    
    attrList=(lapply(OvGff$attributes[tmptypID],function(x){unlist(str_split(string = x,pattern = ";"))}))
    names(attrList)=OvGff$DB[tmptypID]
    
    attrVec=lst2df(attrList)
    attrVec2=str_c("Start_End_strand=",str_trim(unique_typID[i],side= "both"))
    attrVec3=str_replace_all(string = OvGff$Locus_tag[tmptypID],pattern = "^ID=",replacement = "")
    attrVec3=str_c("Locus_tag=",str_c(str_trim(unique(attrVec3),side = "both"),collapse="&"))
    ovAttr=str_c(attrVec,";",attrVec2,";",attrVec3)
    IntsctGff$attributes=ovAttr
    if(i==1)
    {
      FinalGff=IntsctGff
      
    } else {
      
      FinalGff=rbind(FinalGff,IntsctGff)
    }
    
    
  }
  
  refseqFin=unq_RefSeq[,c("seqid","start","end","strand","Locus_tag","ID")]
  refseqFin=refseqFin[!duplicated(refseqFin),]
  
  
  GenbankFin=unq_GenBank[,c("seqid","start","end","strand","Locus_tag","ID")]
  GenbankFin=GenbankFin[!duplicated(GenbankFin),]
  
  ### FinalGff: It contains information for those genes commonly present in both genbank and refseq with exactly same genome coordinates.
  
  ###################################################################################################
  ##### Category-2: Those with different ID, Joining all the rows with same Locus_tag
  ###################################################################################################
  
  
  stdfrr1=rr1[(rr1$RefSeq_ID %in% stdff1),]
  stdfrr2=rr2[(rr2$Genbank_ID %in% stdff2),]
  
  stdfrr1$RefSeq_GeneID=str_replace_all(string = stdfrr1$RefSeq_attributes,pattern = ".*;locus_tag=",replacement = "")
  stdfrr1$RefSeq_GeneID=str_replace_all(string = stdfrr1$RefSeq_GeneID,pattern = ";.*",replacement = "")
  stdfrr1$Refseq_length=abs(stdfrr1$RefSeq_start-stdfrr1$RefSeq_end)+1
  
  stdfrr2$Genbank_GeneID=str_replace_all(string = stdfrr2$Genbank_attributes,pattern = ".*;locus_tag=",replacement = "")
  stdfrr2$Genbank_GeneID=str_replace_all(string = stdfrr2$Genbank_GeneID,pattern = ";.*",replacement = "")
  stdfrr2$Genbank_length=abs(stdfrr2$Genbank_start-stdfrr2$Genbank_end)+1
  
  unq_Locus_tags=unique(c(stdfrr1$Locus_tag,stdfrr2$Locus_tag))
  
  for(unqLoc1 in 1: length(unq_Locus_tags))
  {
    
    tmpLocusDf1=stdfrr1[which(stdfrr1$Locus_tag==unq_Locus_tags[unqLoc1]),]
    tmpLocusDf2=stdfrr2[which(stdfrr2$Locus_tag==unq_Locus_tags[unqLoc1]),]
    
    if(dim(tmpLocusDf2)[1]==0)
    {
      tmp_merge_std=tmpLocusDf1  
      
    }
    if(dim(tmpLocusDf1)[1]==0)
    {
      tmp_merge_std=tmpLocusDf2
      
      
    }
    
    if(dim(tmpLocusDf1)[1] !=0 &&  dim(tmpLocusDf2)[1]!=0)
    {
      tmp_merge_std=cbind(tmpLocusDf1,tmpLocusDf2)
      
       
      if(any(c(as.character(tmp_merge_std$RefSeq_type),as.character(tmp_merge_std$Genbank_type))=="pseudogene"))
      {
        if(any(tmpLocusDf1$RefSeq_type=="pseudogene") )
        {
          tmp_merge_std$FinalLocus_tag=tmp_merge_std$Genbank_GeneID
          tmp_merge_std=tmp_merge_std[,-which(colnames(tmp_merge_std)=="Locus_tag")]
          colnames(tmp_merge_std)[grep(pattern = "Final_Locus_tag",x=colnames(tmp_merge_std))]="Locus_tag"
          
        }
        if(any(tmpLocusDf2$Genbank_type=="pseudogene") )
        {
          
          tmp_merge_std$Final_Locus_tag=tmp_merge_std$RefSeq_GeneID
          tmp_merge_std=tmp_merge_std[,-which(colnames(tmp_merge_std)=="Locus_tag")]
          colnames(tmp_merge_std)[grep(pattern = "Final_Locus_tag",x=colnames(tmp_merge_std))]="Locus_tag"
          
          
        }
        
      } else {
        
        if(unique(tmpLocusDf1$Refseq_length) > unique(tmpLocusDf2$Genbank_length))
        {
             tmp_merge_std$Final_Locus_tag=unique(tmpLocusDf1$RefSeq_GeneID)[1]
             tmp_merge_std=tmp_merge_std[,-which(colnames(tmp_merge_std)=="Locus_tag")]
             colnames(tmp_merge_std)[grep(pattern = "Final_Locus_tag",x=colnames(tmp_merge_std))]="Locus_tag"
             
             
        }
      
        
        if(unique(tmpLocusDf2$Genbank_length) > unique(tmpLocusDf1$Refseq_length))
        {
      
          tmp_merge_std$Final_Locus_tag=unique(tmpLocusDf2$Genbank_GeneID)[1]          
          tmp_merge_std=tmp_merge_std[,-which(colnames(tmp_merge_std)=="Locus_tag")]
          colnames(tmp_merge_std)[grep(pattern = "Final_Locus_tag",x=colnames(tmp_merge_std))]="Locus_tag"
              
        }
        
        if(unique(tmpLocusDf1$Refseq_length) == unique(tmpLocusDf2$Genbank_length))
        {
          tmp_merge_std$Final_Locus_tag=unique(tmpLocusDf1$Genbank_GeneID)[1]
          tmp_merge_std=tmp_merge_std[,-which(colnames(tmp_merge_std)=="Locus_tag")]
          colnames(tmp_merge_std)[grep(pattern = "Final_Locus_tag",x=colnames(tmp_merge_std))]="Locus_tag"
          
          
        }
        
        
        }
      
    }
    
    if(unqLoc1 ==1) 
    {
      
      merge_std=tmp_merge_std
     
      
      
    } else {
      
      setDiff_cols=setdiff(colnames(merge_std),colnames(tmp_merge_std))
      
      if(length(setDiff_cols)> 0)
      {
        leftDf=data.frame(matrix(nrow=nrow(tmp_merge_std),ncol=length(setDiff_cols)))
        colnames(leftDf)=setDiff_cols
        tmp_merge_std=cbind(tmp_merge_std,leftDf)
      }
      
      setDiff_cols_merge=setdiff(colnames(tmp_merge_std),colnames(merge_std))
      
      if(length(setDiff_cols_merge)> 0)
      {
        leftDf=data.frame(matrix(nrow=nrow(merge_std),ncol=length(setDiff_cols_merge)))
        colnames(leftDf)=setDiff_cols_merge
        merge_std=cbind(merge_std,leftDf)
      }
      
      
      merge_std=rbind(tmp_merge_std,merge_std)
    }
    
    }
  
  # merge_std=(merge(stdfrr1,stdfrr2,by="Locus_tag"))
  merge_std=merge_std[!duplicated(merge_std),]
  
  na_merge_std=unname(which(apply(merge_std,1,function(x){all(is.na(x))})==T))
  if(length(na_merge_std)> 0 )
  {
    merge_std=merge_std[-na_merge_std,]
    
  }
  
  Cat2unq_RefSeq=merge_std[,grep(pattern = "RefSeq|Locus_tag",x=colnames(merge_std))]
  # Cat2unq_RefSeq$Locus_tag=merge_std$Locus_tag
  na_cat2unq_refseq=which(unname(apply(Cat2unq_RefSeq,1,function(x){
    all(is.na(x))
  }))==TRUE)
  
  if(length(na_cat2unq_refseq) > 0)
  {
    Cat2unq_RefSeq=Cat2unq_RefSeq[-(na_cat2unq_refseq),]
    
  }
  Cat2unq_RefSeq$typ_st_end_str=apply(Cat2unq_RefSeq[,c("RefSeq_type","RefSeq_GeneID")],1,function(x){str_c(str_trim(x,side = "both"),collapse="_")})
  
  Cat2unq_Genbank=merge_std[,grep(pattern = "Genbank|Locus_tag",x=colnames(merge_std))]
  # Cat2unq_Genbank$Locus_tag=merge_std$Locus_tag
  na_cat2unq_Genbank=which(unname(apply(Cat2unq_Genbank,1,function(x){
    all(is.na(x))
  }))==TRUE)
  
  if(length(na_cat2unq_Genbank) > 0)
  {
    Cat2unq_Genbank=Cat2unq_Genbank[-(na_cat2unq_Genbank),]
    
  }
  Cat2unq_Genbank$typ_st_end_str=apply(Cat2unq_Genbank[,c("Genbank_type","Genbank_GeneID")],1,function(x){str_c(str_trim(x,side = "both"),collapse="_")})
  
  
  # Cat2unq_GenBank=merge_std[,grep(pattern = "Genbank",x=colnames(merge_std))]
  # Cat2unq_GenBank$Locus_tag=merge_std$Locus_tag
  # Cat2unq_GenBank$typ_st_end_str=apply(Cat2unq_GenBank[,c("Genbank_type","Locus_tag")],1,function(x){str_c(str_trim(x,side="both"),collapse="_")})
  # 
  Cat2unq_RefSeq$DB="RefSeq"
  Cat2unq_Genbank$DB="GenBank"
  colnames(Cat2unq_Genbank)=str_replace_all(string = colnames(Cat2unq_Genbank),pattern = "Genbank_",replacement = "")
  colnames(Cat2unq_RefSeq)=str_replace_all(string = colnames(Cat2unq_RefSeq),pattern = "RefSeq_",replacement = "")
  cmmn_cols=intersect(colnames(Cat2unq_Genbank),colnames(Cat2unq_RefSeq))
  Cat2OvGff=rbind(Cat2unq_RefSeq[,cmmn_cols],Cat2unq_Genbank[,cmmn_cols])
  Cat2OvGff=Cat2OvGff[!is.na(Cat2OvGff$typ_st_end_str),]
  unique_typID=unique(as.character(Cat2OvGff$typ_st_end_str))
  Cat2OvGff$score[is.na(Cat2OvGff$score)]=""
  Cat2OvGff$strand[is.na(Cat2OvGff$strand)]=""
  Cat2OvGff$phase[is.na(Cat2OvGff$phase)]=""
  
  if(length(unique_typID)>0)
  {
    for(i in 1:length(unique_typID))
    {
      tmptypID=which(Cat2OvGff$typ_st_end_str==unique_typID[i])
      
      Cat2IntsctGff=data.frame(seqid=str_c(unique(Cat2OvGff$seqid[tmptypID]),collapse ="&"),source=str_c(unique(Cat2OvGff$source[tmptypID]),collapse="&"),type=str_c(unique(Cat2OvGff$type[tmptypID]),collapse="&"))
      Lngs=sort(apply(Cat2OvGff[tmptypID,c("start","end" )],1,function(x){
        x=as.numeric(x);
        abs(x[1]-x[2])+1}),decreasing = T)
      
      Cat2IntsctGff$score=str_c(unique(Cat2OvGff$score[tmptypID]),collapse="&")
      Cat2IntsctGff$strand=Cat2OvGff$strand[tmptypID[1]]
      Cat2IntsctGff$phase=str_c(unique(Cat2OvGff$phase[tmptypID]),collapse="&")
      Cat2IntsctGff$start=min(unlist(unique(Cat2OvGff[tmptypID,c("start","end")])))
      Cat2IntsctGff$end=max(unlist(unique(Cat2OvGff[tmptypID,c("start","end")])))
      
      
      attrList=(lapply(Cat2OvGff$attributes[tmptypID],function(x){unlist(str_split(string = x,pattern = ";"))}))
      names(attrList)=Cat2OvGff$DB[tmptypID]
      
      attrVec=lst2df(attrList)
      # attrVec2=str_c("Start_End_strand=",str_trim(unique_typID[i],side= "both"))
      attrVec2=str_c("Start_End_strand=",str_c(as.character(unique(Cat2IntsctGff$type)),"_",min(Cat2IntsctGff$start),"_",max(Cat2IntsctGff$end),"_",unique(Cat2IntsctGff$strand)))
      
      attrVec3=str_replace_all(string = Cat2OvGff$Locus_tag[tmptypID],pattern = "^ID=",replacement = "")
      attrVec3=str_c("Locus_tag=",str_c(str_trim(unique(attrVec3),side = "both"),collapse="&"))
      
      attrTmp=unique(apply(Cat2OvGff[tmptypID,c(2,4,5)],1,function(x){str_c(x[1],"_start=",x[2],";",x[1],"_end=",x[3])}))
      attrVec4=paste(attrTmp,collapse=";")
      
      ovAttr=str_c(attrVec,";",attrVec2,";",attrVec3,";",attrVec4)
      Cat2IntsctGff$attributes=ovAttr
      if(i==1)
      {
        FinalGff2=Cat2IntsctGff
        
      } else {
        
        FinalGff2=rbind(FinalGff2,Cat2IntsctGff)
      }
      
      
    }
    
  } else {
    
    FinalGff2=data.frame()
  }
  
  refseqFin=rbind(refseqFin,Cat2unq_RefSeq[,colnames(refseqFin)])
  GenbankFin=rbind(GenbankFin,Cat2unq_RefSeq[,colnames(GenbankFin)])
  
  
  ###############################################################################################################
  ##### Category-3: Those with different ID and Locus_tags, Joining all the rows subset of another gene : sbstrr
  ##### Category-4: Those with different ID and Locus_tags, Joining all the rows where genes overlap : overlapGn
  ###############################################################################################################
  
  unq_std=setdiff(unique(c(stdff1,stdff2)),unique(c(merge_std$RefSeq_ID, merge_std$Genbank_ID)))
  unq_std_locus=c(stdfrr1$Locus_tag[match(unq_std,stdfrr1$RefSeq_ID)],
                  stdfrr2$Locus_tag[match(unq_std,stdfrr2$Genbank_ID)])
  names(unq_std_locus)=c(stdfrr1$RefSeq_ID[match(unq_std,stdfrr1$RefSeq_ID)],
                         stdfrr2$Genbank_ID[match(unq_std,stdfrr2$Genbank_ID)])
  
  unq_std_locus=unq_std_locus[!is.na(unq_std_locus)]
  
  
  rr1_unq_std=rr1[(rr1$RefSeq_ID %in% unq_std),]
  rr2_unq_std=rr2[(rr2$Genbank_ID %in% unq_std),]
  
  rr1_unq_std$GnCoord=apply(rr1_unq_std[,c("RefSeq_start", "RefSeq_end")],1,function(x){str_c(as.character(seq(x[1],x[2],by=1)),collapse = ",")})
  rr2_unq_std$GnCoord=apply(rr2_unq_std[,c("Genbank_start", "Genbank_end")],1,function(x){str_c(as.character(seq(x[1],x[2],by=1)),collapse = ",")})
  
  
  sbstrr1= vector("list", length = nrow(rr1_unq_std))
  sbstrr2=vector("list", length = nrow(rr2_unq_std))
  names(sbstrr1)=rr1_unq_std$Locus_tag
  names(sbstrr2)=rr2_unq_std$Locus_tag
  
  overlapGn=list()
  k=1
  if(nrow(rr1_unq_std)==0 && nrow(rr2_unq_std)!=0)
  {
    colnames(rr2_unq_std)=str_replace_all(string = colnames(rr2_unq_std),pattern = "RefSeq_|Genbank_",replacement = "")
    rr2_unq_std$attributes=str_c(rr2_unq_std$attributes,";Locus_tag=",rr2_unq_std$Locus_tag)
    FinalGff2=rbind(FinalGff2,rr2_unq_std[,1:ncol(FinalGff2)])
    FinalGff3=data.frame()
    FinalGff4=data.frame()
  } else {
    
    if(nrow(rr2_unq_std)==0 && nrow(rr1_unq_std) !=0)
    {
      colnames(rr1_unq_std)=str_replace_all(string = colnames(rr1_unq_std),pattern = "RefSeq_|Genbank_",replacement = "")
      rr1_unq_std$attributes=str_c(rr1_unq_std$attributes,";Locus_tag=",rr1_unq_std$Locus_tag)
      FinalGff2=rbind(FinalGff2,rr1_unq_std[,1:ncol(FinalGff2)])
      FinalGff3=data.frame()  
      FinalGff4=data.frame()
    } else {
      if(nrow(rr1_unq_std) !=0 && nrow(rr2_unq_std) != 0)
      {
        
        for(i in 1:nrow(rr1_unq_std))
        {
          gnCrds=as.numeric(unlist(str_split(rr1_unq_std$GnCoord[i],pattern = ",")  ))
          
          
          for(j in 1:nrow(rr2_unq_std))
          {
            
            gnCrds2=as.numeric(unlist(str_split(rr2_unq_std$GnCoord[j],pattern = ",")  ))
            instCrds=intersect(gnCrds,gnCrds2)  
            
            if(length(instCrds)==length(gnCrds) && rr1_unq_std$RefSeq_strand[i] == rr2_unq_std$Genbank_strand[j])
            {
              
              sbstrr1[rr1_unq_std$Locus_tag[i]]=c(unname(unlist(sbstrr1[rr1_unq_std$Locus_tag[i]])),rr2_unq_std$Locus_tag[j])
            }
            
            if(length(instCrds)==length(gnCrds2) && rr1_unq_std$RefSeq_strand[i] == rr2_unq_std$Genbank_strand[j] )
            {
              lctags=unname(unlist(sbstrr2[rr2_unq_std$Locus_tag[j]]))
              lctags[is.null(lctags)]=""
              lctags[is.na(lctags)]=""
              lctags_all=unique(c(lctags,rr1_unq_std$Locus_tag[i]))
              lctags_all=lctags_all[lctags_all!=""]
              sbstrr2[rr2_unq_std$Locus_tag[j]]=lctags_all
            }
            
            if(length(instCrds)>1 && rr1_unq_std$RefSeq_strand[i] == rr2_unq_std$Genbank_strand[j])
            {
              if(length(instCrds) != length(gnCrds) && length(instCrds)!=length(gnCrds2) )
              {
                overlapGn[[k]]=c(rr1_unq_std$Locus_tag[i],rr2_unq_std$Locus_tag[j])
                k=k+1  
              }
              
            }
          }
        }
        
        
        ## Remove NULL lists, then remove redundant lists.
        overlapGn=Filter(Negate(is.null), overlapGn)
        
        res.list <- overlapGn
        # take set difference between contents of list elements and accumulated elements
        res.list[-1] <- mapply("setdiff", res.list[-1],
                               head(Reduce(c, overlapGn, accumulate=TRUE), -1))
        
        overlapGn=res.list[lapply(res.list,length)>0]
        
        
        sbstrr1=Filter(Negate(is.null),sbstrr1)
        res.list <- sbstrr1
        
        # take set difference between contents of list elements and accumulated elements
        res.list[-1] <- mapply("setdiff", res.list[-1],
                               head(Reduce(c, sbstrr1, accumulate=TRUE), -1))
        
        sbstrr1=res.list[lapply(res.list,length)>0]
        
        sbstrr2=Filter(Negate(is.null),sbstrr2)
        res.list <- sbstrr2
        # take set difference between contents of list elements and accumulated elements
        res.list[-1] <- mapply("setdiff", res.list[-1],
                               head(Reduce(c, sbstrr2, accumulate=TRUE), -1))
        sbstrr2=res.list[lapply(res.list,length)>0]
        sbstrr=c(sbstrr1,sbstrr2)
        
        ### 1. Create gff file for these overlapping and subset genes:
        ### 2. Add Rows where the overlapping region and non-overlapping region have separate rows.
        
        if(length(sbstrr)>0)
        {
          for(i in 1: length(sbstrr))
          {
            
            tmpG1rr1=rr1_unq_std[( rr1_unq_std$Locus_tag %in% unname(sbstrr)[i]),]
            tmpG1rr2=rr2_unq_std[( rr2_unq_std$Locus_tag %in% unname(sbstrr)[i]),]
            
            tmpG2rr1=rr1_unq_std[( rr1_unq_std$Locus_tag %in% names(sbstrr)[i]),]
            tmpG2rr2=rr2_unq_std[( rr2_unq_std$Locus_tag %in% names(sbstrr)[i]),]
            
            colnames(tmpG1rr1)=str_replace_all(string = colnames(tmpG1rr1),pattern = "RefSeq_|Genbank_",replacement = "")
            colnames(tmpG1rr2)=str_replace_all(string = colnames(tmpG1rr2),pattern = "RefSeq_|Genbank_",replacement = "")
            colnames(tmpG2rr1)=str_replace_all(string = colnames(tmpG2rr1),pattern = "RefSeq_|Genbank_",replacement = "")
            colnames(tmpG2rr2)=str_replace_all(string = colnames(tmpG2rr2),pattern = "RefSeq_|Genbank_",replacement = "")
            
            tmpG1rr1$Length=apply(tmpG1rr1[,c(4,5)],1,function(x){abs(x[1]-x[2])+1})
            tmpG1rr2$Length=apply(tmpG1rr2[,c(4,5)],1,function(x){abs(x[1]-x[2])+1})
            
            tmpG2rr1$Length=apply(tmpG2rr1[,c(4,5)],1,function(x){abs(x[1]-x[2])+1})
            tmpG2rr2$Length=apply(tmpG2rr2[,c(4,5)],1,function(x){abs(x[1]-x[2])+1})
            
            gh=rbind(tmpG1rr1,tmpG1rr2,tmpG2rr1,tmpG2rr2)
            gh=gh[order(gh$Length),]
            gh$Set=""
            if(length(unique(gh$Locus_tag))==1)
            {
              ghGff3=gh 
            } else { 
              bgLoc=unique(gh$Locus_tag[which(gh$Length==max(gh$Length))])
              smLoc=unique(gh$Locus_tag[which(gh$Length!=max(gh$Length))])
              
              for(jj in 1:length(smLoc))
              {
                smGnCrd=gh$GnCoord[which(gh$Locus_tag==smLoc[jj])]
                bgGnCrd=gh$GnCoord[which(gh$Locus_tag==bgLoc)]
                
                Newbg1=gh[gh$Locus_tag==bgLoc,]
                Newbg1$Set="Overlapping"
                
                Newbg2=gh[gh$Locus_tag==bgLoc,]
                Newbg2$Set="Non_overlapping"
                
                Newbg=rbind(Newbg1,Newbg2)
                
                intsctGnCrd=intersect(as.numeric(unlist(str_split(string = bgGnCrd,pattern = ","))),
                                      as.numeric(unlist(str_split(string = smGnCrd,pattern = ","))))
                
                stdffsctGnCrd=setdiff(as.numeric(unlist(str_split(string = bgGnCrd,pattern = ","))),
                                      as.numeric(unlist(str_split(string = smGnCrd,pattern = ","))))
                
                if(all(Newbg$strand=="+"))
                {
                  
                  Newbg1$start=min(intsctGnCrd)
                  Newbg1$end=max(intsctGnCrd)
                  
                  Newbg2$start=min(stdffsctGnCrd)
                  Newbg2$end=max(stdffsctGnCrd)
                  
                  Newbg=rbind(Newbg1,Newbg2)
                  
                  
                } 
                
                if(all(Newbg$strand=="-"))
                {
                  
                  
                  Newbg1$start=max(intsctGnCrd)
                  Newbg1$end=min(intsctGnCrd)
                  
                  Newbg2$start=max(stdffsctGnCrd)
                  Newbg2$end=min(stdffsctGnCrd)
                  
                  Newbg=rbind(Newbg1,Newbg2)
                  
                  
                  
                  
                  
                }
                
                if(jj==1)
                {
                  ExtraRws=Newbg
                  
                } else {
                  
                  ExtraRws=rbind(ExtraRws,Newbg)
                }
                
              }
              
              ghGff3=rbind(gh,ExtraRws)
            }
            if(i==1)
            {
              FinGff3=ghGff3
            } else {
              
              FinGff3=rbind(FinGff3,ghGff3)
            }
          }
          
          FinalGff3=FinGff3[,1:9]
          
          for(i in 1:nrow(FinGff3))
          {
            tmpLocus=ifelse(FinGff3$Set[i]!="",str_c(FinGff3$Locus_tag[i],"_",FinGff3$Set[i]  ),FinGff3$Locus_tag[i])
            FinalGff3$attributes[i]=str_c(FinalGff3$attributes[i],";Locus_tag=",tmpLocus)
            
          }
          FinalGff3=FinalGff3[!duplicated(FinalGff3),]
          
        } else {
          
          FinalGff3=data.frame()
        }
        
        if(length(overlapGn)>0)
        {
          for(i in 1: length(overlapGn))
          {
            tmpGenes=unlist(overlapGn[[i]])
            
            tmpG1rr1=rr1_unq_std[( rr1_unq_std$Locus_tag %in% tmpGenes[1]),]
            tmpG1rr2=rr2_unq_std[( rr2_unq_std$Locus_tag %in% tmpGenes[2]),]
            
            tmpG2rr1=rr1_unq_std[( rr1_unq_std$Locus_tag %in% tmpGenes[1]),]
            tmpG2rr2=rr2_unq_std[( rr2_unq_std$Locus_tag %in% tmpGenes[2]),]
            
            colnames(tmpG1rr1)=str_replace_all(string = colnames(tmpG1rr1),pattern = "RefSeq_|Genbank_",replacement = "")
            colnames(tmpG1rr2)=str_replace_all(string = colnames(tmpG1rr2),pattern = "RefSeq_|Genbank_",replacement = "")
            colnames(tmpG2rr1)=str_replace_all(string = colnames(tmpG2rr1),pattern = "RefSeq_|Genbank_",replacement = "")
            colnames(tmpG2rr2)=str_replace_all(string = colnames(tmpG2rr2),pattern = "RefSeq_|Genbank_",replacement = "")
            
            
            gh=rbind(tmpG1rr1,tmpG1rr2,tmpG2rr1,tmpG2rr2)
            gh=gh[!duplicated(gh),]
            gh$Set=""
            
            if(length(unique(gh$Locus_tag))==1)
            {
              
              ghGff4=gh
            } else { 
              
              smLoc=unique(gh$Locus_tag)[1]
              bgLoc=unique(gh$Locus_tag)[2]
              
              smGnCrd=unique(gh$GnCoord[which(gh$Locus_tag==smLoc)])
              bgGnCrd=unique(gh$GnCoord[which(gh$Locus_tag==bgLoc)])
              
              Newbg1=gh[gh$Locus_tag==bgLoc,]
              Newbg1$Set="Overlapping"
              
              Newbg2=gh[gh$Locus_tag==bgLoc,]
              Newbg2$Set="Non_overlapping"
              
              
              Newbg3=gh[gh$Locus_tag==smLoc,]
              Newbg3$Set="Overlapping"
              
              Newbg4=gh[gh$Locus_tag==smLoc,]
              Newbg4$Set="Non_overlapping"
              
              Newbg=rbind(Newbg1,Newbg2,Newbg3,Newbg4)
              
              intsctGnCrd=intersect(as.numeric(unlist(str_split(string = bgGnCrd,pattern = ","))),
                                    as.numeric(unlist(str_split(string = smGnCrd,pattern = ","))))
              
              stdffsctGnCrdbg=setdiff(as.numeric(unlist(str_split(string = bgGnCrd,pattern = ","))),
                                      as.numeric(unlist(str_split(string = smGnCrd,pattern = ","))))
              
              stdffsctGnCrdsm=setdiff(as.numeric(unlist(str_split(string = smGnCrd,pattern = ","))),
                                      as.numeric(unlist(str_split(string = bgGnCrd,pattern = ","))))
              
              if(all(Newbg$strand=="+"))
              {
                
                Newbg1$start=min(intsctGnCrd)
                Newbg1$end=max(intsctGnCrd)
                
                Newbg2$start=min(stdffsctGnCrdbg)
                Newbg2$end=max(stdffsctGnCrdbg)
                
                
                Newbg3$start=min(intsctGnCrd)
                Newbg3$end=max(intsctGnCrd)
                
                Newbg4$start=min(stdffsctGnCrdsm)
                Newbg4$end=max(stdffsctGnCrdsm)
                
                
                Newbg=rbind(Newbg1,Newbg2,Newbg3,Newbg4)
                
                
              } 
              
              if(all(Newbg$strand=="-"))
              {
                
                Newbg1$start=min(intsctGnCrd)
                Newbg1$end=max(intsctGnCrd)
                
                Newbg2$start=min(stdffsctGnCrdbg)
                Newbg2$end=max(stdffsctGnCrdbg)
                
                Newbg3$start=min(intsctGnCrd)
                Newbg3$end=max(intsctGnCrd)
                
                Newbg4$start=min(stdffsctGnCrdsm)
                Newbg4$end=max(stdffsctGnCrdsm)
                
                Newbg=rbind(Newbg1,Newbg2,Newbg3,Newbg4)
                
              }
              
              ghGff4=rbind(gh,Newbg)
            }  
            if(i==1)
            {
              FinGff4=ghGff4
            } else {
              
              FinGff4=rbind(FinGff4,ghGff4)
            }
          }
          
          FinalGff4=FinGff4[,1:9]
          
          for(i in 1:nrow(FinGff4))
          {
            tmpLocus=ifelse(FinGff4$Set[i]!="",str_c(FinGff4$Locus_tag[i],"_",FinGff4$Set[i]  ),FinGff4$Locus_tag[i])
            FinalGff4$attributes[i]=str_c(FinalGff4$attributes[i],";Locus_tag=",tmpLocus)
            
          }
          
        } else {
          FinalGff4=data.frame()
          
        }
        
        
        
      } else {
        
        FinalGff3=data.frame()
        FinalGff4=data.frame()
        sbstrr=list()
      }
    }
  }
  colnames(rr1_unq_std)=str_replace_all(string = colnames(rr1_unq_std),pattern = "RefSeq_",replacement = "")
  colnames(rr2_unq_std)=str_replace_all(string = colnames(rr2_unq_std),pattern = "Genbank_",replacement = "")
  refseqFin=rbind(refseqFin,rr1_unq_std[,colnames(refseqFin)])
  GenbankFin=rbind(GenbankFin,rr2_unq_std[,colnames(GenbankFin)])
  
  
  #######################################################################################################################
  ### There are few geneIDs with start and stop information missing.
  #### In order to fix that, grep whether a particular locus tag in "GenbankFin" and "refseqFin" variable is present in FinalGff[.234] variables, if so grab the corresponding
  #### start and stop positions to the "GenbankFin" and "refseqFin" variables.
  #######################################################################################################################
  
  GenbankFin$seqid=str_replace_all(string = GenbankFin$seqid,pattern = "NZ_",replacement = "")
  refseqFin$seqid=str_replace_all(string = refseqFin$seqid,pattern = "NZ_",replacement = "")
  
  tmpFinGff=FinalGff
  tmpFinGff2=FinalGff2
  tmpFinGff3=FinalGff3
  tmpFinGff4=FinalGff4
  
  if(dim(tmpFinGff)[1] > 0)
  {
    tmpFinGff$var="FinalGff"
  }
  
  if(dim(tmpFinGff2)[1] > 0)
  {
    tmpFinGff2$var="FinalGff2"
  }
  
  if(dim(tmpFinGff3)[1] > 0)
  {
    tmpFinGff3$var="FinalGff3"
  }
  
  if(dim(tmpFinGff4)[1] > 0)
  {
    tmpFinGff4$var="FinalGff4"
  }
  
  
  tmpOverallFinGff=rbind(tmpFinGff,tmpFinGff2,tmpFinGff3,tmpFinGff4)
  
  na_gbk_lc=unique(c(which(GenbankFin$Locus_tag==""),which(is.na(GenbankFin$Locus_tag))))
  if(length(na_gbk_lc) > 0)
  {
    GenbankFin=GenbankFin[-(na_gbk_lc),]
    
  }
  
  
  na_rfsq_lc=unique(c(which(refseqFin$Locus_tag==""),which(is.na(refseqFin$Locus_tag))))
  if(length(na_rfsq_lc) > 0)
  {
    refseqFin=refseqFin[-(na_rfsq_lc),]
    
  }
  
  na_genbank_coord=unique(c(which(GenbankFin$start==""),which(is.na(GenbankFin$start))))
  na_refseq_coord=unique(c(which(refseqFin$start==""),which(is.na(refseqFin$start))))
  
  na_genbank_unqLoc=unique(GenbankFin$Locus_tag[na_genbank_coord])
  na_refseq_unqLoc=unique(refseqFin$Locus_tag[na_refseq_coord])
  
  if(length(na_genbank_unqLoc) > 0)
  {
    for(naGbk in 1: length(na_genbank_unqLoc))
    {
      
      tmp_idx_na_gbk  =grep(pattern = na_genbank_unqLoc[naGbk],tmpOverallFinGff$attributes)
      
      tmp_na_gbk_start= unique(tmpOverallFinGff$start[tmp_idx_na_gbk])
      tmp_na_gbk_end= unique(tmpOverallFinGff$end[tmp_idx_na_gbk])
      tmp_na_gbk_strand=unique(tmpOverallFinGff$strand[tmp_idx_na_gbk])
      tmp_na_gbk_seqid=as.character(unique(tmpOverallFinGff$seqid[tmp_idx_na_gbk]))
      
      GenbankFin$start[intersect(na_genbank_coord,which(GenbankFin$Locus_tag==na_genbank_unqLoc[naGbk]))]=tmp_na_gbk_start
      GenbankFin$end[intersect(na_genbank_coord,which(GenbankFin$Locus_tag==na_genbank_unqLoc[naGbk]))]=tmp_na_gbk_end
      GenbankFin$strand[intersect(na_genbank_coord,which(GenbankFin$Locus_tag==na_genbank_unqLoc[naGbk]))]=tmp_na_gbk_strand
      GenbankFin$seqid[intersect(na_genbank_coord,which(GenbankFin$Locus_tag==na_genbank_unqLoc[naGbk]))]=as.character(tmp_na_gbk_seqid)
    
      }
    
  }
  
  GenbankFin=GenbankFin[(!duplicated(GenbankFin)),]
  
  if((length(na_refseq_unqLoc)) > 0)
  {
    
  
  for(naRfsq in 1: length(na_refseq_unqLoc))
  {
    
    tmp_idx_na_rfsq  =grep(pattern = na_refseq_unqLoc[naRfsq],tmpOverallFinGff$attributes)
    
    tmp_na_rfsq_start= tmpOverallFinGff$start[tmp_idx_na_rfsq]
    tmp_na_rfsq_end= tmpOverallFinGff$end[tmp_idx_na_rfsq]
    tmp_na_rfsq_strand=tmpOverallFinGff$strand[tmp_idx_na_rfsq]
    tmp_na_rfsq_seqid=as.character(unique(tmpOverallFinGff$seqid[tmp_idx_na_rfsq]))
    
    
    refseqFin$start[intersect(na_refseq_coord,which(refseqFin$Locus_tag==na_refseq_unqLoc[naRfsq]))]=tmp_na_rfsq_start
    refseqFin$end[intersect(na_refseq_coord,which(refseqFin$Locus_tag==na_refseq_unqLoc[naRfsq]))]=tmp_na_rfsq_end
    refseqFin$strand[intersect(na_refseq_coord,which(refseqFin$Locus_tag==na_refseq_unqLoc[naRfsq]))]=tmp_na_rfsq_strand
    refseqFin$seqid[intersect(na_refseq_coord,which(refseqFin$Locus_tag==na_refseq_unqLoc[naRfsq]))]=tmp_na_rfsq_seqid
  }
  }
    refseqFin=refseqFin[(!duplicated(refseqFin)),]
  source("S:/Abilash/Marbouty_annotation_refseq_genbank_2022/scripts/R_files/remove_duplicate_genes.R")

    overallFinGff=remove_duplicate_genes(tmpOverallFinGff,GenbankFin,refseqFin)

    #######################################################################################################################
    # l1=lapply(unique(GenbankFin$Locus_tag),function(x){
    #   
    #   xx=grep(x = tmpOverallFinGff$attributes,pattern=x);
    #   yy=as.character(tmpOverallFinGff$source[xx]);
    #   zz=table(yy);
    #   
    #   if(length(xx)> 1 && length(which(zz>1)) > 0)
    #   {
    #     zzNm=names(zz[zz>1]);
    #     
    #     for(zzC in 1: length(zzNm))
    #     {
    #       tmpzzNm=zzNm[zzC];
    #       
    #     
    #       }
    #     which(names(zz)[zz>1])
    #     
    #     if(length() > 1)
    #     {
    #       return(xx);
    #       
    #     }
    #   } 
    # })
    # 
  
    #### For intersect, all the attributes in V9 for each unique protein
  # View(rbind(rr1[(rr1$RefSeq_ID %in% intsct),],rr2[(rr2$Genbank_ID %in% intsct),]))
  
  #### For setdiff, check whether they are subset of another protein:
  
  ###### If yes, check whether they both are with same locus_tag in Genbank and old_locus_tag in RefSeq.
  
  
  ######### If yes, then include only the start and end part of the largest region.
  ######### If no, then they both are different. Retain both of them as separate proteins.
  ###### If no, then they both are different. Retain both of them as separate proteins.
  
  FinalGff=overallFinGff[overallFinGff$var=="FinalGff",]  
  FinalGff=FinalGff[,-grep(pattern = "var",x=colnames(FinalGff))]
  FinalGff2=overallFinGff[overallFinGff$var=="FinalGff2",]  
  FinalGff2=FinalGff2[,-grep(pattern = "var",x=colnames(FinalGff2))]
  FinalGff3=overallFinGff[overallFinGff$var=="FinalGff3",]  
  FinalGff3=FinalGff3[,-grep(pattern = "var",x=colnames(FinalGff3))]
  FinalGff4=overallFinGff[overallFinGff$var=="FinalGff4",]  
  FinalGff4=FinalGff4[,-grep(pattern = "var",x=colnames(FinalGff4))]
  
  if(any(is.na(GenbankFin$ID)))
  {
    GenbankFin$ID[which(is.na(GenbankFin$ID))]=unlist(lapply(which(is.na(GenbankFin$ID)),function(x){
      
      str_c(GenbankFin$start[x],"_",GenbankFin$end[x],"_",GenbankFin$strand[x])
      
    }))
    
  }
  
  
  if(any(is.na(refseqFin$ID)))
  {
    refseqFin$ID[which(is.na(refseqFin$ID))]=unlist(lapply(which(is.na(refseqFin$ID)),function(x){
      
      str_c(refseqFin$start[x],"_",refseqFin$end[x],"_",refseqFin$strand[x])
      
    }))
    
  }
  
    
  return(new(Class = "gffAnnotate",FinalGff=FinalGff,FinalGff2=FinalGff2,RefSeq_coord=refseqFin,Genbank_coord=GenbankFin,
             FinalGff3=FinalGff3,FinalGff4=FinalGff4,overlapGn=overlapGn,sbstrr=sbstrr))  
  
  
}



for(i in 1:length(unqGff))
{
  
  tmpGff=which(ovGff==unqGff[i])
  tmpLnk=as.character(ovLnk[tmpGff])
 
  if(length(tmpLnk)>1)
  {
    # flType=unlist(lapply(tmpLnk,function(x){ifelse(length(grep(pattern = "_genomic",x=x))>0,"Genbank","RefSeq")}))
    flType=srcTyp[tmpGff]
    rr1=read.gff(file = tmpLnk[2])
    rr2=read.gff(file = tmpLnk[1])
    
    rr1$ID=apply(rr1[,c(4,5,7)],1,function(x){str_c(str_trim(x,side = "both"),collapse="_")})
    rr2$ID=apply(rr2[,c(4,5,7)],1,function(x){str_c(str_trim(string = x,side = "both"),collapse="_")})
    
    
    overallgff=unique_gffv2(rr1,rr2,flType)  
    
  
    } else {
    
    src1=unlist(str_split(string = tmpLnk,pattern = "/"))[6]
    src1=str_replace_all(string = src1,pattern = "Marbouty_",replacement = "")
    
    overallgff=read.gff(file = tmpLnk[1])
    if(src1=="Genbank")
    {
      Genbank_coord=overallgff[,c("start","end","strand")]
      Genbank_coord$ID=apply(Genbank_coord[,c(1:3)],1,function(x){
        str_c(x,collapse = "_")
        
      })
    Genbank_coord$Locus_tag=overallgff$attributes
   Genbank_coord$Locus_tag=str_replace_all(string = Genbank_coord$Locus_tag,pattern = ".*locus_tag=",replacement = "")   
   Genbank_coord$Locus_tag=str_replace_all(string = Genbank_coord$Locus_tag,pattern = ";.*",replacement = "")   
   RefSeq_coord=data.frame()
   
   }
     
   if(src1=="RefSeq")
    {
      RefSeq_coord=overallgff[,c("start","end","strand")]
      RefSeq_coord$ID=apply(RefSeq_coord[,c(1:3)],1,function(x){
        str_c(x,collapse = "_")
        
      })
      RefSeq_coord$Locus_tag=overallgff$attributes
      RefSeq_coord$Locus_tag=str_replace_all(string = RefSeq_coord$Locus_tag,pattern = ".*locus_tag=",replacement = "")   
      RefSeq_coord$Locus_tag=str_replace_all(string = RefSeq_coord$Locus_tag,pattern = ";.*",replacement = "")   
      Genbank_coord=data.frame()
      
    }
    
    }
  if(class(overallgff)=="gffAnnotate")
  {
    
  saveRDS(object = overallgff,file = str_c(MarboutyFolder,unqGff[i],"_overallgff.rds"))
  saveRDS(object = overallgff@RefSeq_coord,file = str_c(MarboutyFolder,unqGff[i],"_RefSeq_gene_coordinates.rds"))
  saveRDS(object = overallgff@Genbank_coord,file = str_c(MarboutyFolder,unqGff[i],"_Genbank_gene_coordinates.rds"))
  
  } else {
    saveRDS(object = overallgff,file = str_c(MarboutyFolder,unqGff[i],"_overallgff.rds"))
    saveRDS(object = RefSeq_coord,file = str_c(MarboutyFolder,unqGff[i],"_RefSeq_gene_coordinates.rds"))
    saveRDS(object = Genbank_coord,file = str_c(MarboutyFolder,unqGff[i],"_Genbank_gene_coordinates.rds"))
    
    
}
  }


for(jj in 1:length(unqGff))
{
  
  gnAnnotation=readRDS(str_c(MarboutyFolder,unqGff[jj],"_overallgff.rds"))
  
  if(all(! slotNames(gnAnnotation) %in% ".S3Class") == TRUE)
  {
    
    tmpgff=rbind(gnAnnotation@FinalGff,gnAnnotation@FinalGff2)
    tmpoverlap=rbind(gnAnnotation@FinalGff3,gnAnnotation@FinalGff4)
    
  } else {
    
    tmpgff=cbind(as.character(gnAnnotation@.Data[[1]]),as.character(gnAnnotation@.Data[[2]]),as.character(gnAnnotation@.Data[[3]]),as.character(gnAnnotation@.Data[[4]]),as.character(gnAnnotation@.Data[[5]]),
                 as.character(gnAnnotation@.Data[[6]]), as.character( gnAnnotation@.Data[[7]]),as.character(gnAnnotation@.Data[[8]]),as.character(gnAnnotation@.Data[[9]]))
    
    colnames(tmpgff)=c("seqid"   ,   "source"   ,  "type"    ,   "start"    ,  "end"     ,   "score"    ,  "strand"    , "phase"    , "attributes")
    tmpoverlap=data.frame()
  }
  
    if(jj==1)
  {
    ovGffFile=tmpgff
    overlapGff=tmpoverlap
    
  } else {
    
    ovGffFile=rbind(ovGffFile,tmpgff)
    overlapGff=rbind(overlapGff,tmpoverlap)
  }
  
}

locus_tagIdx=grep(pattern = "locus_tag=",x=ovGffFile$attributes)
Locus_tagIdx=grep(pattern = "Locus_tag=",x=ovGffFile$attributes)
stdffloc_Loc=setdiff(locus_tagIdx,Locus_tagIdx)
if(length(stdffloc_Loc)>0)
{
  
  ovGffFile$attributes[stdffloc_Loc]=str_replace_all(string = ovGffFile$attributes[stdffloc_Loc],pattern = "locus_tag=",replacement = "Locus_tag=")
  
}

locus_tagIdx=grep(pattern = "locus_tag=",x=overlapGff$attributes)
Locus_tagIdx=grep(pattern = "Locus_tag=",x=overlapGff$attributes)
stdffloc_Loc=setdiff(locus_tagIdx,Locus_tagIdx)
if(length(stdffloc_Loc)>0)
{
  overlapGff$attributes[stdffloc_Loc]=str_replace_all(string = overlapGff$attributes[stdffloc_Loc],pattern = "locus_tag=",replacement = "Locus_tag=")
  
}

ovGffFile$seqid=str_replace_all(string = ovGffFile$seqid,pattern = "NZ_",replacement = "")
overlapGff$seqid=str_replace_all(string = overlapGff$seqid,pattern = "NZ_",replacement = "")

write.table(x = ovGffFile,file = str_c(MarboutyFolder,"OMM_ecoli_overall_cat12.gff"),quote = F,sep = "\t",row.names = F,col.names = F)
write.table(x = overlapGff,file = str_c(MarboutyFolder,"OMM_ecoli_overall_cat34.gff"),quote = F,sep = "\t",row.names = F,col.names = F)
write.table(x = rbind(ovGffFile,overlapGff),file = str_c(MarboutyFolder,"OMM_ecoli_overall_genbank_and_refseq.gff"),quote = F,sep = "\t",row.names = F,col.names = F)

cdsRgidx=sort(c(which(ovGffFile$type=="region"),which(ovGffFile$type=="CDS")))
ovGffCDS=ovGffFile[cdsRgidx,]
ovGffCDS$attributes=str_replace_all(string = ovGffCDS$attributes,pattern = ".*;Locus_tag=",replacement = "Locus_tag=")
ovGffCDS$attributes=str_replace_all(string = ovGffCDS$attributes,pattern = ";.*",replacement = "")
ovGffCDS$attributes[which(ovGffCDS$type=="region")]=str_replace_all(string = ovGffCDS$attributes[which(ovGffCDS$type=="region")],pattern = "Locus_tag=",replacement = "ID=")



write.table(x = ovGffCDS,file = str_c(MarboutyFolder,"OMM_ecoli_overall_cat12_onlyCDS.gff"),quote = F,sep = "\t",row.names = F,col.names = F)


## gff file with only Locus tag 
ovGffShort=ovGffFile
# ovGffShort$attributesFull=ovGffShort$attributes
ovGffShort$attributes=str_replace_all(string = ovGffShort$attributes,pattern = ".*;Locus_tag=",replacement = "Locus_tag=")
ovGffShort$attributes=str_replace_all(string = ovGffShort$attributes,pattern = ";.*",replacement = "")
ovGffShort$attributes[which(ovGffShort$type=="region")]=str_replace_all(string = ovGffShort$attributes[which(ovGffShort$type=="region")],pattern = "Locus_tag=",replacement = "ID=")
ovGffShort$seqid=str_replace_all(string = ovGffShort$seqid,pattern = ".*&",replacement = "")
ovGffShort$seqid=str_replace_all(string = ovGffShort$seqid,pattern = "NZ_",replacement = "")
ovGffShort$attributes=str_replace_all(string = ovGffShort$attributes,pattern = "ID=.*&",replacement = "ID=")
ovGffShort$attributes=str_replace_all(string = ovGffShort$attributes,pattern = "NZ_",replacement = "")
ovGffShort$attributes=str_replace_all(string = ovGffShort$attributes,pattern = "Locus_tag=&",replacement = "Locus_tag=")
ovGffShort$attributes=str_replace_all(string = ovGffShort$attributes,pattern = "Locus_tag=&",replacement = "Locus_tag=")
ovGffShort$seqid=str_replace_all(string = ovGffShort$seqid,pattern = "CP028714.1",replacement = "NZ_CP028714.1")
ovGffShort$seqid=str_replace_all(string = ovGffShort$seqid,pattern = "CP028715.1",replacement = "NZ_CP028715.1")
ovGffShort$attributes=str_replace_all(string = ovGffShort$attributes,pattern = "CP028714.1",replacement = "NZ_CP028714.1")
ovGffShort$attributes=str_replace_all(string = ovGffShort$attributes,pattern = "CP028715.1",replacement = "NZ_CP028715.1")
ovGffShort$attributes[which(ovGffShort$type!="region")]=str_replace_all(string = ovGffShort$attributes[which(ovGffShort$type!="region")],pattern = "ID=",replacement = "Locus_tag=")
ovGffShort$attributes[which(ovGffShort$type!="region")]=str_replace_all(string = ovGffShort$attributes[which(ovGffShort$type!="region")],pattern = "&.*",replacement = "")

for(i in 1:nrow(ovGffShort))
{
  ovGffShort[i,4]=c(min(as.numeric(ovGffShort[i,c(4,5)])))
  ovGffShort[i,5]=c(max(as.numeric(ovGffShort[i,c(4,5)])))
  
  
}

write.table(x = ovGffShort,file = str_c(MarboutyFolder,"OMM_ecoli_overall_cat12_short.gff"),quote = F,sep = "\t",row.names = F,col.names = F)

locus_tagIdx=grep(pattern = "locus_tag=",x=overlapGff$attributes)
Locus_tagIdx=grep(pattern = "Locus_tag=",x=overlapGff$attributes)
stdffloc_Loc=setdiff(locus_tagIdx,Locus_tagIdx)
if(length(stdffloc_Loc)>0)
{
  
  overlapGff$attributes[stdffloc_Loc]=str_replace_all(string = overlapGff$attributes[stdffloc_Loc],pattern = "locus_tag=",replacement = "Locus_tag=")
  
}






## gff file with only Locus tag 
overlapShort=overlapGff
overlapShort$attributes=str_replace_all(string = overlapShort$attributes,pattern = ".*;Locus_tag=",replacement = "Locus_tag=")
overlapShort$attributes=str_replace_all(string = overlapShort$attributes,pattern = ";.*",replacement = "")
overlapShort$attributes[which(overlapShort$type=="region")]=str_replace_all(string = overlapShort$attributes[which(overlapShort$type=="region")],pattern = "Locus_tag=",replacement = "ID=")
overlapShort$seqid=str_replace_all(string = overlapShort$seqid,pattern = ".*&",replacement = "")
overlapShort$seqid=str_replace_all(string = overlapShort$seqid,pattern = "NZ_",replacement = "")
overlapShort$attributes=str_replace_all(string = overlapShort$attributes,pattern = "ID=.*&",replacement = "ID=")
overlapShort$attributes=str_replace_all(string = overlapShort$attributes,pattern = "NZ_",replacement = "")
overlapShort$attributes=str_replace_all(string = overlapShort$attributes,pattern = "Locus_tag=&",replacement = "Locus_tag=")
overlapShort$attributes=str_replace_all(string = overlapShort$attributes,pattern = "Locus_tag=&",replacement = "Locus_tag=")
overlapShort$seqid=str_replace_all(string = overlapShort$seqid,pattern = "CP028714.1",replacement = "NZ_CP028714.1")
overlapShort$seqid=str_replace_all(string = overlapShort$seqid,pattern = "CP028715.1",replacement = "NZ_CP028715.1")
overlapShort$attributes=str_replace_all(string = overlapShort$attributes,pattern = "CP028714.1",replacement = "NZ_CP028714.1")
overlapShort$attributes=str_replace_all(string = overlapShort$attributes,pattern = "CP028715.1",replacement = "NZ_CP028715.1")
overlapShort$attributes[which(overlapShort$type!="region")]=str_replace_all(string = overlapShort$attributes[which(overlapShort$type!="region")],pattern = "ID=",replacement = "Locus_tag=")
overlapShort$attributes[which(overlapShort$type!="region")]=str_replace_all(string = overlapShort$attributes[which(overlapShort$type!="region")],pattern = "&.*",replacement = "")

for(i in 1:nrow(overlapShort))
{
  crd1=as.numeric(overlapShort[i,c(4,5)])
  overlapShort[i,4]=min(crd1)
  overlapShort[i,5]=max(crd1)
  
  
}

write.table(x = overlapShort,file = str_c(MarboutyFolder,"OMM_ecoli_overall_cat34_short.gff"),quote = F,sep = "\t",row.names = F,col.names = F)

overallShort=rbind(ovGffShort,overlapShort)
overallShort=overallShort[with(overallShort, order(seqid,as.numeric(start),as.numeric(end),strand)), ]
write.table(x = overallShort,file = str_c(MarboutyFolder,"OMM_ecoli_overall_cat1234_short.gff"),quote = F,sep = "\t",row.names = F,col.names = F)
as.character(unique(overallShort$type))

overallRegion=overallShort[overallShort$type=="region",]
overallgene=overallShort[overallShort$type=="gene",]
overallCDS=overallShort[overallShort$type=="CDS",]

type2=sort(unique(as.character(overallShort$type)))
type2=sort(setdiff(type2,c("CDS","gene","region")))

overallShortFin=overallgene

for(i in 1:length(type2))
{
  
  tmpIdx=which(overallShort$type==type2[i])
  
  IdxInt=intersect(overallShort$attributes[tmpIdx],overallShortFin$attributes)
  IdxStDff=setdiff(overallShort$attributes[tmpIdx],overallShortFin$attributes)
  
  if(length(IdxInt)>0)
  {
    overallShortFin$type[match(IdxInt,overallShortFin$attributes)]=type2[i]
    
    
  }
  
  if(length(IdxStDff)> 0)
  {
    overallShortFin=rbind(overallShortFin,overallShort[match(IdxStDff,overallShort$attributes),])
    
  }
  

  
  }

## find those annotations only in CDS and add them.
tmpIdx=which(overallShort$type=="CDS")
IdxStDff=setdiff(overallShort$attributes[tmpIdx],overallShortFin$attributes)

if(length(IdxStDff)> 0)
{
  overallShortFin=rbind(overallShortFin,overallShort[match(IdxStDff,overallShort$attributes),])
  
}

overallShortFin$attributes=apply(overallShortFin[,c(3,9)],1,function(x){str_c(x[2],";biotype=",x[1])})
overallShortFin$type="GeneFeatures"
overallShortFin=rbind(overallRegion,overallShortFin)

headerFile=readLines(str_c(MarboutyFolder,"OMM_ecoli_overall_headers"))
file.remove(str_c(MarboutyFolder,"OMM_ecoli_overallShortFin.gff"))

sink(str_c(MarboutyFolder,"OMM_ecoli_overallShortFin.gff"))
cat(headerFile,sep = "\n")
sink()
write.table(x = overallShortFin[with(overallShortFin, order(seqid,as.numeric(start),as.numeric(end),strand)), ],file = str_c(MarboutyFolder,"OMM_ecoli_overallShortFin.gff"),append = T,quote = F,sep = "\t",row.names = F,col.names = F)


overallShortFinv2=overallShortFin
rm_Fin=grep(pattern = "verlapping",x=overallShortFinv2$attributes)
if(length(rm_Fin) > 0)
{
  overallShortFinv2=overallShortFinv2[-(rm_Fin),]
  
}

file.remove(str_c(MarboutyFolder,"OMM_ecoli_overallShortNoOvlp.gff"))

sink(str_c(MarboutyFolder,"OMM_ecoli_overallShortNoOvlp.gff"))
cat(headerFile,sep = "\n")
sink()
write.table(x = overallShortFinv2[with(overallShortFinv2, order(seqid,as.numeric(start),as.numeric(end),strand)), ],file = str_c(MarboutyFolder,"OMM_ecoli_overallShortNoOvlp.gff"),append = T,quote = F,sep = "\t",row.names = F,col.names = F)

altogetherGff=rbind(ovGffFile,overlapGff)

altogetherGff$seqid=str_replace_all(string = altogetherGff$seqid,pattern = "NZ_",replacement = "")
altogetherGff$seqid=str_replace_all(string = altogetherGff$seqid,pattern = ".*&",replacement = "")

ID2Name2=gnmID2Org
ID2Name2$OrganismFin=ID2Name2$Organism
ID2Name2$OrganismFin=str_replace_all(string = ID2Name2$OrganismFin,pattern = "Akkermansia",replacement = "akkermansia")
ID2Name2$OrganismFin=str_replace_all(string = ID2Name2$OrganismFin,pattern = "Escherichia",replacement = "escherichia")
ID2Name2$OrganismFin=str_replace_all(string = ID2Name2$OrganismFin,pattern = "muribaculum_intestinale_YL27",replacement = "muribaculum_intestinales_YL27")


for(i in 1:nrow(ID2Name2))
{

  tmpGff=altogetherGff[altogetherGff$seqid==ID2Name2$GenbankID[i],]
    
  flName=str_c("S:/Abilash/Marbouty_annotation_refseq_genbank_2022/gffFiles/",ID2Name2$OrganismFin[i],"_genomic.gff")
  
  file.remove(flName)
  write.table(x = tmpGff[with(tmpGff, order(seqid,as.numeric(start),as.numeric(end),strand)), ],file = flName,append = F,quote = F,sep = "\t",row.names = F,col.names = F)
  
}

#### Added the headers later


### bash script for gff final cat12 files
# rm -f  /dss/dssfs02/lwp-dss-0001/u7b03/u7b03-dss-0000/ra72seh/DB/OMM_ecoli_overall_cat12Fin.gff
# touch /dss/dssfs02/lwp-dss-0001/u7b03/u7b03-dss-0000/ra72seh/DB/OMM_ecoli_overall_cat12Fin.gff
# cat /dss/dssfs02/lwp-dss-0001/u7b03/u7b03-dss-0000/ra72seh/DB/OMM_ecoli_overall_headers  >>  /dss/dssfs02/lwp-dss-0001/u7b03/u7b03-dss-0000/ra72seh/DB/OMM_ecoli_overall_cat12Fin.gff
# cat /dss/dssfs02/lwp-dss-0001/u7b03/u7b03-dss-0000/ra72seh/DB/OMM_ecoli_overall_cat12.gff  >>  /dss/dssfs02/lwp-dss-0001/u7b03/u7b03-dss-0000/ra72seh/DB/OMM_ecoli_overall_cat12Fin.gff




