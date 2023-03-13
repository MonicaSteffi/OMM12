library("stringr")
library("RCurl")
library("data.table")
library("seqinr")

kofamPath="S:/Abilash/Marbouty_annotation_refseq_genbank 11_2021/Data/Marbouty_annotations_for_integrated_proteins/KOFAMSCAN_results/"
gffPath="S:/Abilash/Marbouty_annotation_refseq_genbank 11_2021/Data/Marbouty_annotations_for_integrated_proteins/gffFiles/"
OperonPath="S:/Abilash/Marbouty_annotation_refseq_genbank 11_2021/Data/Marbouty_annotations_for_integrated_proteins/OperonMapper/"
cgcPath="S:/Abilash/Marbouty_annotation_refseq_genbank 11_2021/Data/Marbouty_annotations_for_integrated_proteins/dbCAN2/"
eggNOGPath="S:/Abilash/Marbouty_annotation_refseq_genbank 11_2021/Data/Marbouty_annotations_for_integrated_proteins/eggNOG_2022/"
list_koFiles=list.files(path = kofamPath,full.names = T,pattern = "_result_all_tab_edited")
gnmNms=str_replace_all(string = list_koFiles,pattern = ".*\\/",replacement = "")
gnmNms=str_replace_all(string = gnmNms,pattern = "_result.*",replacement = "")


############################################################################################################################
##### PUL annotation
############################################################################################################################
pulFA=read.fasta(file = "S:/Abilash/Marbouty_annotation_refseq_genbank 11_2021/Data/Marbouty_annotations_for_integrated_proteins/dbCAN-PUL/PUL.faa.txt",seqtype = "AA")
nmPULFA=names(pulFA)
strSp=str_split(string = nmPULFA,pattern = "\t")
for(i in 1:length(strSp))
{
chrc=unlist(strSp[[i]])
  Df=data.frame(PULGeneID=chrc[1],PULID=chrc[2],GeneName=chrc[3],Locus_tagSp=chrc[4],
                ProteinIDSp=chrc[5],DBType=chrc[6],Ann=chrc[7])  
  
  if(i==1)
  {PULAnnDf=Df} else {
    
    PULAnnDf=rbind(PULAnnDf,Df)
  }
}

## Downloaded from https://bcb.unl.edu/dbcan_pul/Webserver/static/DBCAN-PUL/
dbcan_pul=read.csv(file ="S:/Abilash/Marbouty_annotation_refseq_genbank 11_2021/Data/Marbouty_annotations_for_integrated_proteins/dbCAN-PUL/dbCAN-PUL_updated_05_09_2021",header = T,sep = "\t" ) # ,fileEncoding="latin1"

overallPUL=merge(dbcan_pul,PULAnnDf,by="PULID",all=TRUE)
write.table(overallPUL,file = "S:/Abilash/Marbouty_annotation_refseq_genbank 11_2021/Data/Marbouty_annotations_for_integrated_proteins/dbCAN-PUL/overallPUL_annotation.csv",quote = F,sep = "\t",row.names = F)
saveRDS(overallPUL,"S:/Abilash/Marbouty_annotation_refseq_genbank 11_2021/Data/Marbouty_annotations_for_integrated_proteins/dbCAN-PUL/overallPUL_annotation.rds")

for(i2 in 1:length(gnmNms))
{
  
  ## Read KOFAM
  kofam=readRDS(str_c(kofamPath,gnmNms[i2],"_kofam_results_eval1e03.rds"))
  ## Read dbCAN2
  dbCAN2=readRDS(str_c("S:/Abilash/Marbouty_annotation_refseq_genbank 11_2021/Data/Marbouty_annotations_for_integrated_proteins/",gnmNms[i2],"_overall_cgc.rds"))
  ## Read OperonMapper
  OperonfileName=str_c("S:/Abilash/Marbouty_annotation_refseq_genbank 11_2021/Data/Marbouty_annotations_for_integrated_proteins/",gnmNms[i2],"_operonList_cleaned")
  OperonMapper=readRDS(file = str_c(OperonfileName,".rds"))
  OperonMapper$str_Crd=apply(OperonMapper[,c("PosLeft"  , "postRight", "Strand")],1,function(x){str_c(str_trim(x,side = "both"),collapse = "_")})

  
  ## Read PULs
  PULRes=tryCatch(read.csv(str_c("S:/Abilash/Marbouty_annotation_refseq_genbank 11_2021/Data/Marbouty_annotations_for_integrated_proteins/PULs/",gnmNms[i2],"V10_pul.tsv"),header = T,sep = "\t"),error=function(x){NA})
  ## Read gffFiles
  ov_ann_file=str_c("S:/Abilash/Marbouty_annotation_refseq_genbank 11_2021/Data/Marbouty_annotations_for_integrated_proteins/",gnmNms[i2],"_gff_results.txt")
  gffFiles=readRDS(file = str_c(ov_ann_file,".rds"))
  ovlpIdx=grep(pattern = "verlapping",x=gffFiles$locus_tag)
  
  if(length(ovlpIdx)>0)
  {
    gffFiles=gffFiles[-ovlpIdx,]
    
  }
  
  colnames(gffFiles)[which(colnames(gffFiles)=="GeneID")]="RefSeqGeneID"
  colnames(gffFiles)[which(colnames(gffFiles)=="locus_tag")]="GeneID"
  
  ## Join kofam and dbCAN2
  kofam_dbCAN2=merge(kofam,dbCAN2,by="GeneID",all=TRUE)
  
  ## Join kofam_dbCAN2 and gffFiles
  ko_dbCAN_gff=merge(kofam_dbCAN2,gffFiles,by="GeneID",all=TRUE)
  ko_dbCAN_gff$str_Crd=apply(ko_dbCAN_gff[,c("Start.y","Stop","Strand.y")],1,function(x){str_c(str_trim(x,side = "both"),collapse = "_")})
  ko_dbCAN_gff_op=merge(ko_dbCAN_gff,OperonMapper,by="str_Crd",all=TRUE)
  ko_dbCAN_gff_op=ko_dbCAN_gff_op[!is.na(ko_dbCAN_gff_op$str_Crd),]
  
  if(!is.na(PULRes))
  {
    PULRes$str_Crd=apply(PULRes[,c("start","end","strand")],1,function(x){str_c(str_trim(as.character(x),side = "both"),collapse = "_")})
    FinalAnn=merge(PULRes,ko_dbCAN_gff_op,by="str_Crd",all=TRUE)
    FinalAnn=FinalAnn[!is.na(FinalAnn$str_Crd),]
  } else {
    
    FinalAnn=ko_dbCAN_gff_op
  }
  
  FinalAnn$DB[is.na(FinalAnn$DB)]=""
  substrates=rep("",nrow(FinalAnn))
  geneNames=rep("",nrow(FinalAnn))
  PULID=rep("",nrow(FinalAnn))
  GHID=rep("",nrow(FinalAnn))
  
  tmpDescList=lapply(FinalAnn$Description,function(x){str_replace_all(pattern =   "DB=",string = x,replacement = "")})
  tmpDescList=unlist(lapply(tmpDescList,function(x){str_replace_all(string = x,pattern = ";.*",replacement = "")}))
  # tmpDescList=unlist(lapply(tmpDescList,function(x){unlist(str_split(string = x,pattern = "\\|"))}))
  ###### PUL Gene ID to dbCAN ID:
  PUL_subList=as.list(overallPUL$substrate_final)
  names(PUL_subList)=as.character(overallPUL$PULGeneID)
  
  PUL_prd_dbCAN=as.list(as.character(overallPUL$cazymes_predicted_dbcan))
  names(PUL_prd_dbCAN)=as.character(overallPUL$PULGeneID)
  
  PUL_prd_dbCAN_v2=lapply(PUL_prd_dbCAN,function(x){unlist(str_split(unlist(str_split(string = x,pattern = ",")),pattern = "\\|"))})
  unqGH=unique(unlist(PUL_prd_dbCAN_v2))
  
  
  PUL_names=unique(names(PUL_prd_dbCAN_v2))
  
  PUL_GH_membership=lapply(unqGH,function(x){nn=names(which(PUL_prd_dbCAN_v2==x)); nn[!is.na(nn)]})
  names(PUL_GH_membership)=as.character(unqGH)
  
 GH_subList=lapply(1:length(PUL_GH_membership),function(x){
   x1=unlist(PUL_GH_membership[[x]]);
   
   xx=unique(unname(as.character(unlist(PUL_subList[x1]))));
   xx
   })
  names(GH_subList)=names(PUL_GH_membership)
  
  
  for(ii in 1: nrow(FinalAnn))
  {
    if(FinalAnn$DB[ii]=="TC")
    {
      tmpDsc=tmpDescList[ii]
    } else {
      if(FinalAnn$DB[ii]=="CAZyme")
      {
        tmpDsc=unique(unlist(str_split(tmpDescList[ii],pattern = "\\|")))
        
        
      } else {
        tmpDsc=NA
      } 
      
    }
    
    if(!is.na(tmpDsc))
    {
      
    #### Old:  
    # tmpDsc=str_replace_all(string = tmpDsc,pattern = ";.*",replacement = "")  
    # tmpDsc=str_replace_all(string = tmpDsc,pattern = "DB=",replacement = "")  
    # tmpDsc=unlist(str_split(string = tmpDsc,pattern = "\\|"))
    
    # tt1=table(unlist(lapply(tmpDsc,function(x){grep(overallPUL$Ann,pattern=x)})))
    # IdxFin=names(tt1)[which(tt1==length(tmpDsc))]
    
      #### New: 
      
    # if(length(IdxFin)>0)
    # {
    # uuSbs=unique(overallPUL$substrate_final[as.numeric(IdxFin)])
    # uugNms=unique(overallPUL$GeneName[as.numeric(IdxFin)])
    # uuSbs=uuSbs[uuSbs!=""]
    # uugNms=uugNms[uugNms!=""]
    # 
    # substrates[ii]=str_c(uuSbs,collapse = ";")
    # geneNames[ii]=str_c(uugNms,collapse = ";")
    #   
    #   
    # }
    # 
    ############# For substrate
      if(FinalAnn$DB[ii]=="CAZyme")
      {
       mt_nms=match(tmpDsc,names(GH_subList)  )
       mt_nms=mt_nms[!is.na(mt_nms)]
       
       if(length(mt_nms) > 0)
       {
         
         tmpSub=tryCatch(unique(unlist(GH_subList[mt_nms])),error=function(x){""})
         tmpSub=str_c(tmpSub,collapse = ",")
         tmpSub=str_replace_all(string = tmpSub,pattern = ",",replacement = ";")
         substrates[ii]=ifelse(length(tmpSub) > 0, tmpSub, "")
         PULID[ii]=str_c(unname(unlist(PUL_GH_membership[intersect(names(PUL_GH_membership),tmpDsc)])),collapse = ";")
         
         
         }
        
      }
      
      if(FinalAnn$DB[ii]=="TC")
      {
        
        tmpSub=tryCatch(unique(overallPUL$substrate_final[which(as.character(overallPUL$Ann)==tmpDescList[ii])]),error=function(x){""})
        tmpSub=str_replace_all(string = tmpSub,pattern = ",",replacement = ";")
        substrates[ii]=ifelse(length(tmpSub) > 0, tmpSub, "")
        
        PULID[ii]=tryCatch(str_c(overallPUL$PULGeneID[which(as.character(overallPUL$Ann)==tmpDescList[ii])],collapse = ";"),error=function(x){""})
        
      
        }
      ############# For geneNames
      mt_nms=match(tmpDsc,names(GH_subList)  )
      mt_nms=mt_nms[!is.na(mt_nms)]
      
      intVar=intersect(names(GH_subList)[mt_nms],names(PUL_GH_membership))
      
      Gn1=as.character(unlist(overallPUL$GeneName[match(unique(unlist(PUL_GH_membership[intVar])) , overallPUL$PULGeneID )]))
      Gn1=unique(Gn1[Gn1!=""])
      geneNames[ii]=str_c(sort(Gn1,decreasing = F),collapse = ";")
      
      }
    
  }
  FinalAnn$substrate=substrates
  FinalAnn$geneNames=geneNames
  FinalAnn$PULGeneID=PULID
  write.table(FinalAnn,str_c("S:/Abilash/Marbouty_annotation_refseq_genbank 11_2021/Data/Marbouty_annotations_for_integrated_proteins/",gnmNms[i2],"_allAnnJoint.csv"),quote=F,sep = "\t",row.names = F)
  saveRDS(FinalAnn,str_c("S:/Abilash/Marbouty_annotation_refseq_genbank 11_2021/Data/Marbouty_annotations_for_integrated_proteins/",gnmNms[i2],"_allAnnJoint.rds"))
  
  eggNOGFl=read.csv(file = str_c(eggNOGPath,gnmNms[i2],"/out.emapper.annotations"),header = F,sep = "\t",comment.char = "#")
  colnames(eggNOGFl)=c("query",	"seed_ortholog",	"evalue",	"score",	"eggNOG_OGs",	"max_annot_lvl",	"COG_category",	"Description",
                       "Preferred_name",	"GOs",	"EC",	"KEGG_ko",	"KEGG_Pathway",	"KEGG_Module",	"KEGG_Reaction",	"KEGG_rclass",
                       "BRITE",	"KEGG_TC",	"CAZy",	"BiGG_Reaction",	"PFAMs")
  
  
eggNOGFl$COGID=str_replace_all(string = eggNOGFl$eggNOG_OGs,pattern = "@.*",replacement = "")
FinalAnn$eggNOG_COGID=eggNOGFl$COGID[match(FinalAnn$GeneID.x,as.character(eggNOGFl$query))]
FinalAnn$eggNOG_COG_category=eggNOGFl$COG_category[match(FinalAnn$GeneID.x,as.character(eggNOGFl$query))]  
FinalAnn$eggNOG_COG_Description=eggNOGFl$Description[match(FinalAnn$GeneID.x,as.character(eggNOGFl$query))]
FinalAnn$eggNOG_COG_PFAMs=eggNOGFl$PFAMs[match(FinalAnn$GeneID.x,as.character(eggNOGFl$query))]
FinalAnn$eggNOG_COG_BiGG_Reaction=eggNOGFl$BiGG_Reaction[match(FinalAnn$GeneID.x,as.character(eggNOGFl$query))]


write.table(FinalAnn,str_c("S:/Abilash/Marbouty_annotation_refseq_genbank 11_2021/Data/Marbouty_annotations_for_integrated_proteins/",gnmNms[i2],"_allAnnJoint_withEGGNOG.csv"),quote=F,sep = "\t",row.names = F)
saveRDS(FinalAnn,str_c("S:/Abilash/Marbouty_annotation_refseq_genbank 11_2021/Data/Marbouty_annotations_for_integrated_proteins/",gnmNms[i2],"_allAnnJoint_withEGGNOG.rds"))


}
