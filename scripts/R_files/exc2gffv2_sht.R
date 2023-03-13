exc2gffv2_sht=function()
{
  
  library("dplyr")
  library("stringr")
  annExcl=list.files(path ="S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/",pattern = "_new_FinMatSht.rds",full.names = T)
  gnmNames=str_replace_all(string = annExcl,pattern = ".*/",replacement = "")
  gnmNames=str_replace_all(string = gnmNames,pattern = "_new_FinMat.*",replacement = "")
  # PULFldr="/dss/dssfs02/lwp-dss-0001/u7b03/u7b03-dss-0000/ra72seh/DB/Marbouty_joint/PULs/"
  gnmFile=read.csv(file = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/GenomeID_to_Organism_accession.csv",sep = ";",header = T)
  start_cols=c("NCBIStart","dbCAN2Start","eggNOGStart","operonStart","RGIStart","PULStart")
  
  for(i in 1:length(annExcl))
  {
    
    # pulFile=str_c(PULFldr,gnmNames[i],"V10_pul.tsv")
    
    
    tmpExcl=readRDS(annExcl[i])
    # tmpExcl=(tmpExcl[!duplicated(tmpExcl$str_Crd),])
    
    # tmpExcl$Gnm=unname(apply(tmpExcl[,c("GenomeID","Organism")],1,function(x){unique(x[x!=""])[1]})) # str_c(unique(x[x!=""]),collapse = "|")
    
    
    gff_columns=c("GenomeID","Source","Type","Start","End","Score","Strand","Phase","Attributes")
    srcVec=apply(tmpExcl[,intersect(start_cols,colnames(tmpExcl))],1,function(x){
      
      nms=names(x);
      str_c(nms[which(x!="")],collapse = "&")
      
      
    })
    
    gffDf=data.frame(GenomeID=tmpExcl$Overall_GenomeID, Source=srcVec )
    
    gffDf$Type="GeneFeatures"
    gffDf[c('Start','End','Strand')]=str_split_fixed(tmpExcl$str_Crd,pattern = "_",3)
    gffDf$Score="."
    gffDf$Phase="."
    gffDf$Attributes=apply(tmpExcl[,c("str_Crd","locus_tag")],1,function(x){
      
      if(x[2]=="" || is.na(x[2]))
      {
        return(str_c("locus_tag=",x[1],";"))
          } else {
          x[2]=str_replace_all(string = x[2],pattern = "///$",replacement = "")
          x[2]=str_replace_all(string = x[2],pattern = "=///",replacement = "=")  
          return(str_c("locus_tag=",x[2],";ref_coordinatesID=",x[1],";"))
            
      }
    })
    
    gffDf=gffDf[,gff_columns]
    ####
    str_crd_vec=apply(gffDf[,c("Start","End","Strand")],1,function(x){
      str_c(x,collapse = "_")
    })
    
    refseq_coords=tryCatch(readRDS(file = str_c(MarboutyFolder,gnmNames[i],"_RefSeq_gene_coordinates.rds")),error=function(x){NA})
    genbank_coords=tryCatch(readRDS(file = str_c(MarboutyFolder,gnmNames[i],"_Genbank_gene_coordinates.rds")),error=function(x){NA})
    
    if(dim(genbank_coords)[1] > 0)
    {
      
      genbank_coords$ID=str_replace_all(string = genbank_coords$ID,pattern = " ",replacement = "")
      # saveRDS(genbank_coords,str_c(MarboutyFolder,gnmNames[i],"_Genbank_gene_coordinates.rds"))
      
    }
    
    
    if(dim(refseq_coords)[1] > 0)
    {
      refseq_coords$ID=str_replace_all(string = refseq_coords$ID,pattern = " ",replacement = "")
      # saveRDS(refseq_coords,str_c(MarboutyFolder,gnmNames[i],"_RefSeq_gene_coordinates.rds"))
      
      
      
    }
    unqLocIdx=grep(pattern = "///",x = gffDf$Attributes)
    locVec=vector()
    for(u2 in 1: length(unqLocIdx))
    {
      
      gbkIdx=which(genbank_coords$ID==str_crd_vec[unqLocIdx[u2]])
      
      if(length(gbkIdx) > 0 && is.na(locVec[u2]))
      {
        
        locVec[u2]=str_c("locus_tag=",unique(genbank_coords$Locus_tag[gbkIdx]))
      
      } else {
          
        rfsqIdx=which(refseq_coords$ID==str_crd_vec[unqLocIdx[u2]])
        
        if(length(rfsqIdx) > 0 && is.na(locVec[u2]))
        {
        
          locVec[u2]=str_c("locus_tag=",unique(refseq_coords$Locus_tag[rfsqIdx]))
          
        } else {
          locVec[u2]=str_c("locus_tag=",str_crd_vec[unqLocIdx[u2]])
          
        }
        }
      
      }
    
    gffDf$Attributes[unqLocIdx]=locVec
  gffDf$Source=str_replace_all(string = gffDf$Source,pattern = "Start",replacement = "")
    newGff=str_c("S:/Abilash/Marbouty_annotation_refseq_genbank_2022/gffFiles/For_RNAseq_count/",gnmNames[i],"_new_JointoverallAnnotationNew.gff")
    # unique(gffDf$GenomeID)
    
    write.table(file = newGff,x = gffDf,append = F,quote = F,sep = "\t",row.names = F,col.names = F)
  }
  
}

exc2gffv2_sht()


####
gff_files=list.files(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/gffFiles/For_RNAseq_count/",pattern = "_new_JointoverallAnnotationNew.gff",full.names = T)
gff_files_omm=gff_files[-grep(pattern = "escherichia|OMM",x=gff_files)]
gnmNamesv2=str_replace_all(string = gff_files_omm,pattern = ".*/",replacement = "")
gnmNamesv2=str_replace_all(string = gnmNamesv2,pattern = "_new.*",replacement = "")

for(ii in 1: length(gff_files_omm))
{
  gff_fileOMM=read.csv(file = gff_files_omm[ii],header = F,sep = "\t")
  
  gnmIdx=intersect(which(gff_fileOMM$V5==(gnmFile$Genome_size[which(as.character(gnmFile$Genome_compatible)==as.character(gnmNamesv2)[ii])])),which(gff_fileOMM$V4=="1"))
  
  if(length(gnmIdx) > 0)
  {
    
    gff_fileOMM=gff_fileOMM[-(gnmIdx),]
    
  }
  
  
  gff_fileOMM$V4=as.numeric(gff_fileOMM$V4)
  gff_fileOMM$V5=as.numeric(gff_fileOMM$V5)
  gff_fileOMM=gff_fileOMM[order(gff_fileOMM$V4,gff_fileOMM$V5),]
  
  gff_fileOMM=gff_fileOMM[with(gff_fileOMM, order(V4, V5)),  ]
  v1_col=unique(gff_fileOMM$V1)
  v2_col="genome"
  v3_col="genome"
  v4_col=1
  v5_col=gnmFile$Genome_size[which(gnmFile$GenbankID==as.character(v1_col))]
  
  v6_col="."
  v7_col="+"
  v8_col="."
  v9_col=str_c("Genome: ",gnmFile$Genome_compatible[which(gnmFile$GenbankID==as.character(v1_col))])
  gnm_col=data.frame(V1=v1_col,V2=v2_col,V3=v3_col,V4=v4_col,V5=v5_col,V6=v6_col,V7=v7_col,V8=v8_col,V9=v9_col)
  
  gff_fileOMM=rbind(gnm_col,gff_fileOMM)
  gff_fileOMM$V9=str_replace_all(string = gff_fileOMM$V9,pattern = ";.*",replacement = "")
  gff_fileOMM$V2=str_replace_all(string = gff_fileOMM$V2,pattern = "Start",replacement = "")
  
  if(ii==1)
  {
    
    overalldf=gff_fileOMM
  } else {
    
    overalldf=rbind(overalldf,gff_fileOMM)
    
    }
  
  gff_file_name=str_c("S:/Abilash/Marbouty_annotation_refseq_genbank_2022/gffFiles/For_RNAseq_count/",gnmNamesv2[ii],"_new_JointoverallAnnotationNew.gff")
  sink(file = gff_file_name)
  cat(readLines(str_c("S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Extra/gff_headers/",gnmNamesv2[ii],"_header.txt")),sep = "\n")
  sink()

  write.table(x = gff_fileOMM,file = gff_file_name,quote = F,sep = "\t",row.names = F,append = T,col.names = F)

  }


###### Overall annotation:
gff_file_name="S:/Abilash/Marbouty_annotation_refseq_genbank_2022/gffFiles/For_RNAseq_count/OMM_Overallannotations.gff"
sink(file = gff_file_name)
cat(readLines("S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Extra/gff_headers/OMM_ecoli_overall_headers"),sep = "\n")
sink()
write.table(overalldf,file = gff_file_name,quote = F,sep = "\t",row.names = F,col.names = F,append = T)

###### annotations with only RGI info:
gff_file_name="S:/Abilash/Marbouty_annotation_refseq_genbank_2022/gffFiles/For_RNAseq_count/OMM_Overallannotations_RGI.gff"
sink(file = gff_file_name)
cat(readLines("S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Extra/gff_headers/OMM_ecoli_overall_headers"),sep = "\n")
sink()
write.table(overalldf[which(overalldf$V2=="RGI"),],file = gff_file_name,quote = F,sep = "\t",row.names = F,append = T,col.names = F)

###### Overall annotation without genes from RGI alone:
gff_file_name="S:/Abilash/Marbouty_annotation_refseq_genbank_2022/gffFiles/For_RNAseq_count/OMM_Overallannotations_without_genes_from_RGI_alone.gff"
sink(file = gff_file_name)
cat(readLines("S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Extra/gff_headers/OMM_ecoli_overall_headers"),sep = "\n")
sink()
write.table(overalldf[-which(overalldf$V2=="RGI"),],file = gff_file_name,quote = F,sep = "\t",row.names = F,append = T,col.names = F)
