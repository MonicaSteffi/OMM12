exc2gffv2=function()
{
  
  library("dplyr")
  library("stringr")
  annExcl=list.files(path ="S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/",pattern = "_FinMatSht.rds",full.names = T)
  gnmNames=str_replace_all(string = annExcl,pattern = ".*/",replacement = "")
  gnmNames=str_replace_all(string = gnmNames,pattern = "_FinDf.*",replacement = "")
  # PULFldr="/dss/dssfs02/lwp-dss-0001/u7b03/u7b03-dss-0000/ra72seh/DB/Marbouty_joint/PULs/"
  
  
  for(i in 1:length(annExcl))
  {
    
    # pulFile=str_c(PULFldr,gnmNames[i],"V10_pul.tsv")
    
    
    tmpExcl=readRDS(annExcl[i])
    tmpExcl$GeneStart=as.numeric(tmpExcl$GeneStart)
    tmpExcl=tmpExcl %>% group_by(Overall_GenomeID,GeneStart) %>% arrange(Overall_GenomeID,GeneStart)
    tmpExcl=(tmpExcl[!duplicated(tmpExcl$str_Crd),])
    tmpExcl$GenomeID=as.character(tmpExcl$GenomeID)
    
    # tmpExcl$Gnm=unname(apply(tmpExcl[,c("GenomeID","Organism")],1,function(x){unique(x[x!=""])[1]})) # str_c(unique(x[x!=""]),collapse = "|")
    
    koidx=which(!is.na(tmpExcl$KO))
    gffidx=which(!is.na(tmpExcl$ProteinID))
    cgcidx=which(!is.na(tmpExcl$CGCID))
    operonidx=which(!is.na(tmpExcl$Operon))
    
    tmpExcl[is.na(tmpExcl)]=""
    ######### KOFAM
    koMat=tmpExcl[koidx,]  
    koMatAttr=as.data.frame(koMat[,c("KO_definition","GeneID.x","F_score","KO","threshold","score","evalue")])
    
    for(i3 in 1:ncol(koMatAttr))
    {
      koMatAttr[,i3]=str_c(colnames(koMatAttr)[i3],":",koMatAttr[,(i3)])
      
    }
    koMatAttrVec=apply(koMatAttr,1,function(x){str_c(x,collapse = ";")})
    
    koMatgff=data.frame(Organism=koMat$Gnm  )
    koMatgff$Source="KofamScan"
    koMatgff$FeatureType="KO"
    # koMatgff$FeatureType=koMat$KO_definition
    koMatgff$Start=unlist(apply(koMat[,c("Start.x","Start.y")],1,function(x){str_c(unique(x[x!=""]),collapse = "|")}))
    koMatgff$Stop=unlist(apply(koMat[,c("Stop","End")],1,function(x){str_c(unique(x[x!=""]),collapse = "|")}))
    koMatgff$Score="."
    koMatgff$Strand=unlist(apply(koMat[,c("Strand","Strand.x","Strand.y")],1,function(x){str_c(unique(x[x!=""]),collapse = "|")}))
    koMatgff$Phase="."
    koMatgff$Attributes=koMatAttrVec
    
    ######### Genbank
    gffMat=as.data.frame(tmpExcl[gffidx,])
    gffCol=c("Product","str_Crd","GeneID.x","ProteinID","Name", "RptFamily" ,"gbkey" ,"Dbxref","locus_tag" , "gene"  , "Note" ,  "regulatory_class", "bound_moiety")
    gffMatAttr=as.data.frame(gffMat[,intersect(gffCol,colnames(gffMat))])
    for(i3 in 1:ncol(gffMatAttr))
    {
        gffMatAttr[,i3]=str_c(colnames(gffMatAttr)[i3],":",gffMatAttr[,i3])  
      
      
      
    }
    gffMatAttrVec=apply(gffMatAttr,1,function(x){str_c(x,collapse = ";")})
    
    gffMatgff=data.frame(Organism=gffMat$Gnm  )
    gffMatgff$Source="Genbank"
    gffMatgff$FeatureType="GFF"
    # gffMatgff$FeatureType=gffMat$Product
    gffMatgff$Start=unlist(apply(gffMat[,c("Start.x","Start.y")],1,function(x){str_c(unique(x[x!=""]),collapse = "|")}))
    gffMatgff$Stop=unlist(apply(gffMat[,c("Stop","End")],1,function(x){str_c(unique(x[x!=""]),collapse = "|")}))
    gffMatgff$Score="."
    gffMatgff$Strand=unlist(apply(gffMat[,c("Strand","Strand.x","Strand.y")],1,function(x){str_c(unique(x[x!=""]),collapse = "|")}))
    gffMatgff$Phase="."
    gffMatgff$Attributes=gffMatAttrVec
    gffMatgff$Attributes=str_replace_all(string = gffMatgff$Attributes,pattern = "GeneID.x:",replacement = "GeneID:")
    
    ######### dbCAN
    cgcMat=as.data.frame(tmpExcl[cgcidx,]  )
    
    cgcMatAttr=cgcMat[,intersect(c("Product_Name","str_Crd","DB","CGCID","GenomeID","GeneID.x","Description","Product_CGC","substrate","geneNames"),colnames(gffMat))]
    for(i3 in 1:ncol(cgcMatAttr))
    {
      
      if(colnames(cgcMatAttr)[i3]=="Product_Name")
      {
        cgcMatAttr[,i3]=str_c("Product_Name",":",cgcMatAttr[,i3])
      } else {
        cgcMatAttr[,i3]=str_c(colnames(cgcMatAttr)[i3],":",cgcMatAttr[,i3])  
      }
      
      
    }
    cgcMatAttrVec=apply(cgcMatAttr,1,function(x){str_c(x,collapse = ";")})
    
    cgcMatgff=data.frame(Organism=cgcMat$Gnm )
    cgcMatgff$Source="dbCAN_V10"
    cgcMatgff$FeatureType="PUL"
    # cgcMatgff$FeatureType=cgcMat$Product_Name
    cgcMatgff$Start=unlist(apply(cgcMat[,c("Start.x","Start.y")],1,function(x){str_c(unique(x[x!=""]),collapse = "|")}))
    cgcMatgff$Stop=unlist(apply(cgcMat[,c("Stop","End")],1,function(x){str_c(unique(x[x!=""]),collapse = "|")}))
    cgcMatgff$Score="."
    cgcMatgff$Strand=unlist(apply(cgcMat[,c("Strand","Strand.x","Strand.y")],1,function(x){str_c(unique(x[x!=""]),collapse = "|")}))
    cgcMatgff$Phase="."
    cgcMatgff$Attributes=cgcMatAttrVec
    
    
    
    ######### operonMapper
    operonMat=as.data.frame(tmpExcl[operonidx,]  )
    operonMatAttr=operonMat[,c("Operon","str_Crd", "GeneID.y" , "Type" , "COGgene" , "PosLeft" , "postRight" ,"Strand.y" , "Function")]
    for(i3 in 1:ncol(operonMatAttr))
    {
      operonMatAttr[,i3]=str_c(colnames(operonMatAttr)[i3],":",operonMatAttr[,i3])
      
    }
    operonMatAttrVec=apply(operonMatAttr,1,function(x){str_c(x,collapse = ";")})
    
    operonMatgff=data.frame(Organism=operonMat$Gnm  )
    operonMatgff$Source="OperonMapper"
    operonMatgff$FeatureType="Operon"
    # operonMatgff$FeatureType=operonMat$Function
    operonMatgff$Start=unlist(apply(operonMat[,c("Start.x","Start.y")],1,function(x){str_c(unique(x[x!=""]),collapse = "|")}))
    operonMatgff$Stop=unlist(apply(operonMat[,c("Stop","End")],1,function(x){str_c(unique(x[x!=""]),collapse = "|")}))
    operonMatgff$Score="."
    operonMatgff$Strand=unlist(apply(operonMat[,c("Strand","Strand.x","Strand.y")],1,function(x){str_c(unique(x[x!=""]),collapse = "|")}))
    operonMatgff$Phase="."
    operonMatgff$Attributes=operonMatAttrVec
    
    overallMatgff=rbind(koMatgff,gffMatgff,cgcMatgff,operonMatgff)
    
    if(file.size(pulFile)>0)
    {
      pulFl=read.csv(file = pulFile,header = T,sep = "\t") 
      pulcl=c("pulid","protein_id","dist","protein_name","sus","hmm","active","pattern")
      pulMat=pulFl 
      pulMatAttr=pulMat[,intersect(pulcl,colnames(pulMat))]
      pulFl$genome=str_replace_all(string = pulFl$genome,pattern = "V10",replacement = "")
      pulMatID=apply(pulFl[,c("genome","start","end")],1,function(x){str_c(x,collapse = "_")})
      pulMatAttr=cbind(pulMatAttr,pulMatID)
      colnames(pulMatAttr)[ncol(pulMatAttr)]="str_Crd"
      pulMatAttr[is.na(pulMatAttr)]=""
      for(i3 in 1:ncol(pulMatAttr))
      {
        pulMatAttr[,i3]=str_c(colnames(pulMatAttr)[i3],":",as.character(pulMatAttr[,i3]))
        
      }
      pulMatAttrVec=apply(pulMatAttr,1,function(x){str_c(as.character(x),collapse = ";")})
      
      pulMatgff=data.frame(Organism=pulMat$contig  )
      pulMatgff$Source="PULpy"
      pulMatgff$FeatureType="susC_susD"
      # pulMatgff$FeatureType=pulMat$sus
      pulMatgff$Start=pulMat$start
      pulMatgff$Stop=pulMat$end
      pulMatgff$Score="."
      pulMatgff$Strand=pulMat$strand
      pulMatgff$Phase="."
      pulMatgff$Attributes=pulMatAttrVec
      
      overallMatgff=rbind(overallMatgff,pulMatgff)
      
    }
    overallMatgff$Attributes=str_replace_all(string = overallMatgff$Attributes,pattern = "GeneID.y:",replacement = "GeneID:")
    overallMatgff$Start=as.numeric(overallMatgff$Start)
    overallMatgff$Stop=as.numeric(overallMatgff$Stop)
    overallMatgff=overallMatgff %>% arrange(Start)
    saveRDS(overallMatgff,file = str_c("/dss/dssfs02/lwp-dss-0001/u7b03/u7b03-dss-0000/ra72seh/DB/Marbouty_joint/",gnmNames[i],"_overallAnnotation.gff.rds"))
    write.table(overallMatgff,file = str_c("/dss/dssfs02/lwp-dss-0001/u7b03/u7b03-dss-0000/ra72seh/DB/Marbouty_joint/",gnmNames[i],"_overallAnnotation.gff"),quote = F,sep = "\t",row.names = F,col.names = F)
    saveRDS(overallMatgff,file = str_c("/dss/dssfs02/lwp-dss-0001/u7b03/u7b03-dss-0000/ra72seh/DB/Marbouty_joint/",gnmNames[i],"_overallAnnotation.gff.rds"))
    
    ####
    
    defGff=str_c("/dss/dssfs02/lwp-dss-0001/u7b03/u7b03-dss-0000/ra72seh/DB/Marbouty_joint/gffFiles/",gnmNames[i],"_genomic.gff")
    
    
    con <- file(defGff, "r", blocking = FALSE)
    singleString <- readLines(con) # empty
    close(con)
    
    newGff=str_c("/dss/dssfs02/lwp-dss-0001/u7b03/u7b03-dss-0000/ra72seh/DB/Marbouty_joint/",gnmNames[i],"_JointoverallAnnotationNew.gff")
    
    txt <- str_c(singleString[1:9],collapse = "\n")
    conFile=file(newGff)
    writeLines(txt, conFile)
    close(conFile)
    
    write.table(file = newGff,x = overallMatgff,append = T,quote = F,sep = "\t",row.names = F,col.names = F)
  }
  
}

exc2gffv2()