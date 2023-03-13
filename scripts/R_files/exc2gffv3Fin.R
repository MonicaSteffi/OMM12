exc2gffv3Fin=function(AnnFile,genomeID,type1,start1,end1,strand1,attrVec,annType,newGff_folder,gnmNames)
{
  
  library("dplyr")
  library("stringr")
  annExcl=AnnFile
  
  gff_columns=c("GenomeID","Source","Type","Start","End","Score","Strand","Phase","Attributes")
  gnmFile=read.csv(file = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/GenomeID_to_Organism_accession.csv",sep = ";",header = T)
  
  if(annType =="RGI")
  {
    tmpExcl=read.csv(file = annExcl,header = T,sep = "\t")
    tmpExcl[c('Organism','GNm')]    =str_split_fixed(string = tmpExcl$Contig,pattern = "_",n = 2)
    tmpExcl$Type="Antibiotic_resistance_gene"
  } else {
    if(annType == "PULs")
    {
      tmpExcl=read.csv(file = annExcl,header = T,sep = "\t")
      
    } else {
      tmpExcl=readRDS(annExcl)
      
    }
    
  }
  
  if(annType=="eggNOGDB")
  {
    tmpExcl[c('Start','End','Strand')]=str_split_fixed(string = tmpExcl$str_Crd,pattern = "_",n = 3)
    tmpExcl$COG=str_replace_all(string = tmpExcl$eggNOG_OGs,pattern = "@.*",replacement = "")
    tmpExcl$GeneType="GeneFeatures"
    tmpExcl[c('GenomeID','GNm')]=str_split_fixed(string = tmpExcl$eggNOG_query,pattern = "_",n = 2)
  }
  
  if(annType=="NovelFamDB")
  {
    tmpExcl[c('Start','End','Strand')]=str_split_fixed(string = tmpExcl$str_Crd,pattern = "_",n = 3)
    tmpExcl$GeneType="GeneFeatures"
    tmpExcl[c('GenomeID','GNm')]=str_split_fixed(string = tmpExcl$query,pattern = "_",n = 2)
  }
  
  
  if(annType== "kofam")
  {
    gff3=readRDS(str_c("S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/",gnmNames,"_gff_results.txt.rds"))
    
    tmpExcl=as.data.frame(tmpExcl)
    tmpExcl$str_Crd=gff3$str_Crd[match(tmpExcl$GeneID,as.character(gff3$locus_tag))]
    tmpExcl$GeneType="GeneFeatures"
    na_genes=which(is.na(tmpExcl$str_Crd))
    
    if(length(na_genes) > 0)
    {
      
      tmpExcl$str_Crd[na_genes]=gff3$str_Crd[match(tmpExcl$GeneID[na_genes],as.character(gff3$locus_tag_refseq))]
      
    }
    
    na_genes=which(is.na(tmpExcl$str_Crd))
    
    if(length(na_genes) > 0)
    {
      
      tmpExcl$str_Crd[na_genes]=gff3$str_Crd[match(tmpExcl$GeneID[na_genes],as.character(gff3$locus_tag_genbank))]
      
    }
    
    tmpExcl[c('Start','End','Strand')]=str_split_fixed(string = tmpExcl$str_Crd,pattern = "_",n = 3)
    tmpExcl$Organism=gff3$Organism[match(tmpExcl$str_Crd,as.character(gff3$str_Crd))]
    
    
  }
  
  
  if(annType == "PULs")
  {
    
    tmpExcl[,type1]="PULs"
    
  }
  
  tmpExcl$DBType=annType
  
  gffDf=data.frame(Start=tmpExcl[,start1],End=tmpExcl[,end1],Strand=tmpExcl[,strand1],Type=(tmpExcl[,type1]),GenomeID=tmpExcl[,genomeID])
  
  gffDf$Source=annType
  gffDf$Score="."
  gffDf$Phase="."
  
  attrDf=as.data.frame(tmpExcl[,attrVec])
  colnames(attrDf)[1]="Name"
  attrPrfx=str_c(unlist(lapply(1:ncol(attrDf),function(x){colnames(attrDf)[x]})),"=")
  attrList=lapply(1:length(attrPrfx),function(x){
    
    str_c(attrPrfx[x],(attrDf[,x]))
  })
  names(attrList)=colnames(attrDf)
  attrDfFin=t(do.call("rbind",attrList))
  
  gffDf$Attributes=apply(attrDfFin,1,function(x){
    x=x[!is.na(x)];
    x=x[x!=""];
    str_c(x,collapse = ";")
  })
  
  gffDf=gffDf[,gff_columns]
  
  
  ####
  
  header_File=readLines(con = str_c("S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Extra/gff_headers/",gnmNames,"_header.txt"))
  
  newGff=str_c(newGff_folder,gnmNames,"_",annType,"_genious_annotationv2.gff")
  # unique(gffDf$GenomeID)
  # sink(newGff)
  # print(str_c("S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Extra/gff_headers/",gnmNames,"_header.txt"))
  # sink()
  # 
  
  sink(newGff)
  cat(header_File,sep = "\n") 
  sink()
  write.table(file = newGff,x = gffDf,quote = F,sep = "\t",row.names = F,col.names = F,append = T)
  saveRDS(object = gffDf,file   = str_c(newGff,".rds"))
  
  
}

