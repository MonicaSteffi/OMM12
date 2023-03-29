Genbank_RefSeq_Locus_tags_coordinates=function(x)
{
  library("ape")
  library("stringr")
  MarboutyFolder="S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/"
  list_files_genbank=list.files(path="S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_genbank_and_refseq_gffFiles/Marbouty_Genbank/",full.names = T)
  gnkNm=str_replace_all(string = list_files_genbank,pattern = ".*Marbouty_Genbank/",replacement = "")
  gnkNm=str_replace_all(string = gnkNm,pattern = "_genomic.gff",replacement = "")
  
  setClass(Class = "gffAnnotate",representation("FinalGff"="data.frame","FinalGff2"="data.frame",
                                                "FinalGff3"="data.frame","FinalGff4"="data.frame",
                                                overlapGn="list",sbstrr="list"))
  
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
      
    
        
      
    } else {
      
      overallgff=read.gff(file = tmpLnk[1])
    }
    
    
  }
  
  
}




############
reqgffFile=gbk_rfs_int_gff[which(gbk_rfs_int_gff_genomes==unqGnNames[pr1])]   
gffFile=read.csv(file = reqgffFile,header = F,sep = "\t",comment.char = "#")$V9
gffV9=str_split(string = gffFile,pattern = ";")


locus_refseq=str_trim(str_replace_all(string = gffFile,pattern = ".*;locus_tag_RefSeq=",replacement = ""),side = "both")
locus_refseq=str_trim(str_replace_all(string = locus_refseq,pattern = ";.*",replacement = ""),side = "both")

locus_genbank=str_trim(str_replace_all(string = gffFile,pattern = ".*;locus_tag_GenBank=",replacement = ""),side = "both")
locus_genbank=str_trim(str_replace_all(string = locus_genbank,pattern = ";.*",replacement = ""),side = "both")


st_end_strnd=str_trim(str_replace_all(string = gffFile,pattern = ".*;Start_End_strand=",replacement = ""),side = "both")
st_end_strnd=str_trim(str_replace_all(string = st_end_strnd,pattern = ";.*",replacement = ""),side = "both")
names(st_end_strnd)=locus_refseq

GeneIDs=data.frame(refseq=locus_refseq,genbank=locus_genbank)
GeneIDs=GeneIDs[!duplicated(GeneIDs),]
GeneIDs$st_end=st_end_strnd[match(GeneIDs$refseq,names(st_end_strnd))]


newpr=fl1
# leftOverGeneIDs=GeneIDs[setdiff(1:nrow(GeneIDs) ,match(names(fl1),GeneIDs$genbank  )),]
ll=unlist(lapply(names(fl1),function(x){which(GeneIDs$genbank==x)}))
mm=unlist(lapply(GeneIDs$refseq[ll],function(x){which(GeneIDs$refseq==x)}))

leftOverGeneIDs=GeneIDs[-unique(c(ll,mm)),]


leftOverGeneIDs$refseq=str_replace_all(string = leftOverGeneIDs$refseq,pattern = ".*=",replacement = "")
leftOverGeneIDs$refseq=str_replace_all(string = leftOverGeneIDs$refseq,pattern = ".*-",replacement = "")

mtch=match(leftOverGeneIDs$refseq,names(fl2))
mtch=mtch[!is.na(mtch)]
leftOverfl2=fl2[mtch]
finPr=c(fl1,leftOverfl2)
