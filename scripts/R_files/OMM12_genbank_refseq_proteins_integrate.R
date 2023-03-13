BedFilesCreate=function(altogetherGff)
{
library("stringr")
library("seqinr")

ResultFolder="S:/Abilash/Marbouty_annotation_refseq_genbank 11_2021/Data/Marbouty_genbank_and_refseq_proteins_integrated/"

if(!(dir.exists(ResultFolder)))
{
  dir.create(ResultFolder)
  
}


BedFolder="S:/Abilash/Marbouty_annotation_refseq_genbank 11_2021/Data/Marbouty_annotations_proteins_integrated_bed/"

if(!(dir.exists(BedFolder)))
{
  dir.create(BedFolder)
  
}


proteinFiles=list.files(path = "S:/Abilash/Marbouty_annotation_refseq_genbank 11_2021/Data/Marbouty_genbank_and_refseq_proteins/",full.names = T)
unqGenomes=str_replace_all(string = proteinFiles,pattern = ".*/",replacement = "")
unqGenomes=unique(str_replace_all(string = unqGenomes,pattern = ".protein.*",replacement = ""))


for(i in 1:length(unqGenomes))
{
  GnmFiles=proteinFiles[grep(pattern = unqGenomes[i],x=proteinFiles)]
  FinFileNm=str_c(ResultFolder,unqGenomes[i],"_proteins_integrated.faa")  
  
  if(length(GnmFiles) >1)
  {
    #### Find the unique proteins from RefSeq and append it to the Genbank proteins.
    ##### How?: check with the location tag (start to stop with strand-specificity)
    ##### Change the name of the protein sequences to its locus tags
    ###### How? First assign the Annot attribute as the name and then edit based on 'locus_tag' tag
    
    
    GenbankFile=GnmFiles[grep(pattern = "genbank",x=GnmFiles)]
    RefSeqFile=GnmFiles[grep(pattern = "refseq",x=GnmFiles)]
  
    GenbankFasta=read.fasta(file = GenbankFile,seqtype = "AA")  
    RefSeqFasta=read.fasta(file = RefSeqFile,seqtype = "AA")  
    
    GenbankAttributes=unlist(lapply(GenbankFasta,function(x){attr(x = x,which = "Annot",exact = TRUE)}))
    RefSeqAttributes=unlist(lapply(RefSeqFasta,function(x){attr(x = x,which = "Annot",exact = TRUE)}))
    
    names(GenbankFasta)=GenbankAttributes
    names(RefSeqFasta)=RefSeqAttributes
    
    GenbankAttributesShort=str_replace_all(string = GenbankAttributes,pattern = ".*location=",replacement = "")
    GenbankAttributesShort=str_replace_all(string = GenbankAttributesShort,pattern = "].*",replacement = "")
    
    RefSeqAttributesShort=str_replace_all(string = RefSeqAttributes,pattern = ".*location=",replacement = "")
    RefSeqAttributesShort=str_replace_all(string = RefSeqAttributesShort,pattern = "].*",replacement = "")
    
    RefSeqUnq=RefSeqAttributesShort[match(GenbankAttributesShort,RefSeqAttributesShort)]
    RefSeqUnq=which(is.na(RefSeqUnq))
    
    names(GenbankFasta)=str_replace_all(string = names(GenbankFasta),pattern = ".*locus_tag=",replacement = "")
    names(GenbankFasta)=str_replace_all(string = names(GenbankFasta),pattern = "].*",replacement = "")
    
    names(RefSeqFasta)=str_replace_all(string = names(RefSeqFasta),pattern = ".*locus_tag=",replacement = "")
    names(RefSeqFasta)=str_replace_all(string = names(RefSeqFasta),pattern = "].*",replacement = "")
    
    
    RefSeqUnqFasta=RefSeqFasta[RefSeqUnq]
    
    NewFasta=c(GenbankFasta,RefSeqUnqFasta)
    NewFasta=lapply(NewFasta, function(x){gsub(pattern = "[[:punct:]]",replacement = "",x = x)})

    
    write.fasta(sequences = NewFasta,file.out = FinFileNm,names = names(NewFasta))
    
  } else {
    FinVersion=read.fasta(file = GnmFiles)
    FinVersionAttributes=unlist(lapply(FinVersion,function(x){attr(x = x,which = "Annot",exact = TRUE)}))
    names(FinVersion)=FinVersionAttributes
    
    
    names(FinVersion)=str_replace_all(string = names(FinVersion),pattern = ".*locus_tag=",replacement = "")
    names(FinVersion)=str_replace_all(string = names(FinVersion),pattern = "].*",replacement = "")
    
    FinVersion=lapply(FinVersion, function(x){gsub(pattern = "[[:punct:]]",replacement = "",x = x)})
    
    write.fasta(sequences = FinVersion,file.out = FinFileNm,names = names(FinVersion))
  }

}

## Create bed files:



FinFileNm=list.files(path = ResultFolder,pattern = "_proteins_integrated.faa",full.names = T)  
unqGenomes=str_replace_all(string = FinFileNm,pattern = ".*/",replacement = "")
unqGenomes=str_replace_all(string = unqGenomes,pattern = "_proteins.*",replacement = "")

for(i in 1:length(FinFileNm))
{
  fFile=read.fasta(file = FinFileNm[i],seqtype = "AA")
  nmsFile=names(fFile)
  message(str_c(i,":::",unqGenomes[i],"::: is started!!"))
  
  for(j in 1:length(nmsFile))
  {
    message(str_c(j,":::::",nmsFile[j],"::: is started!!"))
    
    tmpIdx=grep(pattern = nmsFile[j],x=altogetherGff$attributes)[1]
    Df1=data.frame(Scaffold=altogetherGff$seqid[tmpIdx],GeneID=nmsFile[j],
                   start=altogetherGff$start[tmpIdx],end=altogetherGff$end[tmpIdx],strand=altogetherGff$strand[tmpIdx])
    
    
    if(j==1)
    {
      ovDf=Df1
      
    } else {
      
      ovDf=rbind(ovDf,Df1)
    }
    
    message(str_c(j,":::::",nmsFile[j],"::: is done!!"))
    
  }
  message(str_c(i,":::",unqGenomes[i],"::: is done!!"))
  ovDf$Scaffold=str_replace_all(string = ovDf$Scaffold,pattern = ".*&",replacement = "")
  write.table(x = ovDf,file = str_c(BedFolder,unqGenomes[i],"_bedfile.txt"),quote = F,sep = "\t",row.names = F,col.names = F)
}

}

BedFilesCreate(altogetherGff = altogetherGff)