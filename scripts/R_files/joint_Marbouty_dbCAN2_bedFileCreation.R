

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