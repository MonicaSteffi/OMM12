library("seqinr")

operonProtein=list.files(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/OperonMapper_new/",pattern = "predicted_protein_sequences",full.names = T,recursive = T)
GenomeNames=str_replace_all(string = operonProtein,pattern = ".*//",replacement = "")
GenomeNames=str_replace_all(string = GenomeNames,pattern = "/.*",replacement = "")
totLen=0
for(i in 1: length(operonProtein))
{

  protFile=read.fasta(file = operonProtein[i],seqtype = "AA",whole.header = T)  
  names(protFile)=str_c(GenomeNames[i],";",names(protFile))
  
  totLen=totLen+length(protFile)
  print(str_c(GenomeNames[i]," : ", length(protFile)))
  if(i==1)
  {
    overallProt=protFile
    
    
  } else {
    
    
    overallProt=c(overallProt,protFile)
  }
  write.fasta(sequences = protFile,names = names(protFile),file.out = str_c("S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/OperonMapper_new/",GenomeNames[i],"_operonMapper_proteins.fasta"))
}

write.fasta(sequences = overallProt,names = names(overallProt),file.out = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/OperonMapper_new/OperonMapper_overallproteins.fasta")
