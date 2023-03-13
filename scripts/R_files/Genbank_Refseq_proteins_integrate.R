# install.packages("MASS")
# install.packages("ade4")
# install.packages("seqinr")
#### Note: Made folder names with Escherichia_coli_Mt1B1 to Escherichia_coli_strain_Mt1B1

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


library("seqinr")
library("stringr")
list_genbank_proteins=list.files(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Downloaded_from_genbank_and_refseq/Marbouty_Genbank/",
                                 pattern = "_translated_cds.faa",full.names = T,recursive = T)


list_refseq_proteins=list.files(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Downloaded_from_genbank_and_refseq/Marbouty_refseq/",
                                pattern = "_translated_cds.faa",full.names = T,recursive = T)

gbk_rfs_int_gff=list.files(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/gffFiles/",pattern = "_genomic.gff",full.names = T)
gbk_rfs_int_gff_genomes=str_replace_all(string = gbk_rfs_int_gff,pattern = ".*gffFiles/",replacement = "")
gbk_rfs_int_gff_genomes=str_replace_all(string = gbk_rfs_int_gff_genomes,pattern = "_genomic.gff",replacement = "")


all_proteins=c(list_genbank_proteins,list_refseq_proteins)
gnNames=str_replace_all(string = all_proteins,pattern = ".*Marbouty_Genbank//",replacement = "")
gnNames=str_replace_all(string = gnNames,pattern = ".*Marbouty_refseq//",replacement = "")
# gnNames=unlist(lapply(str_split(gnNames,pattern = "/"),function(x){unlist(x)[1]}))
gnNames=str_replace_all(string = gnNames,pattern = "\\/.*",replacement = "")

unqGnNames=unique(gnNames)





for(pr1 in 1: length(unqGnNames))
{
  
  tmpIdx=which(gnNames==unqGnNames[pr1])  
  finPrName=str_c("S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_genbank_and_refseq_proteins_integrated/",unqGnNames[pr1],"_proteins_integrated.faa")
  
  if(length(tmpIdx)==1)
  {
    rr=read.fasta(file = all_proteins[tmpIdx],seqtype = "AA",whole.header = T)
    
    nn1=names(rr)
    
    nn1=str_replace_all(string = nn1,pattern = ".*locus_tag=",replacement = "")
    nn1=str_replace_all(string = nn1,pattern = "].*",replacement = "")
    names(rr)=nn1
    
    write.fasta(sequences = rr,names = names(rr),
               file.out = finPrName)
    
  
    } else {
    
    fl1=read.fasta(file=all_proteins[tmpIdx[1]],seqtype = "AA",whole.header = T)
    fl2=read.fasta(file=all_proteins[tmpIdx[2]],seqtype = "AA",whole.header = T)
    nn1=names(fl1)
    nn2=names(fl2)
    
    nn1=str_replace_all(string = nn1,pattern = ".*locus_tag=",replacement = "")
    nn1=str_replace_all(string = nn1,pattern = "].*",replacement = "")
    
    nn2=str_replace_all(string = nn2,pattern = ".*locus_tag=",replacement = "")
    nn2=str_replace_all(string = nn2,pattern = "].*",replacement = "")
    
    names(fl1)=nn1
    names(fl2)=nn2
    
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
    
    finPr=finPr[match(unique(names(finPr)),names(finPr))]
    write.fasta(sequences = finPr,names = names(finPr),
             file.out = finPrName)
    
    
    }
  

  
 
}


###### Checked whether there is any change in # of proteins happen before & after september: No changes, hence continue with the analysis
# [1] "####################################################################"
# [1] "New method proteins: 3886"
# [1] "Old method protiens: 3886"
# [1] "Common #of proteins:3886"
# [1] "####################################################################"
# [1] "####################################################################"
# [1] "New method proteins: 2252"
# [1] "Old method protiens: 2252"
# [1] "Common #of proteins:2252"
# [1] "####################################################################"
# [1] "####################################################################"
# [1] "New method proteins: 3933"
# [1] "Old method protiens: 3933"
# [1] "Common #of proteins:3933"
# [1] "####################################################################"
# [1] "####################################################################"
# [1] "New method proteins: 1608"
# [1] "Old method protiens: 1608"
# [1] "Common #of proteins:1608"
# [1] "####################################################################"
# [1] "####################################################################"
# [1] "New method proteins: 4728"
# [1] "Old method protiens: 4728"
# [1] "Common #of proteins:4728"
# [1] "####################################################################"
# [1] "####################################################################"
# [1] "New method proteins: 4263"
# [1] "Old method protiens: 4263"
# [1] "Common #of proteins:4263"
# [1] "####################################################################"
# [1] "####################################################################"
# [1] "New method proteins: 2790"
# [1] "Old method protiens: 2790"
# [1] "Common #of proteins:2790"
# [1] "####################################################################"
# [1] "####################################################################"
# [1] "New method proteins: 5226"
# [1] "Old method protiens: 5226"
# [1] "Common #of proteins:5226"
# [1] "####################################################################"
# [1] "####################################################################"
# [1] "New method proteins: 3689"
# [1] "Old method protiens: 3689"
# [1] "Common #of proteins:3689"
# [1] "####################################################################"
# [1] "####################################################################"
# [1] "New method proteins: 1890"
# [1] "Old method protiens: 1890"
# [1] "Common #of proteins:1890"
# [1] "####################################################################"
# [1] "####################################################################"
# [1] "New method proteins: 2675"
# [1] "Old method protiens: 2675"
# [1] "Common #of proteins:2675"
# [1] "####################################################################"
# [1] "####################################################################"
# [1] "New method proteins: 2527"
# [1] "Old method protiens: 2527"
# [1] "Common #of proteins:2527"
# [1] "####################################################################"