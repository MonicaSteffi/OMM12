Genbank_RefSeq_proteins_integrate_based_on_coordinates=function()
{
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
gn_coord_files=list.files(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/",pattern = "gene_coordinates.rds",full.names = T)
gn_coord_files_genomes=str_replace_all(string = gn_coord_files,pattern = ".*/",replacement = "")
gn_coord_files_genomes=str_replace_all(string = gn_coord_files_genomes,pattern = "_gene_co.*",replacement = "")

gn_coordList=lapply(gn_coord_files,readRDS)
names(gn_coordList)=gn_coord_files_genomes


list_genbank_proteins=list.files(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Downloaded_from_genbank_and_refseq/Marbouty_Genbank/",
                                 pattern = "_translated_cds.faa",full.names = T,recursive = T)


list_refseq_proteins=list.files(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Downloaded_from_genbank_and_refseq/Marbouty_refseq/",
                                pattern = "_translated_cds.faa",full.names = T,recursive = T)

refseq_proteinsList=lapply(list_refseq_proteins,function(x){
  
  read.fasta(file = x,seqtype = "AA",whole.header = T,as.string = T,set.attributes = F)
})
rf_prt_nms=str_replace_all(string = list_refseq_proteins,pattern = ".*_refseq/",replacement = "")
rf_prt_nms=str_replace_all(string = rf_prt_nms,pattern = "/GC.*",replacement = "")
rf_prt_nms=str_replace_all(string = rf_prt_nms,pattern = "/",replacement = "")
names(refseq_proteinsList)=rf_prt_nms

genbank_proteinsList=lapply(list_genbank_proteins,function(x){
  
  read.fasta(file = x,seqtype = "AA",whole.header = T,as.string = T,set.attributes = F)
})
gbk_prt_nms=str_replace_all(string = list_genbank_proteins,pattern = ".*_Genbank/",replacement = "")
gbk_prt_nms=str_replace_all(string = gbk_prt_nms,pattern = "/GC.*",replacement = "")
gbk_prt_nms=str_replace_all(string = gbk_prt_nms,pattern = "/",replacement = "")
names(genbank_proteinsList)=gbk_prt_nms

gnmID2Org=read.csv(file = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/GenomeID_to_Organism.csv",header = T,sep = ";")
gnmID2Org$GenbankID_edited=str_replace_all(string = gnmID2Org$GenbankID,pattern = "\\..*",replacement = "")

ID2Name2=gnmID2Org
ID2Name2$OrganismFin=ID2Name2$Organism
ID2Name2$OrganismFin=str_replace_all(string = ID2Name2$OrganismFin,pattern = "Akkermansia",replacement = "akkermansia")
ID2Name2$OrganismFin=str_replace_all(string = ID2Name2$OrganismFin,pattern = "Escherichia",replacement = "escherichia")
ID2Name2$OrganismFin=str_replace_all(string = ID2Name2$OrganismFin,pattern = "muribaculum_intestinale_YL27",replacement = "muribaculum_intestinales_YL27")

ID2Name2$GenomeFin=ID2Name2$Genome
ID2Name2$GenomeFin=str_replace_all(string = ID2Name2$GenomeFin,pattern = "Akkermansia",replacement = "akkermansia")
ID2Name2$GenomeFin=str_replace_all(string = ID2Name2$GenomeFin,pattern = "Escherichia",replacement = "escherichia")
ID2Name2$GenomeFin=str_replace_all(string = ID2Name2$GenomeFin,pattern = "muribaculum_intestinale_YL27",replacement = "muribaculum_intestinales_YL27")
unqGenomes=unique(ID2Name2$GenomeFin)
overallGenomes=list()
for(ii in 1: length(unqGenomes))
{
  
  Genbank_coords=gn_coordList[[str_c(unqGenomes[ii],"_Genbank")]]
  RefSeq_coords=gn_coordList[[str_c(unqGenomes[ii],"_RefSeq")]]
  Genbank_coords=Genbank_coords[!duplicated(Genbank_coords),]
  RefSeq_coords=RefSeq_coords[!duplicated(RefSeq_coords),]
  
  
  if(dim(RefSeq_coords)[1] > 0)
  {
    Genbank_coords$seqid=str_replace_all(string = Genbank_coords$seqid,pattern = "NZ_",replacement = "")
    RefSeq_coords$seqid=str_replace_all(string = RefSeq_coords$seqid,pattern = "NZ_",replacement = "")
    
    Genbank_coords$ID=apply(Genbank_coords[,c("seqid","ID")],1,function(x){
      str_c(x,collapse = "_")
    })
    
    RefSeq_coords$ID=apply(RefSeq_coords[,c("seqid","ID")],1,function(x){
      str_c(x,collapse = "_")
    })
    
    joint_coords=merge(Genbank_coords,RefSeq_coords,by="ID",all=TRUE)
    
    genbank_proteins=genbank_proteinsList[[unqGenomes[ii]]]
    
    if(length(genbank_proteins) > 0)
    {
      names(genbank_proteins)=str_replace_all(string = names(genbank_proteins),pattern = ".*\\[locus_tag=",replacement = "")
      names(genbank_proteins)=str_replace_all(string = names(genbank_proteins),pattern = "\\].*",replacement = "")
      
    }
    
    refseq_proteins=refseq_proteinsList[[unqGenomes[ii]]]
    if(length(refseq_proteins) > 0)
    {
      names(refseq_proteins)=str_replace_all(string = names(refseq_proteins),pattern = ".*\\[locus_tag=",replacement = "")
      names(refseq_proteins)=str_replace_all(string = names(refseq_proteins),pattern = "\\].*",replacement = "")
      
    }
    
    
    genbankIDs_in_genbank_proteinFile=intersect(joint_coords$Locus_tag.x,names(genbank_proteins))
    genbankIDs_NOTin_genbank_proteinFile=setdiff(joint_coords$Locus_tag.x,names(genbank_proteins))
    
    genbankIDS_in_refseq_proteinFile=intersect(genbankIDs_NOTin_genbank_proteinFile,joint_coords$Locus_tag.y)
    refseqIDs_unique=setdiff(unique(joint_coords$Locus_tag.y),unique(joint_coords$Locus_tag.x))
    
    genbankIDsFin=genbankIDs_in_genbank_proteinFile
    refseqIDsFin=unique(c(genbankIDS_in_refseq_proteinFile,refseqIDs_unique))
    
    tmpGenbank_proteins=genbank_proteins[intersect(genbankIDsFin,names(genbank_proteins))]
    tmpRefSeq_proteins=refseq_proteins[intersect(refseqIDsFin,names(refseq_proteins))]       
    FinProteins=c(tmpGenbank_proteins,tmpRefSeq_proteins)
    
  } else {
    
    genbank_proteins=genbank_proteinsList[[unqGenomes[ii]]]
    
    if(length(genbank_proteins) > 0)
    {
      names(genbank_proteins)=str_replace_all(string = names(genbank_proteins),pattern = ".*\\[locus_tag=",replacement = "")
      names(genbank_proteins)=str_replace_all(string = names(genbank_proteins),pattern = "\\].*",replacement = "")
      
    }
    
    FinProteins=genbank_proteins[intersect(Genbank_coords$Locus_tag,names(genbank_proteins))]
    
    
  }
  finPrfaa=str_c("S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_genbank_and_refseq_proteins_integrated/",unqGenomes[ii],"_proteins_integrated.faa")
  finPrRDS=str_c("S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_genbank_and_refseq_proteins_integrated/",unqGenomes[ii],"_proteins_integrated.rds")
    
    saveRDS(FinProteins,str_c(finPrRDS))
    write.fasta(sequences = FinProteins,names = names(FinProteins),
                file.out = finPrfaa)
    
    
   
}

rdsFile=list.files(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_genbank_and_refseq_proteins_integrated/",pattern = "_proteins_integrated.rds",full.names = T)
rdsGnm=str_replace_all(string = rdsFile,pattern = ".*/",replacement = "")
rdsGnm=str_replace_all(string = rdsGnm,pattern = "_protein.*",replacement = "")
overallGenomes=list()
for(ii in 1: length(rdsFile))
{
  readFile=readRDS(rdsFile[ii])
  names(readFile)=str_c(rdsGnm[ii],"_",names(readFile))
  if(ii==1)
  {
    
    overallGenomes=readFile
  } else {
    overallGenomes=c(overallGenomes,readFile)
    
  }
  
  
}
write.fasta(sequences = overallGenomes,names = names(overallGenomes),file.out = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/OMM12_Ecoli_proteins.faa")


##### Check if all the unique proteins are present:
####### Our aim, was to first take all the unique locus tags from genbank and refseq (in the "unq_coords" variable), then find common proteins between
######### "unq_coords" variable and the "tmpNames" variable, which is a vector of all the unique protein names.
length_of_integrated_proteins_vec=rep(0,length(unqGenomes))
length_of_Genbank_refseq_proteins_tr_cds_vec=rep(0,length(unqGenomes))

for(ii in 1: length(unqGenomes))
{
   
  finPrRDS=str_c("S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_genbank_and_refseq_proteins_integrated/",unqGenomes[ii],"_proteins_integrated.rds")
  length_of_integrated_proteins=length(readRDS(finPrRDS))
  length_of_integrated_proteins_vec[ii]=length_of_integrated_proteins

  Genbank_coords=gn_coordList[[str_c(unqGenomes[ii],"_Genbank")]]
  RefSeq_coords=gn_coordList[[str_c(unqGenomes[ii],"_RefSeq")]]
  Genbank_coords=Genbank_coords[!duplicated(Genbank_coords),]
  RefSeq_coords=RefSeq_coords[!duplicated(RefSeq_coords),]
  unq_gbk_coords=unique(Genbank_coords$Locus_tag)
  unq_coords=unq_gbk_coords
  unq_coords=c(unq_gbk_coords,setdiff(RefSeq_coords$Locus_tag,unq_gbk_coords))
  
  tmpRefseq=tryCatch(read.fasta(list_refseq_proteins[grep(pattern = unqGenomes[ii],x=list_refseq_proteins)],seqtype = "AA",as.string = T,whole.header = T),error=function(x){""})
  tmpGenbank=tryCatch(read.fasta(list_genbank_proteins[grep(pattern = unqGenomes[ii],x=list_genbank_proteins)],seqtype = "AA",as.string = T,whole.header = T),error=function(x){""})
  tmpNames=c(names(tmpRefseq),names(tmpGenbank))
  tmpNames=str_replace_all(string = tmpNames,pattern = ".*\\[locus_tag=",replacement = "")
  tmpNames=unique(str_replace_all(string = tmpNames,pattern = "\\].*",replacement = ""))
  length_of_Genbank_refseq_proteins_tr_cds=length(intersect(unq_coords,tmpNames))
  length_of_Genbank_refseq_proteins_tr_cds_vec[ii]=length_of_Genbank_refseq_proteins_tr_cds
 
  
  }
 
Nr_proteins_check=data.frame(Genome=unqGenomes,Proteins_in_integrated_file=length_of_integrated_proteins_vec,Proteins_common_bw_int_n_rfsq_gbk=length_of_Genbank_refseq_proteins_tr_cds_vec)  
  
}


#########################################################################
# > Nr_proteins_check
# Genome Proteins_in_integrated_file Proteins_common_bw_int_n_rfsq_gbk
# 1             acutalibacter_muris_KB18                        3900                              3900
# 2         akkermansia_muciniphila_YL44                        2250                              2250
# 3           bacteroides_caecimuris_I48                        3924                              3924
# 4         bifidobacterium_animalis_YL2                        1606                              1606
# 5               blautia_coccoides_YL58                        4722                              4722
# 6             clostridium_innocuum_I46                        4263                              4263
# 7  enterocloster_clostridioformis_YL32                        6923                              6923
# 8            enterococcus_faecalis_KB1                        2791                              2791
# 9          flavonifractor_plautii_YL31                        3687                              3687
# 10     limosilactobacillus_reuteri_I49                        1890                              1890
# 11       muribaculum_intestinales_YL27                        2674                              2674
# 12              turicimonas_muris_YL45                        2520                              2520
# 13       escherichia_coli_strain_Mt1B1                        5221                              5221
#########################################################################
#### This shows the integrated proteins file has all the unique proteins
