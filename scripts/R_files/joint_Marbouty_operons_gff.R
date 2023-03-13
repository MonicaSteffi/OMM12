joint_Marbouty_operons_gff=function()
{
  #----------------- KOFAMSCAN result for joint protein sequence file ----------------------------------------------------------------
  library("stringr")
  library("RCurl")
  library("data.table")
  
  ### Copied the results from "/dss/dssfs02/lwp-dss-0001/u7b03/u7b03-dss-0000/ra72seh/DB/Marbouty_annotations/OperonMapper/" to "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/OperonMapper/"
  kofamPath="S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/KOFAMSCAN_results_new/"
  gffPath="S:/Abilash/Marbouty_annotation_refseq_genbank_2022/gffFiles/"
  OperonPath="S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/OperonMapper_new/"
  cgcPath="S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/dbCAN2_new/"
  
  list_koFiles=list.files(path = kofamPath,full.names = T,pattern = "_kofam_results_eval1e03.rds")
  gnmNms=str_replace_all(string = list_koFiles,pattern = ".*/",replacement = "")
  gnmNms=str_replace_all(string = gnmNms,pattern = "_kofam.*",replacement = "")
  
  
  
  for(i2 in 1:length(gnmNms))
  {
    ##### Prokka_oligoMM_overall_gff:
    
    ov_gff=read.csv(file = str_c(gffPath,gnmNms[i2],"_genomic.gff"),header = F,sep = "\t",comment.char = "#")
    ov_gff=ov_gff[ov_gff$V3!="gene",]
    ov_ann=data.frame(Organism=ov_gff$V1)
    ov_ann$Start=ov_gff$V4
    ov_ann$Stop=ov_gff$V5
    ov_ann$Strand=ov_gff$V7
    ov_gff$V9=str_replace_all(string = ov_gff$V9,pattern = "NZ_",replacement = "")
    
    ov_ann$GeneID=ov_gff$V9
    ov_ann$GeneID=str_replace_all(string = ov_ann$GeneID,pattern = ";.*",replacement = "")
    ov_ann$GeneID=str_replace_all(string = ov_ann$GeneID,pattern = "ID=",replacement = "")
    ov_ann$GeneID=str_replace_all(string = ov_ann$GeneID,pattern = ".*cds-",replacement = "")
    ov_ann$GeneID=str_replace_all(string = ov_ann$GeneID,pattern = ".*rna-",replacement = "")
    ov_ann$GeneID=str_replace_all(string = ov_ann$GeneID,pattern = ".*exon-",replacement = "")
    
    
    productList=str_split(ov_gff$V9,pattern = ";")
    ov_ann$Product=unlist(lapply(productList,function(x){
      str_c(unique(unlist(x)[grep(pattern = "product_RefSeq|product_GenBank",x=unlist(x))]),collapse = ";")
    }))
    
    
    
    ov_ann$ProteinID=ov_gff$V9
    ov_ann$ProteinID=str_replace_all(string = ov_ann$ProteinID,pattern = ".*protein_id=",replacement = "")
    ov_ann$ProteinID=str_replace_all(string = ov_ann$ProteinID,pattern = ";.*",replacement = "")
    ov_ann$ProteinID=str_replace_all(string = ov_ann$ProteinID,pattern = ".*cds-",replacement = "")
    ov_ann$ProteinID=str_replace_all(string = ov_ann$ProteinID,pattern = ".*rna-",replacement = "")
    ov_ann$ProteinID=str_replace_all(string = ov_ann$ProteinID,pattern = ".*exon-",replacement = "")
    ov_ann$ProteinID=str_replace_all(string = ov_ann$ProteinID,pattern = ".*id-",replacement = "")
    
    
    
    ov_ann$gbkey=unlist(lapply(productList,function(x){
      str_c(unique(unlist(x)[grep(pattern = "gbkey_RefSeq|gbkey_GenBank|gbkey",x=unlist(x))]),collapse = ";")
    }))
    
    ov_ann$gene_biotype=unlist(lapply(productList,function(x){
      str_c(unique(unlist(x)[grep(pattern = "gene_biotype",x=unlist(x))]),collapse = ";")
    }))
    
    
    ov_ann$locus_tag=ov_gff$V9
    ov_ann$locus_tag=str_replace_all(string = ov_ann$locus_tag,pattern = ".*Locus_tag=",replacement = "")
    ov_ann$locus_tag=str_replace_all(string = ov_ann$locus_tag,pattern = ";.*",replacement = "")
    ov_ann$locus_tag=str_replace_all(string = ov_ann$locus_tag,pattern = ".*&",replacement = "")
    
    ov_ann$gene=ov_gff$V9
    ov_ann$gene=str_replace_all(string = ov_ann$gene,pattern = ".*gene=",replacement = "")
    ov_ann$gene=str_replace_all(string = ov_ann$gene,pattern = ";.*",replacement = "")
    ov_ann$gene=str_replace_all(string = ov_ann$gene,pattern = ".*cds-",replacement = "")
    ov_ann$gene=str_replace_all(string = ov_ann$gene,pattern = ".*rna-",replacement = "")
    ov_ann$gene=str_replace_all(string = ov_ann$gene,pattern = ".*id-",replacement = "")
    
    ov_ann$Note=unlist(lapply(productList,function(x){
      str_c(unique(unlist(x)[grep(pattern = "Note_RefSeq|Note_GenBank",x=unlist(x))]),collapse = ";")
    }))
    
    
    
    
    
    ov_ann$locus_tag_refseq=ov_gff$V9
    ov_ann$locus_tag_refseq=str_replace_all(string = ov_ann$locus_tag_refseq,pattern = ".*;locus_tag_RefSeq=",replacement = "")
    ov_ann$locus_tag_refseq=str_replace_all(string = ov_ann$locus_tag_refseq,pattern = ";.*",replacement = "")
    ov_ann$locus_tag_refseq=str_replace_all(string = ov_ann$locus_tag_refseq,pattern = ".*&",replacement = "")
    
    ov_ann$locus_tag_genbank=ov_gff$V9
    ov_ann$locus_tag_genbank=str_replace_all(string = ov_ann$locus_tag_genbank,pattern = ".*;locus_tag_GenBank=",replacement = "")
    ov_ann$locus_tag_genbank=str_replace_all(string = ov_ann$locus_tag_genbank,pattern = ";.*",replacement = "")
    ov_ann$locus_tag_genbank=str_replace_all(string = ov_ann$locus_tag_genbank,pattern = ".*&",replacement = "")
    
    
    
    for(i in 1:ncol(ov_ann))
    {
      ov_ann[,i]=str_replace_all(string = ov_ann[,i],pattern = "ID=",replacement = "")
      
    }
    
    # ov_ann$Name=ov_gff$V9
    # ov_ann$Name=str_replace_all(string = ov_ann$Name,pattern = ".*Name=",replacement = "")
    # ov_ann$Name=str_replace_all(string = ov_ann$Name,pattern = ";.*",replacement = "")
    # ov_ann$Name=str_replace_all(string = ov_ann$Name,pattern = "ID=",replacement = "")
    # ov_ann$Name[grep(pattern = "^ID_",x=ov_ann$Name)]=""
    
    
    ov_ann=ov_ann[which(!is.na(ov_ann$Start)),]
    ov_ann=ov_ann[-grep(pattern = "gene-",x=ov_ann$GeneID),]
    
    
    ov_ann$GeneID=str_replace_all(string = ov_ann$GeneID,pattern = ".*-",replacement = "")
    ov_ann=ov_ann[ov_ann$GeneID!="1",]
    ov_ann$str_Crd=apply(ov_ann[,c("Start"  ,"Stop","Strand")],1,function(x){str_c(str_trim(x,side = "both"),collapse = "_")})
    
    ov_ann_fin=ov_ann
    
    for(ov1 in 1: ncol(ov_ann_fin))
    {
      ov_ann_fin[,ov1]=str_replace_all(string = ov_ann_fin[,ov1],pattern = ".*=",replacement = "")
      
    }
    ov_ann_fin$locus_tag_refseq=str_replace_all(string = ov_ann_fin$locus_tag_refseq,pattern = ".*-",replacement = "")
    ov_ann_fin$locus_tag_genbank=str_replace_all(string = ov_ann_fin$locus_tag_genbank,pattern = ".*-",replacement = "")
    ov_ann_fin$locus_tag=str_replace_all(string = ov_ann_fin$locus_tag,pattern = ".*-",replacement = "")
    
    
    ov_ann_file=str_c("S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/",gnmNms[i2],"_gff_results.txt")
    write.table(ov_ann_fin,ov_ann_file,quote = F,row.names = F,sep = "\t")
    saveRDS(ov_ann_fin,file = str_c(ov_ann_file,".rds"))
    
    ######################################
    
    ##### The output folders are unzipped and named and kept in H:/PROKKA_04282020/Operon-mapper/
    library("stringr")
    
    list_operon_files=list.files(path=str_c(OperonPath,gnmNms[i2],"/"),recursive=T,full.names=T,pattern="list_of_operons")
    orf_coordinates_files=list.files(path = str_c(OperonPath,gnmNms[i2],"/"),recursive = T,full.names = T,pattern = "ORFs_coordinate")
    
    orf_coordinates=lapply(orf_coordinates_files,function(x){
      
      read.csv(x,header = F,sep = "\t")
    })
    orf_coordinates=do.call(rbind,orf_coordinates)
    colnames(orf_coordinates)=c("GenomeID","Source","Type","Start","End","Score","strand","Phase","Attributes")
    orf_coordinates[c('GeneID','Product')]=str_split_fixed(orf_coordinates$Attributes, ';', 2)
    orf_coordinates$GeneID=str_replace_all(string = orf_coordinates$GeneID,pattern = "ID=",replacement = "")
    orf_coordinates$Product=str_replace_all(string = orf_coordinates$Product,pattern = "product=",replacement = "")
    
    j=1
    
    operon_insdc=read.csv(file = list_operon_files[j],header = T,sep = "\t")
    operon_id_list=which(operon_insdc$Operon!="")
    
    for(i in 1:length(operon_id_list))
    {
      operon_insdc$Operon[operon_id_list[i]:nrow(operon_insdc)]=i
      
    }
    
    to_be_removed=intersect(which(is.na(operon_insdc$PosLeft)),which(is.na(operon_insdc$postRight)))
    operon_insdc=operon_insdc[-(to_be_removed),]
    operon_insdc$ID=str_c(gnmNms[i2],"_",operon_insdc$PosLeft,"_",operon_insdc$postRight)
    operon_insdc$Operon=str_c(gnmNms[i2],"_",operon_insdc$Operon)
    colnames(operon_insdc)[grep(pattern = "IdGene",x=colnames(operon_insdc))]="GeneID"
    operon_insdc$GeneID=str_replace_all(string = operon_insdc$GeneID,pattern = "cds-",replacement = "")
    operon_insdc$str_Crd=apply(operon_insdc[,c("PosLeft"  , "postRight" ,"Strand")],1,function(x){str_c(str_trim(x,side = "both"),collapse = "_")})
    operon_insdc$Operon_GenomeID=orf_coordinates$GenomeID[match(operon_insdc$GeneID,as.character(orf_coordinates$GeneID))]
    operon_insdc$Operon_product=orf_coordinates$Product[match(operon_insdc$GeneID,as.character(orf_coordinates$GeneID))]
    
    OperonfileName=str_c("S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/",gnmNms[i2],"_operonList_cleaned")
    write.table(file = OperonfileName,x = operon_insdc,quote = F,sep = "\t",row.names = F)
    saveRDS(operon_insdc,file = str_c(OperonfileName,".rds"))
    
  }

  
  
  
  }

joint_Marbouty_operons_gff()


