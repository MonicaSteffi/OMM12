Marbouty_ann=function()
{
  #----------------- KOFAMSCAN result for joint protein sequence file ----------------------------------------------------------------
  library("stringr")
  library("RCurl")
  library("data.table")
  
  kofamPath="S:/Abilash/Marbouty_annotation_refseq_genbank 11_2021/Data/Marbouty_annotations_for_integrated_proteins/KOFAMSCAN_results/"
gffPath="S:/Abilash/Marbouty_annotation_refseq_genbank 11_2021/Data/Marbouty_annotations_for_integrated_proteins/gffFiles/"
OperonPath="S:/Abilash/Marbouty_annotation_refseq_genbank 11_2021/Data/Marbouty_annotations_for_integrated_proteins/OperonMapper/"
cgcPath="S:/Abilash/Marbouty_annotation_refseq_genbank 11_2021/Data/Marbouty_annotations_for_integrated_proteins/dbCAN2/"

  
  list_koFiles=list.files(path = kofamPath,full.names = T,pattern = "_result_all_tab_edited")
  gnmNms=str_replace_all(string = list_koFiles,pattern = ".*//",replacement = "")
  gnmNms=str_replace_all(string = gnmNms,pattern = "_result.*",replacement = "")
  
  
  for(i2 in 1:length(gnmNms))
  {
    library("data.table")
    kofam_b_coc=fread(file =  list_koFiles[i2],header = F,sep = "\t")
    colnames(kofam_b_coc)=c("F_score","GeneID","KO","threshold","score","evalue","KO_definition")
    
    star_genes=unique(kofam_b_coc$GeneID[which(kofam_b_coc$F_score=="*")])
    kofam_star=kofam_b_coc[match(star_genes,kofam_b_coc$GeneID),]
    
    all_unq_genes=unique(kofam_b_coc$GeneID)
    nostar_genes=setdiff(all_unq_genes,star_genes)
    kofam_nostar=kofam_b_coc[match(nostar_genes,kofam_b_coc$GeneID),]
    unq_gn=unique(kofam_nostar$GeneID)
  
    for(i in 1:length(unq_gn))
    {
      tmp_idx=which(kofam_nostar$GeneID==unq_gn[i])
      
      if(length(tmp_idx)>1)
      {
        fn_idx=tmp_idx[which(kofam_nostar$evalue[tmp_idx]==min(kofam_nostar$evalue[tmp_idx]))]
        
        
      } else {
        fn_idx=tmp_idx
        
      }
      
      if(i==1)
      {
        kofam_nostar_fin=kofam_nostar[fn_idx,]
        
      } else {
        
        kofam_nostar_fin=rbind(kofam_nostar_fin,kofam_nostar[fn_idx,])
      }
    }
    
    kofam_nostar_fin=kofam_nostar_fin[kofam_nostar_fin$evalue<=1e-03,]
    
    unq_gn=unique(kofam_star$GeneID)
    for(i in 1:length(unq_gn))
    {
      tmp_idx=which(kofam_star$GeneID==unq_gn[i])
      
      if(length(tmp_idx)>1)
      {
        fn_idx=tmp_idx[which(kofam_star$evalue[tmp_idx]==min(kofam_star$evalue[tmp_idx]))]
        
        
      } else {
        fn_idx=tmp_idx
        
      }
      
      if(i==1)
      {
        kofam_star_fin=kofam_star[fn_idx,]
        
      } else {
        
        kofam_star_fin=rbind(kofam_star_fin,kofam_star[fn_idx,])
      }
    }
    
    kofam_star_fin=kofam_star_fin[kofam_star_fin$evalue<=1e-03,]
    
    kofam_b_coc_fin=rbind(kofam_star_fin,kofam_nostar_fin)
    fin_file_name=str_c(kofamPath,gnmNms[i2],"_kofam_results_eval1e03")
    write.table(x = kofam_b_coc_fin,quote=F,sep="\t",file = fin_file_name,row.names = F)  
    saveRDS(kofam_b_coc_fin,file = str_c(fin_file_name,".rds"))
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
    
    
    ov_ann$Product=ov_gff$V9
    ov_ann$Product=str_replace_all(string = ov_ann$Product,pattern = ".*product=",replacement = "")
    ov_ann$Product=str_replace_all(string = ov_ann$Product,pattern = ";.*",replacement = "")
    ov_ann$Product=str_replace_all(string = ov_ann$Product,pattern = ".*cds-",replacement = "")
    ov_ann$Product=str_replace_all(string = ov_ann$Product,pattern = ".*rna-",replacement = "")
    ov_ann$Product=str_replace_all(string = ov_ann$Product,pattern = ".*exon-",replacement = "")
    ov_ann$Product=str_replace_all(string = ov_ann$Product,pattern = ".*id-",replacement = "")
    
    
    ov_ann$ProteinID=ov_gff$V9
    ov_ann$ProteinID=str_replace_all(string = ov_ann$ProteinID,pattern = ".*protein_id=",replacement = "")
    ov_ann$ProteinID=str_replace_all(string = ov_ann$ProteinID,pattern = ";.*",replacement = "")
    ov_ann$ProteinID=str_replace_all(string = ov_ann$ProteinID,pattern = ".*cds-",replacement = "")
    ov_ann$ProteinID=str_replace_all(string = ov_ann$ProteinID,pattern = ".*rna-",replacement = "")
    ov_ann$ProteinID=str_replace_all(string = ov_ann$ProteinID,pattern = ".*exon-",replacement = "")
    ov_ann$ProteinID=str_replace_all(string = ov_ann$ProteinID,pattern = ".*id-",replacement = "")
    
    
    ov_ann$RptFamily=ov_gff$V9
    ov_ann$RptFamily=str_replace_all(string = ov_ann$RptFamily,pattern = ".*rpt_family=",replacement = "")
    ov_ann$RptFamily=str_replace_all(string = ov_ann$RptFamily,pattern = ";.*",replacement = "")
    ov_ann$RptFamily=str_replace_all(string = ov_ann$RptFamily,pattern = ".*cds-",replacement = "")
    ov_ann$RptFamily=str_replace_all(string = ov_ann$RptFamily,pattern = ".*rna-",replacement = "")
    ov_ann$RptFamily=str_replace_all(string = ov_ann$RptFamily,pattern = ".*exon-",replacement = "")
    ov_ann$RptFamily=str_replace_all(string = ov_ann$RptFamily,pattern = ".*id-",replacement = "")
    
    ov_ann$gbkey=ov_gff$V9
    ov_ann$gbkey=str_replace_all(string = ov_ann$gbkey,pattern = ".*gbkey=",replacement = "")
    ov_ann$gbkey=str_replace_all(string = ov_ann$gbkey,pattern = ";.*",replacement = "")
    ov_ann$gbkey=str_replace_all(string = ov_ann$gbkey,pattern = ".*cds-",replacement = "")
    ov_ann$gbkey=str_replace_all(string = ov_ann$gbkey,pattern = ".*rna-",replacement = "")
    ov_ann$gbkey=str_replace_all(string = ov_ann$gbkey,pattern = ".*exon-",replacement = "")
    ov_ann$gbkey=str_replace_all(string = ov_ann$gbkey,pattern = ".*id-",replacement = "")
    
    ov_ann$Dbxref=ov_gff$V9
    ov_ann$Dbxref=str_replace_all(string = ov_ann$Dbxref,pattern = ".*Dbxref=",replacement = "")
    ov_ann$Dbxref=str_replace_all(string = ov_ann$Dbxref,pattern = ";.*",replacement = "")
    ov_ann$Dbxref=str_replace_all(string = ov_ann$Dbxref,pattern = ".*cds-",replacement = "")
    ov_ann$Dbxref=str_replace_all(string = ov_ann$Dbxref,pattern = ".*rna-",replacement = "")
    ov_ann$Dbxref=str_replace_all(string = ov_ann$Dbxref,pattern = ".*exon-",replacement = "")
    ov_ann$Dbxref=str_replace_all(string = ov_ann$Dbxref,pattern = ".*id-",replacement = "")
    
    
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
    
    ov_ann$Note=ov_gff$V9
    ov_ann$Note=str_replace_all(string = ov_ann$Note,pattern = ".*Note=",replacement = "")
    ov_ann$Note=str_replace_all(string = ov_ann$Note,pattern = ";.*",replacement = "")
    ov_ann$Note[grep(pattern = "^ID_",x=ov_ann$Note)]=""
    ov_ann$Note=str_replace_all(string = ov_ann$Note,pattern = ".*cds-",replacement = "")
    ov_ann$Note=str_replace_all(string = ov_ann$Note,pattern = ".*rna-",replacement = "")
    ov_ann$Note=str_replace_all(string = ov_ann$Note,pattern = ".*exon-",replacement = "")
    ov_ann$Note=str_replace_all(string = ov_ann$Note,pattern = ".*id-",replacement = "")
    
    ov_ann$regulatory_class=ov_gff$V9
    ov_ann$regulatory_class=str_replace_all(string = ov_ann$regulatory_class,pattern = ".*regulatory_class=",replacement = "")
    ov_ann$regulatory_class=str_replace_all(string = ov_ann$regulatory_class,pattern = ";.*",replacement = "")
    ov_ann$regulatory_class[grep(pattern = "^ID_",x=ov_ann$regulatory_class)]=""
    
    
    ov_ann$bound_moiety=ov_gff$V9
    ov_ann$bound_moiety=str_replace_all(string = ov_ann$bound_moiety,pattern = ".*bound_moiety=",replacement = "")
    ov_ann$bound_moiety=str_replace_all(string = ov_ann$bound_moiety,pattern = ";.*",replacement = "")
    ov_ann$bound_moiety[grep(pattern = "^ID_",x=ov_ann$bound_moiety)]=""
    # ov_ann$OtherGeneID=ov_gff$V9
    # ov_ann$OtherGeneID=str_replace_all(string = ov_ann$OtherGeneID,pattern = ".*ID=",replacement = "")
    # ov_ann$OtherGeneID=str_replace_all(string = ov_ann$OtherGeneID,pattern = ";.*",replacement = "")
    
    for(i in 1:ncol(ov_ann))
    {
      ov_ann[,i]=str_replace_all(string = ov_ann[,i],pattern = "ID=.*",replacement = "")
      
    }
    
    # ov_ann$Name=ov_gff$V9
    # ov_ann$Name=str_replace_all(string = ov_ann$Name,pattern = ".*Name=",replacement = "")
    # ov_ann$Name=str_replace_all(string = ov_ann$Name,pattern = ";.*",replacement = "")
    # ov_ann$Name=str_replace_all(string = ov_ann$Name,pattern = "ID=",replacement = "")
    # ov_ann$Name[grep(pattern = "^ID_",x=ov_ann$Name)]=""
    
    
    ov_ann=ov_ann[!is.na(ov_ann$Start),]
    ov_ann=ov_ann[-grep(pattern = "gene-",x=ov_ann$GeneID)]
    
    
    ov_ann$GeneID=str_replace_all(string = ov_ann$GeneID,pattern = ".*-",replacement = "")
    ov_ann=ov_ann[ov_ann$GeneID!="1",]
    ov_ann_kofam_merge=merge(ov_ann,kofam_b_coc_fin,by="GeneID",all = TRUE)
    ov_ann_kofam_merge$ID=str_c(gnmNms[i2],"_",ov_ann_kofam_merge$Start,"_",ov_ann_kofam_merge$Stop)
    ov_ann_kofam_merge_file=str_c("/dss/dssfs02/lwp-dss-0001/u7b03/u7b03-dss-0000/ra72seh/DB/Marbouty_joint/",gnmNms[i2],"_kofam_gff_info_merged_results.txt")
    write.table(ov_ann_kofam_merge,ov_ann_kofam_merge_file,quote = F,row.names = F,sep = "\t")
    saveRDS(ov_ann_kofam_merge,file = str_c(ov_ann_kofam_merge_file,".rds"))
    
    ######################################
    # --------- Operon-mapper annotations ----------------------------------------
    
    ##### Create separate genome files from the aggregated multi-fasta file to upload it to operon-mapper portal.
    # library("seqinr")
    # proteinFiles=list.files(path = ProteinPath,pattern = "")
    # seq_genomes=read.fasta(file = "H:/PROKKA_04282020/PROKKA_04282020.fna",seqtype = "DNA",as.string = TRUE, set.attributes = FALSE)
    # unq_gnm_names=names(seq_genomes)
    
    # for(i in 1:length(unq_gnm_names))
    # {
    #  tmp_seq=seq_genomes[[unq_gnm_names[i]]]
    #  tmp_seq=toupper(tmp_seq) # [unq_gnm_names[i]]
    #  write.fasta(sequences = tmp_seq,names = unq_gnm_names[i],file.out = str_c("H:/PROKKA_04282020/PROKKA_",unq_gnm_names[i],".fasta"))
    
    
    # }
    
    ##### The output folders are unzipped and named and kept in H:/PROKKA_04282020/Operon-mapper/
    library("stringr")
    
    list_operon_files=list.files(path=str_c(OperonPath,gnmNms[i2],"/"),recursive=T,full.names=T,pattern="list_of_operons")
    
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
    OperonfileName=str_c("/dss/dssfs02/lwp-dss-0001/u7b03/u7b03-dss-0000/ra72seh/DB/Marbouty_joint/",gnmNms[i2],"_operonList_cleaned")
    write.table(file = OperonfileName,x = operon_insdc,quote = F,sep = "\t",row.names = F)
    saveRDS(operon_insdc,file = str_c(OperonfileName,".rds"))
    
    kofam_operon_ann=merge(ov_ann_kofam_merge,operon_insdc,by="ID",all=TRUE)
    kofam_operon_gff_mergedFileName=str_c("/dss/dssfs02/lwp-dss-0001/u7b03/u7b03-dss-0000/ra72seh/DB/Marbouty_joint/",gnmNms[i2],"_kofam_op_gff_merged.csv")
    write.table(kofam_operon_ann,kofam_operon_gff_mergedFileName,quote = F,sep = "\t",row.names = F)
    saveRDS(kofam_operon_ann,file = str_c(kofam_operon_gff_mergedFileName,".rds"))
    
    list_cgc_files=list.files(path=str_c(cgcPath,gnmNms[i2],"/"),recursive=T,full.names=T,pattern="cgc.out.txt|cgc_cazy_tc_tf.out.txt|cgc_cazy_tf.out.txt")
    message(str_c("dbCAN is starting!!"))
    source("/dss/dssfs02/lwp-dss-0001/u7b03/u7b03-dss-0000/ra72seh/joint_Marbouty_dbCAN.R")
    overall_cgc=joint_Marbouty_dbCAN(list_cgc_files = list_cgc_files,gnmNmsThis = gnmNms[i2])
    message(str_c("dbCAN is done!!"))
    
    saveRDS(overall_cgc,str_c("/dss/dssfs02/lwp-dss-0001/u7b03/u7b03-dss-0000/ra72seh/DB/Marbouty_joint/",gnmNms[i2],"_overall_cgc.rds"))
    
    overall_ann=merge(overall_cgc,kofam_operon_ann,by="ID",all=TRUE)
    message("merging all files done!")
    
    overall_annFileName=str_c("/dss/dssfs02/lwp-dss-0001/u7b03/u7b03-dss-0000/ra72seh/DB/Marbouty_joint/",gnmNms[i2],"_kofam_op_gff_cgc_merged.csv")
    message("naming final file done!!")
    write.table(overall_ann,overall_annFileName,quote = F,sep = "\t",row.names = F)
    message("writing the overall annotation file done!!")
    saveRDS(overall_ann,file = str_c(overall_annFileName,".rds"))
    message("saved overall annotation rds file!")
    
    
     }
}

Marbouty_ann()


