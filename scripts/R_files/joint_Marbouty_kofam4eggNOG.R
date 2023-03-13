joint_Marbouty_kofam4eggNOG=function()
{
  #----------------- KOFAMSCAN result for joint protein sequence file ----------------------------------------------------------------
  ### This is run after S:\Abilash\Marbouty_annotation_refseq_genbank 11_2021\scripts\Bash_files\02_kofamscan_scripts\Marbouty_joint_kofamscan_analysis_editFormat.sh
  ### The results from this script is used as input for this script
  library("stringr")
  library("RCurl")
  library("data.table")
  library("foreach")
  library("dplyr")
  
  kofamPath="S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/eggNOG_2022_new/kofamscan_for_eggNOG_pred_proteins/"
  
  list_koFiles=list.files(path = kofamPath,full.names = T,pattern = "_result_all_tab_edited")
  gnmNms=str_replace_all(string = list_koFiles,pattern = ".*/",replacement = "")
  gnmNms=str_replace_all(string = gnmNms,pattern = "_result.*",replacement = "")
  
  
  for(i2 in 1:length(gnmNms))
  {
    
    kofam_b_coc=fread(file =  list_koFiles[i2],header = F,sep = "\t")
    colnames(kofam_b_coc)=c("F_score","GeneID","KO","threshold","score","evalue","KO_definition")
    kofam_b_coc$evalue=as.numeric(kofam_b_coc$evalue)
    
    star_genes=unique(kofam_b_coc$GeneID[which(kofam_b_coc$F_score=="*")])
    kofam_star=kofam_b_coc[kofam_b_coc$F_score=="*",]
    
    all_unq_genes=unique(kofam_b_coc$GeneID)
    nostar_genes=setdiff(all_unq_genes,star_genes)
    kofam_nostar=kofam_b_coc[kofam_b_coc$GeneID %in% nostar_genes,]
    unq_gn=unique(kofam_nostar$GeneID)
    
    ## Take those KOs with most minimal evalue, for each gene
    
    kofam_nostar_fin=(kofam_nostar %>% group_by(GeneID) %>% slice(which.min(evalue)))
    kofam_nostar_fin$evalue=as.numeric(kofam_nostar_fin$evalue)
    
    ## Remove those rows with evalues < 1e-03
    kofam_nostar_fin=kofam_nostar_fin[kofam_nostar_fin$evalue<=1e-03,]
    
    kofam_star_fin=(kofam_star %>% group_by(GeneID) %>% slice(which.min(evalue)))
    kofam_star_fin$evalue=as.numeric(kofam_star_fin$evalue)
    
    
    
    kofam_b_coc_fin=rbind(kofam_star_fin,kofam_nostar_fin)
    kofam_b_coc_fin=kofam_b_coc_fin[kofam_b_coc_fin$KO!="-",]
    kofam_b_coc_fin=kofam_b_coc_fin %>% arrange(GeneID)
    
    fin_file_name=str_c(kofamPath,gnmNms[i2],"_kofam_results_eval1e03")
    write.table(x = kofam_b_coc_fin,quote=F,sep="\t",file = fin_file_name,row.names = F)  
    saveRDS(kofam_b_coc_fin,file = str_c(fin_file_name,".rds"))
    
    
    
    
    
    
  }
}

joint_Marbouty_kofam4eggNOG()


