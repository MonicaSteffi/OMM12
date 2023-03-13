createGenbank=function()
{
  
  # if (!require("BiocManager"))
  #  install.packages("BiocManager")
  # BiocManager::install("GenomicRanges")
  library("GenomicRanges")
  # BiocManager::install("genbankr")
  library("genbankr")  
  library("stringr")
  # devtools::install_github("gschofl/biofiles")
  library("biofiles")
  library("stringr")
  library("seqinr")
  library("data.table")
  library("tidyverse")
  
  
  gnmFile=read.csv(file = "S:/Abilash/Marbouty_annotation_refseq_genbank 11_2021/GenomeID_to_Organism.csv",sep = ";",header = T)
  gbkFolder="S:/Abilash/Marbouty_annotation_refseq_genbank_2022//Data/Downloaded_from_genbank_and_refseq/Marbouty_genbank_format_files/"
  AnnFolder="S:/Abilash/Marbouty_annotation_refseq_genbank_2022//newJointAnnotation/"
  fin_colnames=c("seqnames" , "start" , "end" , "width" , "strand" , "type" , "locus_tag" , "loctype" ,"gene" ,"gene_synonym","gene_id","database","str_Crd")
  ## first, merge the gbk files from genbank and refseq.
  ## second, add the annotations to the merged gbk files.
  
  ############################################################################
  # Merging gbk files:--------------------------------------------------------
  ############################################################################
  gnmFile$Organism=str_replace_all(string = gnmFile$Organism,pattern = "Escheri",replacement = "escheri")
  gnmFile$Genome=str_replace_all(string = gnmFile$Genome,pattern = "Escheri",replacement = "escheri")
  feature_table_files=list.files(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Downloaded_from_genbank_and_refseq/Marbouty_Genbank/",pattern = "_feature_table.txt$",full.names = T,recursive = T)
  
  
  ##### Cleaned files.
  dbCAN2_files=list.files(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/",pattern = "_overall_cgc.rds",full.names = T)
  eggNOGDB_files=list.files(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/eggNOG_2022_new/eggNOGDBRes/",pattern = ".rds",full.names = T)
  NovelFamDB_files=list.files(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/eggNOG_2022_new/NovelFamDBRes/",pattern = ".rds",full.names = T)
  gff_files=list.files(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/",pattern = "_gff_results.txt.rds",full.names = T)
  kofamscan_files=list.files(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/KOFAMSCAN_results_new/",pattern = "_kofam_results_eval1e03.rds",full.names = T)
  operonMapper_files=list.files(pattern = "_operonList_cleaned.rds",path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/",full.names = T)
  PULs_files=list.files(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/PULs_new/",pattern = "pul.tsv$",full.names = T)
  RGI_files=list.files(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/RGI_5.1.1/RGI_for_genome/",pattern = "_genome_ann.txt",full.names = T)
  
  
  GenomeDesc=as.character(gnmFile$Organism)
  for(i in 1: length(GenomeDesc))
  {
    
    ### Read gbk files:
    gbkFile=tryCatch(readGenBank(file = str_c(gbkFolder,as.character(gnmFile$GenbankID[which(gnmFile$Organism==GenomeDesc[i])]),".gb")),error=function(x){NA})
    rfsqFile=tryCatch(readGenBank(file = str_c(gbkFolder,str_c("NZ_",as.character(gnmFile$GenbankID[which(gnmFile$Organism==GenomeDesc[i])])),".gb")),error=function(x){NA})
    
    # saveRDS(gbkFile,str_c("H:/",gnmFile$Organism[i],"_ori_gbkobj.rds"))
    # saveRDS(rfsqFile,str_c("H:/",gnmFile$Organism[i],"_ori_rfsqobj.rds"))
    
    ### Read annotation files:
    dbCAN2Res=readRDS(file = dbCAN2_files[grep(gnmFile$Genome[which(gnmFile$Organism==GenomeDesc[i])],x = dbCAN2_files)])
    kofamscanRes=readRDS(file = kofamscan_files[grep(gnmFile$Genome[which(gnmFile$Organism==GenomeDesc[i])],x = kofamscan_files)])
    gffRes=readRDS(file = gff_files[grep(gnmFile$Genome[which(gnmFile$Organism==GenomeDesc[i])],x = gff_files)])
    kofamscanRes$str_Crd=gffRes$str_Crd[match(kofamscanRes$GeneID,as.character(gffRes$locus_tag))]
    
    eggNOGDBRes=readRDS(file = eggNOGDB_files[grep(gnmFile$Genome[which(gnmFile$Organism==GenomeDesc[i])],x = eggNOGDB_files)])
    novelRes=tryCatch(readRDS(file = NovelFamDB_files[grep(gnmFile$Genome[which(gnmFile$Organism==GenomeDesc[i])],x = NovelFamDB_files)]),error=function(x){NA})
    
    operonRes=readRDS(file =operonMapper_files[grep(gnmFile$Genome[which(gnmFile$Organism==GenomeDesc[i])],x = operonMapper_files)] )
    RGIRes=read.csv(file = RGI_files[grep(gnmFile$Genome[which(gnmFile$Organism==GenomeDesc[i])],x = RGI_files)],header = T,sep = "\t")
    RGIRes$str_Crd=apply(RGIRes[,c("Start","Stop","Orientation")],1,function(x){str_c(str_trim(x,side = "both"),collapse = "_")})
    PULRes=tryCatch(read.csv(file = PULs_files[grep(gnmFile$Genome[which(gnmFile$Organism==GenomeDesc[i])],x = PULs_files)],header = T,sep = "\t"),error=function(x){NA})
    
    original_df=as.data.frame(genes(gbkFile))
    original_df$database="Genbank"
    original_df$str_Crd=apply(original_df[,c("start","end","strand")],1,function(x){str_c(str_trim(x,side = "both"),collapse = "_")})
    ### The gff file has the genes with two differnet coding regions are present in gbk(gbff)  and feature_table files, but NOT in gff file.
    #### To overcome this, we load feature_table and include these info to the gff files.
    ftRes=read.csv(file = feature_table_files[grep(gnmFile$Genome[which(gnmFile$Organism==GenomeDesc[i])],x = feature_table_files)],header = T,sep = "\t")
    
    extra_codingIdx_pos=c(intersect(which(apply(ftRes[,c("start","end")],1,function(x){x[1]>x[2]})),which(ftRes$strand=="+")),
                          intersect(which(apply(ftRes[,c("start","end")],1,function(x){x[1]>x[2]})),which(ftRes$strand=="-")))
    
    extra_coding_pos_locus_tag=as.character(unique(ftRes$locus_tag[extra_codingIdx_pos]))
    
    if(length(extra_coding_pos_locus_tag)>0)
    {
      
      for(exCoding in 1: length(extra_coding_pos_locus_tag))
      {
        locIdx=which(gffRes$locus_tag==extra_coding_pos_locus_tag[exCoding])
        if(length(locIdx)==1)
        {
          
          if(as.numeric(as.character(gffRes$Stop[locIdx]))==as.numeric(as.character(gffRes$Stop[1])))
          {
            newStart=1
            newStop=ftRes$end[extra_codingIdx_pos[exCoding]]
            newLocus=str_c(extra_coding_pos_locus_tag,"_Edge_start")  
            
            newgffDf=setNames(data.frame(matrix(ncol = ncol(gffRes), nrow = 1)), colnames(gffRes))
            newgffDf$Organism=gffRes$Organism[1]
            newgffDf$Start=1
            newgffDf$Stop=newStop
            newgffDf$Strand=ftRes$strand[extra_codingIdx_pos]
            newgffDf$GeneID=newLocus
            newgffDf$locus_tag=newLocus
            newgffDf$str_Crd=apply(newgffDf[,c("Start","Stop","Strand")],1,function(x){str_c(str_trim(x,side = "both"),collapse = "_")})
            } else {
            
              if((as.numeric(as.character(gffRes$Stop[1]))-as.numeric(as.character(gffRes$Stop[locIdx]))+1)>0)
              {
                newStart1=as.numeric(as.character(gffRes$Stop[1]))-as.numeric(as.character(gffRes$Stop[locIdx]))+1
                newStop2=as.numeric(as.character(gffRes$Stop[1]))
                newLocus1=str_c(extra_coding_pos_locus_tag,"_Edge_end")  
                
                newStart=1
                newStop=ftRes$end[extra_codingIdx_pos[exCoding]]
                newLocus=str_c(extra_coding_pos_locus_tag,"_Edge_start")  
                
                newgffDf=setNames(data.frame(matrix(ncol = ncol(gffRes), nrow = 2)), colnames(gffRes))
                
                newgffDf$Organism=rep(gffRes$Organism[1],2)
                newgffDf$Start=c(newStart1,newStart)
                newgffDf$Stop=c(newStop2,newStop)
                newgffDf$Strand=rep(ftRes$strand[extra_codingIdx_pos[1]],2)
                newgffDf$GeneID=rep(newLocus,2)
                newgffDf$locus_tag=rep(newLocus,2)
            
                newgffDf$str_Crd=apply(newgffDf[,c("Start","Stop","Strand")],1,function(x){str_c(str_trim(x,side = "both"),collapse = "_")})
                
              } else {
                
                newStart=1
                newStop=ftRes$end[extra_codingIdx_pos[exCoding]]
                newLocus=str_c(extra_coding_pos_locus_tag,"_Edge_start")  
                
                newgffDf=setNames(data.frame(matrix(ncol = ncol(gffRes), nrow = 1)), colnames(gffRes))
                
                newgffDf$Organism=gffRes$Organism[1]
                newgffDf$Start=newStart
                newgffDf$Stop=newStop
                newgffDf$Strand=ftRes$strand[extra_codingIdx_pos[1]]
                newgffDf$GeneID=newLocus
                newgffDf$locus_tag=newLocus
                newgffDf$str_Crd=apply(newgffDf[,c("Start","Stop","Strand")],1,function(x){str_c(str_trim(x,side = "both"),collapse = "_")})
                
              }
            
            
          }
            }
        
        
        
      }
      
    }
    
    gffRes=rbind(gffRes,newgffDf)
    gffRes$Stop[which(as.numeric(as.character(gffRes$Stop))>as.numeric(as.character(gffRes$Stop[1])))]=as.numeric(as.character(gffRes$Stop[1]))
    gffRes$database="gff_gb_n_rfsq"
    
    saveRDS(object = gffRes,
    file = str_c("S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/",gnmFile$Genome[which(gnmFile$Organism==GenomeDesc[i])],"_gffRes_Final.rds"))
    
    
    ### Format data to make GmakeRanges object:
    
    ##### 1. gff:
    
    if(!is.na(gffRes))
    {
      unqgff_GeneID=str_c(unique(gffRes$str_Crd))
      gffSht=gffRes[match(unqgff_GeneID,as.character(gffRes$str_Crd)),]
      # gffSht[c('start', 'end','strand')] <- str_split_fixed(as.character(gffSht$str_Crd), '_', 3)
      gffdf=data.frame(start=as.numeric(as.character(gffSht$Start)), end=as.numeric(as.character(gffSht$Stop)),strand=gffSht$Strand,
                       locus_tag=gffSht$locus_tag,gene_id=gffSht$str_Crd,
                       str_Crd=gffSht$str_Crd)
      gffdf$seqnames=rep(unique(original_df$seqnames)[1],nrow(gffdf))
      gffdf$type="gene"
      gffdf$loctype="normal"
      gffdf$pseudo=""
      gffdf$gene=""
      gffdf$gene_synonym=""
      gffdf$width=abs(as.numeric(as.character(gffdf$start))-as.numeric(as.character(gffdf$end)))+1
      gffdf$database="gff_gb_n_rfsq"
      
      original_df=original_df[match(setdiff(original_df$str_Crd,gffdf$str_Crd),as.character(original_df$str_Crd)),]
      dfFin=rbind(original_df[,fin_colnames],gffdf[,fin_colnames])
      
      
    }
    
    
    ##### 2. dbCAN2:
    if(!is.na(dbCAN2Res))
    {
      
      unqdbCAN2_GeneID=str_c(unique(dbCAN2Res$str_Crd))
      dbCAN2Sht=dbCAN2Res[match(unqdbCAN2_GeneID,as.character(dbCAN2Res$str_Crd)),]
      df=data.frame(start=dbCAN2Sht$Cluster_Start, end=dbCAN2Sht$Cluster_End,strand=dbCAN2Sht$Strand,
                    locus_tag=dbCAN2Sht$str_Crd,gene_id=dbCAN2Sht$str_Crd,str_Crd=dbCAN2Sht$str_Crd)
      df$seqnames=rep(unique(original_df$seqnames)[1],nrow(df))
      df$type="gene"
      df$loctype="normal"
      df$pseudo=""
      df$gene=""
      df$gene_synonym=""
      df$width=abs(df$start-df$end)+1
      df$database="dbCAN2"
      dfFin=rbind(dfFin,df[,fin_colnames])
      
    }
    
    ##### 3. kofamscan:
    
    if(!is.na(kofamscanRes))
    {
      
      unqkofamscan_GeneID=str_c(unique(kofamscanRes$str_Crd))
      kofamscanSht=kofamscanRes[match(unqkofamscan_GeneID,as.character(kofamscanRes$str_Crd)),]
      kofamscanSht[c('start', 'end','strand')] <- str_split_fixed(kofamscanSht$str_Crd, '_', 3)
      kofamdf=data.frame(start=kofamscanSht$start, end=kofamscanSht$end,strand=kofamscanSht$strand,
                         locus_tag=kofamscanSht$str_Crd,gene_id=kofamscanSht$str_Crd,str_Crd=kofamscanSht$str_Crd)
      kofamdf$seqnames=rep(unique(original_df$seqnames)[1],nrow(kofamdf))
      kofamdf$type="gene"
      kofamdf$loctype="normal"
      kofamdf$pseudo=""
      kofamdf$gene=""
      kofamdf$gene_synonym=""
      kofamdf$width=abs(as.numeric(as.character(kofamdf$start))-as.numeric(as.character(kofamdf$end)))+1
      kofamdf$database="kofamscan"
      dfFin=rbind(dfFin,kofamdf[,fin_colnames])
      
    }
    
    
    
    ##### 4. eggNOG:
    
    if(!is.na(eggNOGDBRes))
    {
      
      unqeggNOGDB_GeneID=str_c(unique(eggNOGDBRes$str_Crd))
      eggNOGDBSht=eggNOGDBRes[match(unqeggNOGDB_GeneID,as.character(eggNOGDBRes$str_Crd)),]
      eggNOGDBSht[c('start', 'end','strand')] <- str_split_fixed(eggNOGDBSht$str_Crd, '_', 3)
      eggNOGdf=data.frame(start=eggNOGDBSht$start, end=eggNOGDBSht$end,strand=eggNOGDBSht$strand,
                          locus_tag=eggNOGDBSht$str_Crd,gene_id=eggNOGDBSht$str_Crd,str_Crd=eggNOGDBSht$str_Crd)
      eggNOGdf$seqnames=rep(unique(original_df$seqnames)[1],nrow(eggNOGdf))
      eggNOGdf$type="gene"
      eggNOGdf$loctype="normal"
      eggNOGdf$pseudo=""
      eggNOGdf$gene=""
      eggNOGdf$gene_synonym=""
      eggNOGdf$width=abs(as.numeric(as.character(eggNOGdf$start))-as.numeric(as.character(eggNOGdf$end)))+1
      eggNOGdf$database="eggNOG5"
      dfFin=rbind(dfFin,eggNOGdf[,fin_colnames])
      
    }
    
    
    ##### 5. PULs:
    
    if(!is.na(PULRes))
    {
      PULRes$str_Crd=apply(PULRes[,c("start","end","strand")],1,function(x){str_c(x,collapse = "_")})
      unqPULRes_GeneID=str_c(unique(PULRes$str_Crd))
      PULResSht=PULRes[match(unqPULRes_GeneID,as.character(PULRes$str_Crd)),]
      PULdf=data.frame(start=PULResSht$start, end=PULResSht$end,strand=PULResSht$strand,
                       locus_tag=PULResSht$str_Crd,gene_id=PULResSht$str_Crd,str_Crd=PULResSht$str_Crd)
      PULdf$seqnames=rep(unique(original_df$seqnames)[1],nrow(PULdf))
      PULdf$type="gene"
      PULdf$loctype="normal"
      PULdf$pseudo=""
      PULdf$gene=""
      PULdf$gene_synonym=""
      PULdf$width=abs(as.numeric(as.character(PULdf$start))-as.numeric(as.character(PULdf$end)))+1
      PULdf$database="PUL"
      dfFin=rbind(dfFin,PULdf[,fin_colnames])
      
      
    }
    
    
    ##### 6. RGI:
    if(!is.na(RGIRes))
    {
      
      unqRGIRes_GeneID=str_c(unique(RGIRes$str_Crd))
      RGIResSht=RGIRes[match(unqRGIRes_GeneID,as.character(RGIRes$str_Crd)),]
      RGIResSht[c('start', 'end','strand')] <- str_split_fixed(RGIResSht$str_Crd, '_', 3)
      RGIdf=data.frame(start=RGIResSht$start, end=RGIResSht$end,strand=RGIResSht$strand,
                       locus_tag=RGIResSht$str_Crd,gene_id=RGIResSht$str_Crd,str_Crd=RGIResSht$str_Crd)
      
      RGIdf$seqnames=rep(unique(original_df$seqnames)[1],nrow(RGIdf))
      RGIdf$type="gene"
      RGIdf$loctype="normal"
      RGIdf$pseudo=""
      RGIdf$gene=""
      RGIdf$gene_synonym=""
      RGIdf$width=abs(as.numeric(as.character(RGIdf$start))-as.numeric(as.character(RGIdf$end)))+1
      RGIdf$database="RGI"
      dfFin=rbind(dfFin,RGIdf[,fin_colnames])
      
    }
    
    
    ##### 7. novelRes:
    
    if(!is.na(novelRes))
    {
      
      unqnovelRes_GeneID=str_c(unique(novelRes$str_Crd))
      novelResSht=novelRes[match(unqnovelRes_GeneID,as.character(novelRes$str_Crd)),]
      novelResSht[c('start', 'end','strand')] <- str_split_fixed(novelResSht$str_Crd, '_', 3)
      novelResdf=data.frame(start=novelResSht$start, end=novelResSht$end,strand=novelResSht$strand,
                            locus_tag=novelResSht$str_Crd,gene_id=novelResSht$str_Crd,str_Crd=novelResSht$str_Crd)
      novelResdf$seqnames=rep(unique(original_df$seqnames)[1],nrow(novelResdf))
      novelResdf$type="gene"
      novelResdf$loctype="normal"
      novelResdf$pseudo=""
      novelResdf$gene=""
      novelResdf$gene_synonym=""
      novelResdf$width=abs(as.numeric(as.character(novelResdf$start))-as.numeric(as.character(novelResdf$end)))+1
      novelResdf$database="novelfam"
      dfFin=rbind(dfFin,novelResdf[,fin_colnames])
      
    }
    
    
    ##### 8. operon:
    if(!is.na(operonRes))
    {
      
      unqoperon_GeneID=str_c(unique(operonRes$str_Crd))
      operonSht=operonRes[match(unqoperon_GeneID,as.character(operonRes$str_Crd)),]
      operonSht[c('start', 'end','strand')] <- str_split_fixed(operonSht$str_Crd, '_', 3)
      operondf=data.frame(start=operonSht$start, end=operonSht$end,strand=operonSht$strand,
                          locus_tag=operonSht$str_Crd,gene_id=operonSht$str_Crd,str_Crd=operonSht$str_Crd)
      operondf$seqnames=rep(unique(original_df$seqnames)[1],nrow(operondf))
      operondf$type="gene"
      operondf$loctype="normal"
      operondf$pseudo=""
      operondf$gene=""
      operondf$gene_synonym=""
      operondf$width=abs(as.numeric(as.character(operondf$start))-as.numeric(as.character(operondf$end)))+1
      operondf$database="operonMapper"
      dfFin=rbind(dfFin,operondf[,fin_colnames])
      
    }
    
    # rfsqGenes=tryCatch(expr = as.data.frame(genes(rfsqFile)),error=function(x){NA})
    
  
    dfFin=dfFin [!duplicated(dfFin[c(2,3,5,7,13)]),]
    dfFin=dfFin[!duplicated(dfFin[,"str_Crd",])]
    dfFin$start=as.numeric(as.character(dfFin$start))
    dfFin$end=as.numeric(as.character(dfFin$end))
    dfFin$width=as.numeric(as.character(dfFin$width))
    
    na_rows=which(is.na(dfFin$start))
    if(length(na_rows)>0)
    {
      dfFin=dfFin[-(na_rows),]
      
    }
    
    # dfFin$str_Crd=apply(dfFin[,c("start","end","strand")],1,function(x){str_c(str_trim(x,side = "both"),collapse = "_")})
    
    dfFin=(dfFin[
with(dfFin, order(as.numeric(as.character(start)), as.numeric(as.character(end)))),
])
  
########## There are many genes with same end coordinate (stop codon) , but slightly different start codon (we consider the difference <=50bp) & same strand to be the same gene.
########## We take one column out of the two:
    ########## But in the one column that is retained, we replaced the end,str_crd with the information from longest gene and geneID, locus_tag with the information with alphabets
    
    stpTable=table(dfFin$end)
    sameStopIdx=which(stpTable>1)
    
    
    
    for(stCnt in 1: length(sameStopIdx))
    {
      
      tmpStp=names(stpTable)[sameStopIdx[stCnt]]  
      
      dfFin_idx=which(dfFin$end==tmpStp)
      
      if(length(dfFin_idx)==2)
      {
        
        if(abs(dfFin$start[dfFin_idx[1]]-dfFin$start[dfFin_idx[2]]) <=50 && dfFin$strand[dfFin_idx[1]] == dfFin$strand[dfFin_idx[2]])
        {
          
          
          retLocus_tag=grep(pattern = "^[A-z]",x=dfFin$locus_tag[dfFin_idx])  
          retloc_fin=dfFin$locus_tag[dfFin_idx[retLocus_tag]]
          retStart=as.numeric(min(dfFin$start[dfFin_idx]))
          retLength=as.numeric(max(dfFin$width[dfFin_idx]))
          retDB=str_c(dfFin$database[dfFin_idx],collapse = "/")
          retStrCrd=str_c(retStart,"_",dfFin$end[dfFin_idx[1]],"_",dfFin$strand[dfFin_idx[1]])
          
          dfFin$locus_tag[dfFin_idx[retLocus_tag]]=retloc_fin
          dfFin$gene_id[dfFin_idx[retLocus_tag]]=retloc_fin
          
          dfFin$start[dfFin_idx[retLocus_tag]]=retStart
          dfFin$width[dfFin_idx[retLocus_tag]]=retLength
          dfFin$database[dfFin_idx[retLocus_tag]]=retDB
          dfFin$str_Crd[dfFin_idx[retLocus_tag]]=retStrCrd
          
          
          tbr_rows=c(tbr_rows,setdiff(dfFin_idx,dfFin_idx[retLocus_tag])) 
          
          
          
        }
      } else {
        
        startDiff=dfFin$start[dfFin_idx]-50
        
        strndInfo=as.character(dfFin$strand[dfFin_idx])
        tt=table(strndInfo)
        dfFin_idx=dfFin_idx[grep(names(which(tt==max(tt))),x=strndInfo)]
        
        dfFin_idx=intersect(dfFin_idx[which(startDiff<=0)],dfFin_idx)
        
        if(length(dfFin_idx)>1)
        {
          retLocus_tag=grep(pattern = "^[A-z]",x=dfFin$locus_tag[dfFin_idx])  
          retloc_fin=dfFin$locus_tag[dfFin_idx[retLocus_tag]]
          retStart=as.numeric(min(dfFin$start[dfFin_idx]))
          retLength=as.numeric(max(dfFin$width[dfFin_idx]))
          retDB=str_c(dfFin$database[dfFin_idx],collapse = "/")
          retStrCrd=str_c(retStart,"_",dfFin$end[dfFin_idx[1]],"_",dfFin$strand[dfFin_idx[1]])
          
          dfFin$locus_tag[dfFin_idx[retLocus_tag]]=retloc_fin
          dfFin$gene_id[dfFin_idx[retLocus_tag]]=retloc_fin
          
          dfFin$start[dfFin_idx[retLocus_tag]]=retStart
          dfFin$width[dfFin_idx[retLocus_tag]]=retLength
          dfFin$database[dfFin_idx[retLocus_tag]]=retDB
          dfFin$str_Crd[dfFin_idx[retLocus_tag]]=retStrCrd
          
          
          tbr_rows=c(tbr_rows,setdiff(dfFin_idx,dfFin_idx[retLocus_tag])) 
          
          
          
        }
        
        
      }
      
      
    }
    
    dfFin=dfFin[-tbr_rows,]
    
    
    
    subMat=matrix(0,nrow=nrow(dfFin),ncol=nrow(dfFin))
    rownames(subMat)=as.character(dfFin$str_Crd)
    colnames(subMat)=as.character(dfFin$str_Crd)
    genomeIdx=intersect(which(dfFin$start==1),which(dfFin$end==max(dfFin$end)))
    dfFinv2
    
    for(dfCnt in 1: nrow(dfFin))
    {
      
      
      
    }
    
    
    
    
    # GRanges_obj=makeGRangesFromDataFrame(dfFin,
    #                          keep.extra.columns=TRUE,
    #                          ignore.strand=FALSE,
    #                          seqinfo=NULL,
    #                          seqnames.field=c("seqnames", "seqname",
    #                                           "chromosome", "chrom",
    #                                           "chr", "chromosome_name",
    #                                           "seqid"),
    #                          start.field="start",
    #                          end.field="end",
    #                          strand.field="strand",
    #                          starts.in.df.are.0based=FALSE)
    # 
    # 
    # gbkFile@genes=GRanges_obj
    # gbkFile@cds=GRanges_obj
    # gbkFile@exons=GRanges_obj
    # gbkFile@transcripts=GRanges_obj
    # 
    #  GRanges_obj=makeGRangesFromDataFrame(dfFin,
    #                                       keep.extra.columns=TRUE,
    #                                       ignore.strand=FALSE,
    #                                       seqinfo=NULL,
    #                                       seqnames.field=c("seqnames", "seqname",
    #                                                        "chromosome", "chrom",
    #                                                        "chr", "chromosome_name",
    #                                                        "seqid"),
    #                                       start.field="start",
    #                                       end.field=c("end", "stop"),
    #                                       strand.field="strand",
    #                                       starts.in.df.are.0based=FALSE)
    
    ########################################################################################################################################
    ## Now, remove the subset genes: Ex: the genes starting from 3end genome to 5end beginning is present as two sets- 
    ### (i) 5end region of the genome.
    ### (ii) 3end region of the genome.
    ## eggNOG and operonMapper predicted the scenario-(i) to be different protein.
    ## It is safe to remove these  proteins from eggNOG and operonmapper which are subset of the scenario-(i)
    ########################################################################################################################################    
    
    remSubsets=function(dfFin)
    {
      coord_list=list()
      coord_vec=vector()
      for(crdCnt in 1: nrow(dfFin))
      {
        CrdVec=as.vector(seq(as.numeric(as.character(dfFin$start[crdCnt])),as.numeric(as.character(dfFin$end[crdCnt])),by=1))
        
        if(dfFin$strand[crdCnt]=="-")
        {
          
          CrdVec=rev(CrdVec)
        }
        coord_list[[crdCnt]]=CrdVec
        coord_vec[crdCnt]=str_c(CrdVec,collapse = "_")
        }
      
      
      
      dfFin$coord_vec=coord_vec
    
      for(y1 in 1: (nrow(dfFin)-1))
      {
        
        
        
      }
        
    
      }
    
    
    
    saveRDS(dfFin,str_c("H:/",gnmFile$Organism[i],"_fin_featuretable.rds"))
    ############################################################################
    # Add annotations to merged gbk file:---------------------------------------
    ############################################################################
  }
  
  
  
  }