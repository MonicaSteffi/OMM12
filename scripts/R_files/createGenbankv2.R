createGenbankv2=function()
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
  
  ## first, merge the gbk files from genbank and refseq.
  ## second, add the annotations to the merged gbk files.
  
  ############################################################################
  # Merging gbk files:--------------------------------------------------------
  ############################################################################
  for(i in 1: nrow(gnmFile))
  {
    ### Read gbk files:
    gbkFile=tryCatch(readGenBank(file = str_c(gbkFolder,gnmFile$GenbankID[i],".gb")),error=function(x){NA})
    rfsqFile=tryCatch(readGenBank(file = str_c(gbkFolder,gnmFile$RefSeqID[i],".gb")),error=function(x){NA})
    
    
    ### Read annotation files:
    dbCAN2Res=readRDS(file = dbCAN2_files[grep(GenomeNames[i],x = dbCAN2_files)])
    kofamscanRes=readRDS(file = kofamscan_files[grep(GenomeNames[i],x = kofamscan_files)])
    gffRes=readRDS(file = gff_files[grep(GenomeNames[i],x = gff_files)])
    kofamscanRes$str_Crd=gffRes$str_Crd[match(kofamscanRes$GeneID,as.character(gffRes$locus_tag))]
    
    eggNOGDBRes=readRDS(file = eggNOGDB_files[grep(GenomeNames[i],x = eggNOGDB_files)])
    novelRes=tryCatch(readRDS(file = NovelFamDB_files[grep(GenomeNames[i],x = NovelFamDB_files)]),error=function(x){NA})
    
    operonRes=readRDS(file =operonMapper_files[grep(GenomeNames[i],x = operonMapper_files)] )
    RGIRes=read.csv(file = RGI_files[grep(GenomeNames[i],x = RGI_files)],header = T,sep = "\t")
    RGIRes$str_Crd=apply(RGIRes[,c("Start","Stop","Orientation")],1,function(x){str_c(str_trim(x,side = "both"),collapse = "_")})
    PULRes=tryCatch(read.csv(file = PULs_files[grep(GenomeNames[i],x = PULs_files)],header = T,sep = "\t"),error=function(x){NA})
    
    original_df=as.data.frame(genes(gbkFile))
    
    ### Format data to make GmakeRanges object:
    ##### 1. dbCAN2:
    
    if(!is.na(dbCAN2Res))
    {
      
      unqdbCAN2_GeneID=str_c(unique(dbCAN2Res$str_Crd))
      dbCAN2Sht=dbCAN2Res[match(unqdbCAN2_GeneID,as.character(dbCAN2Res$str_Crd)),]
      df=data.frame(start=dbCAN2Sht$Cluster_Start, end=dbCAN2Sht$Cluster_End,strand=dbCAN2Sht$Strand,
                    locus_tag=dbCAN2Sht$str_Crd,gene_id=dbCAN2Sht$str_Crd)
      
      df$seqnames=apply(dbCAN2Sht[,c("str_Crd","dbCAN2_annotation")],1,function(x){str_c("Gene-ID=",x[1],";","annotation=",x[2],";")})
      df$type="gene"
      df$loctype="normal"
      df$pseudo=""
      df$gene=""
      df$gene_synonym=""
      df$width=abs(df$start-df$end)-1
      
      dfFin=rbind(original_df,df[,colnames(original_df)])
      
    }
    
    ##### 2. kofamscan:
    
    if(!is.na(kofamscanRes))
    {
      
      unqkofamscan_GeneID=str_c(unique(kofamscanRes$str_Crd))
      kofamscanSht=kofamscanRes[match(unqkofamscan_GeneID,as.character(kofamscanRes$str_Crd)),]
      kofamscanSht[c('start', 'end','strand')] <- str_split_fixed(kofamscanSht$str_Crd, '_', 3)
      kofamdf=data.frame(start=kofamscanSht$start, end=kofamscanSht$end,strand=kofamscanSht$strand,
                         locus_tag=kofamscanSht$str_Crd,gene_id=kofamscanSht$str_Crd)
      kofamdf$seqnames=apply(kofamscanSht[,c("str_Crd","KO","KO_definition")],1,function(x){str_c("Gene-ID=",x[1],";","KO=",x[2],";","KO-Description=",x[3],";")})
      kofamdf$type="gene"
      kofamdf$loctype="normal"
      kofamdf$pseudo=""
      kofamdf$gene=""
      kofamdf$gene_synonym=""
      kofamdf$width=abs(as.numeric(as.character(kofamdf$start))-as.numeric(as.character(kofamdf$end)))-1
      dfFin=rbind(dfFin,kofamdf[,colnames(original_df)])
      
    }
    
    ##### 3. gff:
    
    if(!is.na(gffRes))
    {
      unqgff_GeneID=str_c(unique(gffRes$str_Crd))
      gffSht=gffRes[match(unqgff_GeneID,as.character(gffRes$str_Crd)),]
      gffSht[c('start', 'end','strand')] <- str_split_fixed(gffSht$str_Crd, '_', 3)
      gffdf=data.frame(start=gffSht$start, end=gffSht$end,strand=gffSht$strand,
                       locus_tag=gffSht$str_Crd,gene_id=gffSht$str_Crd)
      gffdf$seqnames=apply(gffSht[,c("str_Crd","locus_tag","Product","gene")],1,function(x){str_c("Gene-ID=",x[1],";","locus_tag=",x[2],";","Product=",x[3],";","Protein-ID=",x[4],";")})
      gffdf$type="gene"
      gffdf$loctype="normal"
      gffdf$pseudo=""
      gffdf$gene=""
      gffdf$gene_synonym=""
      gffdf$width=abs(as.numeric(as.character(gffdf$start))-as.numeric(as.character(gffdf$end)))-1
      dfFin=rbind(dfFin,gffdf[,colnames(original_df)])
      
      
    }
    
    
    ##### 4. eggNOG:
    
    if(!is.na(eggNOGDBRes))
    {
      
      unqeggNOGDB_GeneID=str_c(unique(eggNOGDBRes$str_Crd))
      eggNOGDBSht=eggNOGDBRes[match(unqeggNOGDB_GeneID,as.character(eggNOGDBRes$str_Crd)),]
      eggNOGDBSht[c('start', 'end','strand')] <- str_split_fixed(eggNOGDBSht$str_Crd, '_', 3)
      eggNOGdf=data.frame(start=eggNOGDBSht$start, end=eggNOGDBSht$end,strand=eggNOGDBSht$strand,
                          seqid=eggNOGDBSht$str_Crd)
      
      eggNOGdf$group=apply(eggNOGDBSht[,c("str_Crd","Description","COG_category")],1,function(x){str_c("Gene-ID=",x[1],";","Description=",x[2],";","COG_caterogy=",x[3],";")})
      
      eggNOGdf$type="gene"
      eggNOGdf$score=NA
      eggNOGdf$phase=NA
      eggNOGdf$source="Tool:eggNOG2.1.9,DB:eggNOGv5"
      dfFin=rbind(dfFin,eggNOGdf[,colnames(original_df)])
      
    }
    
    
    ##### 5. PULs:
    
    if(!is.na(PULRes))
    {
      PULRes$str_Crd=apply(PULRes[,c("start","end","strand")],1,function(x){str_c(x,collapse = "_")})
      unqPULRes_GeneID=str_c(unique(PULRes$str_Crd))
      PULResSht=PULRes[match(unqPULRes_GeneID,as.character(PULRes$str_Crd)),]
      PULdf=data.frame(start=PULResSht$start, end=PULResSht$end,strand=PULResSht$strand,
                       locus_tag=PULResSht$str_Crd,gene_id=PULResSht$str_Crd)
      PULdf$seqnames=rep(unique(original_df$seqnames)[1],nrow(PULdf))
      PULdf$type="gene"
      PULdf$loctype="normal"
      PULdf$pseudo=""
      PULdf$gene=""
      PULdf$gene_synonym=""
      PULdf$width=abs(as.numeric(as.character(PULdf$start))-as.numeric(as.character(PULdf$end)))-1
      dfFin=rbind(dfFin,PULdf[,colnames(original_df)])
      
      
    }
    
    
    ##### 6. RGI:
    if(!is.na(RGIRes))
    {
      
      unqRGIRes_GeneID=str_c(unique(RGIRes$str_Crd))
      RGIResSht=RGIRes[match(unqRGIRes_GeneID,as.character(RGIRes$str_Crd)),]
      RGIResSht[c('start', 'end','strand')] <- str_split_fixed(RGIResSht$str_Crd, '_', 3)
      RGIdf=data.frame(start=RGIResSht$start, end=RGIResSht$end,strand=RGIResSht$strand,
                       locus_tag=RGIResSht$str_Crd,gene_id=RGIResSht$str_Crd)
      
      RGIdf$seqnames=rep(unique(original_df$seqnames)[1],nrow(RGIdf))
      RGIdf$type="gene"
      RGIdf$loctype="normal"
      RGIdf$pseudo=""
      RGIdf$gene=""
      RGIdf$gene_synonym=""
      RGIdf$width=abs(as.numeric(as.character(RGIdf$start))-as.numeric(as.character(RGIdf$end)))-1
      dfFin=rbind(dfFin,RGIdf[,colnames(original_df)])
      
    }
    
    
    ##### 7. novelRes:
    
    if(!is.na(novelRes))
    {
      
      unqnovelRes_GeneID=str_c(unique(novelRes$str_Crd))
      novelResSht=novelRes[match(unqnovelRes_GeneID,as.character(novelRes$str_Crd)),]
      novelResSht[c('start', 'end','strand')] <- str_split_fixed(novelResSht$str_Crd, '_', 3)
      novelResdf=data.frame(start=novelResSht$start, end=novelResSht$end,strand=novelResSht$strand,
                            locus_tag=novelResSht$str_Crd,gene_id=novelResSht$str_Crd)
      novelResdf$seqnames=rep(unique(original_df$seqnames)[1],nrow(novelResdf))
      novelResdf$type="gene"
      novelResdf$loctype="normal"
      novelResdf$pseudo=""
      novelResdf$gene=""
      novelResdf$gene_synonym=""
      novelResdf$width=abs(as.numeric(as.character(novelResdf$start))-as.numeric(as.character(novelResdf$end)))-1
      dfFin=rbind(dfFin,novelResdf[,colnames(original_df)])
      
    }
    
    
    ##### 8. operon:
    if(!is.na(operonRes))
    {
      
      unqoperon_GeneID=str_c(unique(operonRes$str_Crd))
      operonSht=operonRes[match(unqoperon_GeneID,as.character(operonRes$str_Crd)),]
      operonSht[c('start', 'end','strand')] <- str_split_fixed(operonSht$str_Crd, '_', 3)
      operondf=data.frame(start=operonSht$start, end=operonSht$end,strand=operonSht$strand,
                          locus_tag=operonSht$str_Crd,gene_id=operonSht$str_Crd)
      operondf$seqnames=rep(unique(original_df$seqnames)[1],nrow(operondf))
      operondf$type="gene"
      operondf$loctype="normal"
      operondf$pseudo=""
      operondf$gene=""
      operondf$gene_synonym=""
      operondf$width=abs(as.numeric(as.character(operondf$start))-as.numeric(as.character(operondf$end)))-1
      dfFin=rbind(dfFin,operondf[,colnames(original_df)])
      
    }
    
    rfsqGenes=tryCatch(expr = as.data.frame(genes(rfsqFile)),error=function(x){NA})
    
    if(!is.na(rfsqGenes))
    {
      dfFin=rbind(df,rfsqGenes[,colnames(original_df)])
      
      
    }
    
    dfFin=dfFin [!duplicated(dfFin[c(2,3,5)]),]
    
    GRanges_obj=makeGRangesFromDataFrame(dfFin,
                                         keep.extra.columns=FALSE,
                                         ignore.strand=FALSE,
                                         seqinfo=NULL,
                                         seqnames.field=c("seqnames", "seqname",
                                                          "chromosome", "chrom",
                                                          "chr", "chromosome_name",
                                                          "seqid"),
                                         start.field="start",
                                         end.field=c("end", "stop"),
                                         strand.field="strand",
                                         starts.in.df.are.0based=FALSE)
    
    gbkFile@genes=GRanges_obj
    gbkFile@cds=GRanges_obj
    gbkFile@exons=GRanges_obj
    gbkFile@transcripts=GRanges_obj
    
    GRanges_obj=makeGRangesFromDataFrame(dfFin,
                                         keep.extra.columns=TRUE,
                                         ignore.strand=FALSE,
                                         seqinfo=NULL,
                                         seqnames.field=c("seqnames", "seqname",
                                                          "chromosome", "chrom",
                                                          "chr", "chromosome_name",
                                                          "seqid"),
                                         start.field="start",
                                         end.field=c("end", "stop"),
                                         strand.field="strand",
                                         starts.in.df.are.0based=FALSE)
    
    if(!is.na(gbkFile) && !is.na(rfsqFile))
    {
      
      
    } else {
      
      if(is.na(gbkFile))
      {
        overallGb=rfsqFile
      } else {
        
        overallGb=gbkFile
      }
      
      
    }
    
    
    ############################################################################
    # Add annotations to merged gbk file:---------------------------------------
    ############################################################################
    if(i!=13)
    {
      AnnotationFl=read.csv(file = str_c(AnnFolder,gnmFile$Genome),"_allAnnJoint.csv")  
      
    }
  }
  
  
  
}