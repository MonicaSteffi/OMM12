annotation_main=function(GenomeNamesHere,dbCAN2_files,eggNOGDB_crd_ko_files,NovelFam_crd_files,gff_files,kofamscan_files,operonMapper_files,PULs_files,RGI_files,kofam_for_eggnog_files,kofam_for_operon_files,eggNOG_for_operon_files,novFam_for_operon_files,eggNOG_for_integrated_proteins,gnmFile)
{
  GenomeID_columns=c("GenomeID","Organism","eggNOG_query","novelfam_query","Operon_GenomeID","Contig")
  setClass(Class = "non_redundant_genes",representation(to_be_keptDf="data.frame",
                                                        UC_stpOverall="data.frame",
                                                        UC_startOverall="data.frame") )
  
  source("S:/Abilash/Marbouty_annotation_refseq_genbank_2022/scripts/R_files/coordinate_based_non_redundant_genes.R")
  
  dbCAN2Res=readRDS(file = dbCAN2_files[grep(GenomeNamesHere,x = dbCAN2_files)])
  dbCAN2Res$str_Crd=str_replace_all(string = dbCAN2Res$str_Crd,pattern = " ",replacement = "")
  
  kofamscanRes=readRDS(file = kofamscan_files[grep(GenomeNamesHere,x = kofamscan_files)])
  gffRes=readRDS(file = gff_files[grep(GenomeNamesHere,x = gff_files)])
  gffRes$locus_tag=str_replace_all(string = gffRes$locus_tag,pattern = "_Non_overlapping",replacement = "")
  gffRes$locus_tag=str_replace_all(string = gffRes$locus_tag,pattern = "_Overlapping",replacement = "")
  gffRes=gffRes[(!duplicated(gffRes)),]
  unqLocus=unique(gffRes$locus_tag)
  
  
  
  kofamscanRes$str_Crd=gffRes$str_Crd[match(kofamscanRes$GeneID,as.character(gffRes$locus_tag))]
  na_kofamscanRes=which(is.na(kofamscanRes$str_Crd))
  if(length(na_kofamscanRes) > 0)
  {
    for(na_kfamCnt in 1: length(na_kofamscanRes))
    {
      kofam_na_idx_rfsq=which(as.character(gffRes$locus_tag_refseq)==kofamscanRes$GeneID[na_kofamscanRes[na_kfamCnt]])
      if(length(kofam_na_idx_rfsq) > 0)
      {
        kofamscanRes$str_Crd[na_kofamscanRes[na_kfamCnt]]=unique(gffRes$str_Crd[kofam_na_idx_rfsq])
        
      }
      
    }
    
  }
  
  na_kofamscanRes=which(is.na(kofamscanRes$str_Crd))
  if(length(na_kofamscanRes) > 0)
  {
    for(na_kfamCnt in 1: length(na_kofamscanRes))
    {
      kofam_na_idx=which(as.character(gffRes$locus_tag_genbank)==kofamscanRes$GeneID[na_kofamscanRes[na_kfamCnt]])
      if(length(kofam_na_idx) > 0)
      {
        kofamscanRes$str_Crd[na_kofamscanRes[na_kfamCnt]]=unique(gffRes$str_Crd[kofam_na_idx])
        
      }
        
    }
    
  }
  
  #### There is a chance that eggNOG and novelfam can predict unique proteins, hence their kofam annotations should be in alignment with other annotations as well.
  eggNOGDBRes=readRDS(file = eggNOGDB_crd_ko_files[grep(GenomeNamesHere,x = eggNOGDB_crd_ko_files)])
  # colnames(eggNOGDBRes)[which(colnames(eggNOGDBRes)=="evalue")]="eggNOG_evalue"
  # colnames(eggNOGDBRes)[which(colnames(eggNOGDBRes)=="score")]="eggNOG_score"
  # colnames(eggNOGDBRes)[which(colnames(eggNOGDBRes)=="query")]="eggNOG_query"
  # 
  # 
  # 
  # kofam4eggNOG=readRDS(file = kofam_for_eggnog_files[grep(GenomeNamesHere,x = kofam_for_eggnog_files)])
  # colnames(kofam4eggNOG)=str_c("ko4eggNOG_",colnames(kofam4eggNOG))
  # colnames(kofam4eggNOG)[2]="eggNOG_query"
  # 
  # 
  # 
  # eggNOGDBRes=merge(eggNOGDBRes,kofam4eggNOG,by="eggNOG_query",all=TRUE)
  # eggNOGDBRes=eggNOGDBRes[,-grep(pattern = "KEGG_|BRITE",x=colnames(eggNOGDBRes))]
  # 
  novelRes=tryCatch(readRDS(file = NovelFam_crd_files[grep(GenomeNamesHere,x = NovelFam_crd_files)]),error=function(x){NA})
  
  if(!is.na(novelRes))
  {
    nov_famIdx=grep(pattern = "novel_fam",x=colnames(novelRes))
    if(length(nov_famIdx) > 0)
    {
      novelRes=novelRes[,-nov_famIdx]
      
    }
    colnames(novelRes)[which(colnames(novelRes)=="target")]="novelfam_target"
    colnames(novelRes)[which(colnames(novelRes)=="score")]="novelfam_score"
    colnames(novelRes)[which(colnames(novelRes)=="KO")]="novelfam_KO"
    colnames(novelRes)[which(colnames(novelRes)=="KO_Desc")]="novelfam_KO_Desc"
    colnames(novelRes)[which(colnames(novelRes)=="query")]="novelfam_query"
    
    
    eggNOGDBRes=merge(eggNOGDBRes,novelRes,by="str_Crd",all=TRUE)
    
  }
  
  eggNOG_for_operon=read.csv(file = eggNOG_for_operon_files[grep(pattern = GenomeNamesHere,x = eggNOG_for_operon_files)],sep = "\t",header = F,comment.char = "#")
  colnames(eggNOG_for_operon)=c("query","seed_ortholog","evalue","score","eggNOG_OGs","max_annot_lvl","COG_category","Description",
                                "Preferred_name",	"GOs",	"EC",	"KEGG_ko",	"KEGG_Pathway","KEGG_Module","KEGG_Reaction",	"KEGG_rclass",
                                "BRITE","KEGG_TC",	"CAZy",	"BiGG_Reaction",	"PFAMs")
  
  novFamGrep=novFam_for_operon_files[grep(pattern = GenomeNamesHere,x = novFam_for_operon_files)]
  if(length(novFamGrep)==1 )
  {
    novFam_for_operon=tryCatch(read.csv(file = novFamGrep,sep = "\t",header = F,comment.char = "#"),error=function(x){NA})
    
    if(!is.na(novFam_for_operon))
    {
      colnames(novFam_for_operon)=c("query","novelFam_target","novelFam_evalue","novelFam_score","novelFamID")
      eggNOG_for_operon=merge(eggNOG_for_operon,novFam_for_operon,by="query",all=TRUE)  
      
    }
  }
  
  operonRes=readRDS(file =operonMapper_files[grep(GenomeNamesHere,x = operonMapper_files)] )
  operonRes$query=str_c(GenomeNamesHere,";",operonRes$GeneID)
  
  operonRes=merge(x = operonRes,y=eggNOG_for_operon,by="query",all=TRUE)
  rm_op_cols=grep(pattern = "KEGG_|BRITE|COGgene|Function",x=colnames(operonRes))
  
  if(length(rm_op_cols) > 0)
  {
    operonRes=operonRes[,-(rm_op_cols)]
    
  }
  operonRes$COG4Operon=operonRes$eggNOG_OGs
  operonRes$COG4Operon=str_replace_all(string = operonRes$COG4Operon,pattern = "@.*",replacement = "")
  
  kofam4operon=readRDS(file = kofam_for_operon_files[grep(GenomeNamesHere,x = kofam_for_operon_files)])
  colnames(kofam4operon)=str_c("ko4Operon_",colnames(kofam4operon))
  colnames(kofam4operon)[2]="query"
  operonRes=merge(x = operonRes,y = kofam4operon,by="query",all=TRUE)
  
  
  RGIRes=read.csv(file = RGI_files[grep(GenomeNamesHere,x = RGI_files)],header = T,sep = "\t")
  RGIRes$str_Crd=apply(RGIRes[,c("Start","Stop","Orientation")],1,function(x){str_c(str_trim(x,side = "both"),collapse = "_")})
  PULRes=tryCatch(read.csv(file = PULs_files[grep(GenomeNamesHere,x = PULs_files)],header = T,sep = "\t"),error=function(x){NA})
  
  
  ### Integrate files step-by-step:
  #### 1. Integrate dbCAN2 with gff:
  dC_gff=merge(x = dbCAN2Res,y = gffRes,by="str_Crd",all=TRUE)
  #### 2. Integrate dbCAN2,gff with eggNOGDB:
  dC_gff_eN=merge(x = dC_gff,y = eggNOGDBRes,by="str_Crd",all=TRUE)
  #### 3. Integrate dbCAN2,gff,eggNOGDB with OperonMapper:
  dC_gff_eN_O=merge(x = dC_gff_eN,y = operonRes,by="str_Crd",all=TRUE)
  
  #### 4. Integrate dbCAN2,gff,eggNOGDB,OperonMapper with NovelFamDB:
  # if(!is.null(dim(novelRes)))
  # {
  #   dC_gff_eN_O_N=merge(x = dC_gff_eN_O,y = novelRes,by="str_Crd",all=TRUE)
  #   
  # } else {
  #   dC_gff_eN_O_N=dC_gff_eN_O
  #   
  # }
  # 
  
  dC_gff_eN_O_N=dC_gff_eN_O
  if(!is.null(dim(PULRes)))
  {
    PULRes$str_Crd=apply(PULRes[,c("start","end","strand")],1,function(x){str_c(str_trim(x,side="both"),collapse = "_")})
    dC_gff_eN_O_N_P=merge(x = dC_gff_eN_O_N,y=PULRes,by="str_Crd",all=TRUE)
    
  } else {
    
    dC_gff_eN_O_N_P=dC_gff_eN_O_N
  }
  
  if(!is.null(dim(RGIRes)))
  {
    
    dC_gff_eN_O_N_P_R=merge(x = dC_gff_eN_O_N_P,y = RGIRes,by="str_Crd",all=TRUE)
  } else {
    dC_gff_eN_O_N_P_R=dC_gff_eN_O_N_P
    
  }
  
  dC_gff_eN_O_N_P_R_KO=merge(x = dC_gff_eN_O_N_P_R,y = kofamscanRes,by="str_Crd",all=TRUE)
  NA_str_crd=unique(c(which(is.na(dC_gff_eN_O_N_P_R_KO$str_Crd),which(dC_gff_eN_O_N_P_R_KO$str_Crd==""))))
  if(length(NA_str_crd) > 0) 
  {
    
  tmpStartVec=vector()
   for(stCrd in 1: length(NA_str_crd))
   {
     tmpStart=dC_gff_eN_O_N_P_R_KO[NA_str_crd[stCrd],c("Cluster_Start","Start.x","PosLeft","Start.y")]
     tmpStart=tmpStart[!is.na(tmpStart)];
     tmpStart=tmpStart[tmpStart!=""];
     tmpStartVec[stCrd]=unique(tmpStart);
     
   }
    
  
    
 
  
  
  tmpStpVec=vector()
    for(stCrd in 1: length(NA_str_crd))
    {
      tmpStp=dC_gff_eN_O_N_P_R_KO[NA_str_crd[stCrd],c("Cluster_End","Stop.x","postRight","Start.y")]
      tmpStp=tmpStp[!is.na(tmpStp)];
      tmpStp=tmpStp[tmpStp!=""];
      tmpStpVec[stCrd]=unique(tmpStp);
      
    }
    
  
  
  tmpStrndVec=vector()
    for(stCrd in 1: length(NA_str_crd))
    {
      tmpStrnd=dC_gff_eN_O_N_P_R_KO[NA_str_crd[stCrd],c("Strand.x", "Strand.y" ,"Strand" )]
      tmpStrnd=tmpStrnd[!is.na(tmpStrnd)];
      tmpStrnd=tmpStrnd[tmpStrnd!=""];
      tmpStrndVec[stCrd]=unique(tmpStrnd);
      
    }
    
  
  
  str_crd_naVal=apply(cbind(tmpStartVec,tmpStpVec,tmpStrndVec),1,function(x){
    
    str_c(x,collapse = "_")
  })
  dC_gff_eN_O_N_P_R_KO$str_Crd[NA_str_crd]=""
  dC_gff_eN_O_N_P_R_KO$str_Crd[NA_str_crd]=as.character(str_crd_naVal)
  dC_gff_eN_O_N_P_R_KO[c('GeneStart', 'GeneEnd','GeneStrand')]=str_split_fixed(dC_gff_eN_O_N_P_R_KO$str_Crd, '_', 3)
  # dC_gff_eN_O_N_P_R_KO$GeneEnd=as.numeric(dC_gff_eN_O_N_P_R_KO$GeneEnd)
  # dC_gff_eN_O_N_P_R_KO$GeneStart=as.numeric(dC_gff_eN_O_N_P_R_KO$GeneStart)
  
  } else {
    
    dC_gff_eN_O_N_P_R_KO[c('GeneStart', 'GeneEnd','GeneStrand')]=str_split_fixed(dC_gff_eN_O_N_P_R_KO$str_Crd, '_', 3)
    
  }
  rm_overlappingAnn=grep(pattern = "verlapping",x=dC_gff_eN_O_N_P_R_KO$locus_tag)
  
  if(length(rm_overlappingAnn) > 0)
  {
    dC_gff_eN_O_N_P_R_KO=dC_gff_eN_O_N_P_R_KO[-(rm_overlappingAnn),]
    
  }
  
  
  
  
  
  dC_gff_eN_O_N_P_R_KO$COG4gbk_rfsq=eggNOG_for_integrated_proteins$COGgene[match(dC_gff_eN_O_N_P_R_KO$locus_tag,as.character(eggNOG_for_integrated_proteins$locus_tag))]
  dC_gff_eN_O_N_P_R_KO$COG4eggNOG=str_replace_all(string = dC_gff_eN_O_N_P_R_KO$eggNOG_OGs.x,pattern = "@.*",replacement = "")
  dC_gff_eN_O_N_P_R_KO$GenomeID=str_replace_all(string = dC_gff_eN_O_N_P_R_KO$GenomeID,pattern = "NZ_",replacement = "")
  dC_gff_eN_O_N_P_R_KO$GenomeID=apply(dC_gff_eN_O_N_P_R_KO[,intersect(GenomeID_columns,colnames(dC_gff_eN_O_N_P_R_KO))],1,function(x){
    x=x[!is.na(x)];
    x=x[x!=""];
    x=str_replace_all(string = x,pattern = "_.*",replacement = "");
    x=str_replace_all(string = x,pattern = "_.*",replacement = "");
    x=unique(x);
    if(length(x)==1)
    {
      return(x)
    } else {
      return(NA)
    }
    
  })
  
  unqGenomeID=unique(dC_gff_eN_O_N_P_R_KO$GenomeID)
  unqGenomeID=unqGenomeID[!is.na(unqGenomeID)]
  unqGenomeID=unqGenomeID[unqGenomeID!=""]
  print(str_c("########################################################################################################"))
  print(str_c("Genome-Name : ", unique(gnmFile$Genome[match(unqGenomeID,as.character(gnmFile$GenbankID))])))
  
  print(str_c("Genome-IDs : ", str_c(unqGenomeID,collapse = ",")))
  print(str_c("Total non-redundant genes in dC_gff_eN_O_N_P_R_KO: ", nrow(dC_gff_eN_O_N_P_R_KO)))
  print(str_c("########################################################################################################"))
  
  genomeSizeVec=vector(length = length(unqGenomeID))
  lastGeneVec=vector(length = length(unqGenomeID))
  
  for(gbk1 in 1: length(unqGenomeID))
  {
    dC_gff_eN_O_N_P_R_KO_Sht=dC_gff_eN_O_N_P_R_KO[which(dC_gff_eN_O_N_P_R_KO$GenomeID==unqGenomeID[gbk1]),]
    genomeSize=as.numeric(as.character(gnmFile$Genome_size[which(gnmFile$GenbankID==unqGenomeID[gbk1])]))
    # lastGene=as.character((dC_gff_eN_O_N_P_R_KO_Sht$GeneEnd[!is.na(dC_gff_eN_O_N_P_R_KO_Sht$GeneEnd)]))
    lastGene=max(as.numeric(as.character((dC_gff_eN_O_N_P_R_KO_Sht$GeneEnd[!is.na(dC_gff_eN_O_N_P_R_KO_Sht$GeneEnd)]))))
    tmpGenomeName=as.character(unique(gnmFile$Genome[which(gnmFile$GenbankID==unqGenomeID[gbk1])]))
    
    
    genomeSizeVec[gbk1]=as.numeric(as.character(genomeSize))
    lastGeneVec[gbk1]=lastGene
    
    
    if(lastGene > genomeSize)
    {
      newLine=dC_gff_eN_O_N_P_R_KO_Sht[(grep(pattern = as.character(lastGene),x= dC_gff_eN_O_N_P_R_KO_Sht$GeneEnd)),]  
      tmpStrnd=unname(apply(newLine[,grep(pattern = "Strand",x=colnames(newLine))],1,function(x){
        x=x[x!=""];
        x=x[!is.na(x)];
        unique(x)
      }))
      newLine$str_Crd=str_c("1_",(lastGene-genomeSize),"_",tmpStrnd)
      newLine$GeneEnd=(lastGene-genomeSize)
      newLine$GeneStart=1
      
      newLine$PosLeft=newLine$GeneStart
      newLine$postRight=newLine$GeneEnd
      newLine$Start.x=newLine$GeneStart
      newLine$Stop.x=newLine$GeneEnd
      newLine$Start.y=newLine$GeneStart
      newLine$Stop.y=newLine$GeneEnd
      newLine$locus_tag=str_c(newLine$locus_tag,"_Edge")
      
      
      dC_gff_eN_O_N_P_R_KO_Sht$GeneEnd=str_replace_all(string = dC_gff_eN_O_N_P_R_KO_Sht$GeneEnd,pattern = as.character(lastGene),replacement = as.character(genomeSize))
      dC_gff_eN_O_N_P_R_KO_Sht$str_Crd=str_replace_all(string = dC_gff_eN_O_N_P_R_KO_Sht$str_Crd,pattern = as.character(lastGene),replacement = as.character(genomeSize))
      dC_gff_eN_O_N_P_R_KO_Sht=rbind(dC_gff_eN_O_N_P_R_KO_Sht,newLine)
      
      int_gnmsize_lastGene=intersect(which(dC_gff_eN_O_N_P_R_KO_Sht$GeneStart=="1"),
                                     which(dC_gff_eN_O_N_P_R_KO_Sht$GeneEnd==genomeSize))
      
      if(length(int_gnmsize_lastGene) > 0)
      {
        dC_gff_eN_O_N_P_R_KO_Sht=dC_gff_eN_O_N_P_R_KO_Sht[-(int_gnmsize_lastGene),]
        
      }
      
      
      
    } 
    
   
    
    if(gbk1==1)
    {
      new_dC_gff_eN_O_N_P_R_KO=dC_gff_eN_O_N_P_R_KO_Sht
     
    } else {
      
      new_dC_gff_eN_O_N_P_R_KO=rbind(new_dC_gff_eN_O_N_P_R_KO,dC_gff_eN_O_N_P_R_KO_Sht)
    }
    
    newAnn=dC_gff_eN_O_N_P_R_KO[intersect(which(dC_gff_eN_O_N_P_R_KO$GeneStart=="1"),
                                          match(genomeSizeVec,dC_gff_eN_O_N_P_R_KO$GeneEnd)),]
    
    print(str_c("########################################################################################################"))
    print(str_c("Genome-Name : ", unique(gnmFile$Genome[match(unqGenomeID[gbk1],as.character(gnmFile$GenbankID))])))
    
    print(str_c("Genome-IDs : ", str_c(unqGenomeID[gbk1],collapse = ",")))
    print(str_c("Total non-redundant genes in dC_gff_eN_O_N_P_R_KO_Sht: ", nrow(dC_gff_eN_O_N_P_R_KO_Sht)))
    print(str_c("########################################################################################################"))
    
   
    print(str_c("########################################################################################################"))
    print(str_c("Genome-Name : ", unique(gnmFile$Genome[match(unqGenomeID[gbk1],as.character(gnmFile$GenbankID))])))
    
    print(str_c("Genome-IDs : ", str_c(unqGenomeID[gbk1],collapse = ",")))
    print(str_c("Total Genome annotation lines: ", nrow(newAnn)))
    print(str_c("########################################################################################################"))
    
    
    
  }
  
  dC_gff_eN_O_N_P_R_KO=new_dC_gff_eN_O_N_P_R_KO
  
  gnmCoord=intersect(which(dC_gff_eN_O_N_P_R_KO$GeneStart=="1"),
                     unlist(lapply(genomeSizeVec,function(x)
                       {
                       which(dC_gff_eN_O_N_P_R_KO$GeneEnd==x)
                     })))
  
  if(length(gnmCoord)>0)
  {
    dC_gff_eN_O_N_P_R_KO=dC_gff_eN_O_N_P_R_KO[-(gnmCoord),]
    
  }
  
  print(str_c("########################################################################################################"))
  print(str_c("Genome-Name : ", unique(gnmFile$Genome[match(unqGenomeID[gbk1],as.character(gnmFile$GenbankID))])))
  
  print(str_c("Genome-IDs : ", str_c(unqGenomeID[gbk1],collapse = ",")))
  print(str_c("Total non-redundant genes in dC_gff_eN_O_N_P_R_KO Final: ", nrow(dC_gff_eN_O_N_P_R_KO)))
  print(str_c("########################################################################################################"))
  
  
  print(str_c("########################################################################################################"))
  print(str_c("Genome-Name : ", unique(gnmFile$Genome[match(unqGenomeID[gbk1],as.character(gnmFile$GenbankID))])))
  
  print(str_c("Genome-IDs : ", str_c(unqGenomeID[gbk1],collapse = ",")))
  print(str_c("Total Genome annotation lines: ", nrow(newAnn)))
  print(str_c("########################################################################################################"))
  
  
  saveRDS(dC_gff_eN_O_N_P_R_KO,str_c("S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/",GenomeNamesHere,"_allAnnotationsCompleteRaw.rds"))
  
  eggNOG_for_integrated_proteins=eggNOG_for_integrated_proteins[which(eggNOG_for_integrated_proteins$Genome==GenomeNamesHere),]
  
  
  for(st_stp_cnter in 1: length(unqGenomeID))
  {
    dC_gff_eN_O_N_P_R_KO_Sht=dC_gff_eN_O_N_P_R_KO[which(dC_gff_eN_O_N_P_R_KO$GenomeID==unqGenomeID[st_stp_cnter]),]
    
    non_rd_genesObj=coordinate_based_non_redundant_genes(dC_gff_eN_O_N_P_R_KO_Sht = dC_gff_eN_O_N_P_R_KO_Sht)
    
    to_be_keptDf=non_rd_genesObj@to_be_keptDf
    UC_stpOverall=non_rd_genesObj@UC_stpOverall
    UC_startOverall=non_rd_genesObj@UC_startOverall
    
    dim_obj=c(dim(to_be_keptDf)[1],dim(UC_stpOverall)[1],dim(UC_startOverall)[1])
    names(dim_obj)=c("to_be_keptDf","UC_stpOverall","UC_startOverall")
    dim_obj=dim_obj[which(dim_obj>0)]
    
    
    if(length(dim_obj)>0)
    {
      for(c1 in 1: length(dim_obj))
      {
        tmpObj1=get(x = names(dim_obj)[c1])
        if(c1 ==1)
        {
          
          FinDf=tmpObj1
          
        } else {
          FinDf=rbind(FinDf,tmpObj1)
          
        }
        
        
      }
      
    }
    
    
    non_rd_genesObjFin=coordinate_based_non_redundant_genes(dC_gff_eN_O_N_P_R_KO_Sht = FinDf)
    
    to_be_keptDf=non_rd_genesObjFin@to_be_keptDf
    UC_stpOverall=non_rd_genesObjFin@UC_stpOverall
    UC_startOverall=non_rd_genesObjFin@UC_startOverall
    
    dim_obj=c(dim(to_be_keptDf)[1],dim(UC_stpOverall)[1],dim(UC_startOverall)[1])
    names(dim_obj)=c("to_be_keptDf","UC_stpOverall","UC_startOverall")
    dim_obj=dim_obj[which(dim_obj>0)]
    
    
    if(length(dim_obj)>0)
    {
      for(c1 in 1: length(dim_obj))
      {
        tmpObj1=get(x = names(dim_obj)[c1])
        if(c1 ==1)
        {
          
          FinDf=tmpObj1
          
        } else {
          FinDf=rbind(FinDf,tmpObj1)
          
        }
        
        
      }
      
    }
    FinDf=FinDf[!duplicated(FinDf),]
    
    
    
    sDff_Columns=setdiff(colnames(FinDf),colnames(newAnn))
    newAnnExtra=as.data.frame(matrix(rep("",length(sDff_Columns)),nrow = nrow(newAnn),ncol=length(sDff_Columns)))
    colnames(newAnnExtra)=sDff_Columns
    newAnn=cbind(newAnn,newAnnExtra)
    FinDf=rbind(newAnn,FinDf)
    unqStr_Crd=unique(FinDf$str_Crd)
   
    for(cll1 in 1: ncol(FinDf))
    {
      FinDf[,cll1]=as.character(FinDf[,cll1])
      
    }
    
    k3=1
    for(f1 in 1: length(unqStr_Crd))
    {
      tmpIdx=which(FinDf$str_Crd==unqStr_Crd[f1])
      
      if(length(tmpIdx)>0)
      {
        tmpMat=FinDf[tmpIdx,]
        
        if(length(tmpIdx) > 1)
        {
          
          tmpMatSht=as.character(apply(tmpMat,2,function(x){
            x=as.character(x);
            x1=unique(x);
            x1=x1[!is.na(x1)];
            x1=x1[x1!=""];
            return(str_c(unique(x1),collapse = "///"))
            
          }))
          
        } else {
          
          tmpMatSht=tmpMat
        }
        
        if(k3==1)
        {
          FinMatSht=tmpMatSht
          k3=k3+1
        } else {
          
          FinMatSht=rbind(FinMatSht,tmpMatSht)
        }
      }
      
    }
    
    if(st_stp_cnter==1)
    {
      new_FinDf=FinDf
      new_FinMatSht=FinMatSht
    } else {
      
      new_FinDf=rbind(new_FinDf,FinDf)
      new_FinMatSht=rbind(new_FinMatSht,FinMatSht)
      
    }
  }
    
    
  if(file.exists(str_c("S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/",GenomeNamesHere,"_FinDf.rds")))
  {
    file.remove(str_c("S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/",GenomeNamesHere,"_FinDf.rds"))
  }
  
  if(file.exists(str_c("S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/",GenomeNamesHere,"_FinMatSht.rds")))
  {
    file.remove(str_c("S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/",GenomeNamesHere,"_FinMatSht.rds"))
  }
  
  
  saveRDS(new_FinDf,str_c("S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/",GenomeNamesHere,"_new_FinDf.rds"))
  saveRDS(new_FinMatSht,str_c("S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/",GenomeNamesHere,"_new_FinMatSht.rds"))
  
  
  
}