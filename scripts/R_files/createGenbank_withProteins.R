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
  
  
  annExcl=list.files(path ="S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/",pattern = "_FinMatSht.rds",full.names = T)
  gnmFile=read.csv(file = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/GenomeID_to_Organism_accession.csv",sep = ";",header = T)
  gbkFolder="S:/Abilash/Marbouty_annotation_refseq_genbank_2022//Data/Downloaded_from_genbank_and_refseq/Marbouty_genbank_format_files/"
  AnnFolder="S:/Abilash/Marbouty_annotation_refseq_genbank_2022/gffFiles/For_RNAseq_count/"
  fin_colnames=c("seqnames" , "start" , "end" , "width" , "strand" , "type" , "locus_tag" , "loctype" ,"gene" ,"gene_synonym","gene_id","database","str_Crd")
  ## first, merge the gbk files from genbank and refseq.
  ## second, add the annotations to the merged gbk files.
  
  ############################################################################
  # Merging gbk files:--------------------------------------------------------
  ############################################################################
  gnmFile$Organism=str_replace_all(string = gnmFile$Organism,pattern = "Escheri",replacement = "escheri")
  gnmFile$Genome=str_replace_all(string = gnmFile$Genome,pattern = "Escheri",replacement = "escheri")
  
  gnmFile$Genome=str_replace_all(string = gnmFile$Genome,pattern = "Akker",replacement = "akker")
  gnmFile$Organism=str_replace_all(string = gnmFile$Organism,pattern = "Akker",replacement = "akker")
  
  gnmFile$Genome=str_replace_all(string = gnmFile$Genome,pattern = "muribaculum_intestinale_",replacement = "muribaculum_intestinales_")
  gnmFile$Organism=str_replace_all(string = gnmFile$Organism,pattern = "muribaculum_intestinale_",replacement = "muribaculum_intestinales_")
  
  
  # feature_table_files=list.files(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Downloaded_from_genbank_and_refseq/Marbouty_Genbank/",pattern = "_feature_table.txt$",full.names = T,recursive = T)
  
  
  ##### protein files.
  operonMapper_files=list.files(pattern = "predicted_protein_sequences",path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/OperonMapper_new/",full.names = T,recursive = T)
  eggNOGDB_files=list.files(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/eggNOG_2022_new/eggNOGDB/",pattern = "out.emapper.genepred.fasta",full.names = T,recursive = T)
  NovelFamDB_files=list.files(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/eggNOG_2022_new/NovelFamDB/",pattern = "out.emapper.genepred.fasta",full.names = T,recursive = T)
  dbCAN2_files=list.files(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/",pattern = "uniInput.txt",full.names = T,recursive = T)
  integrated_proteins_files=list.files("S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_genbank_and_refseq_proteins_integrated/",pattern = ".faa",full.names = T,recursive = T)
  
  
  

  
  GenomeDesc=as.character(gnmFile$Organism)
  GenomeNm=unique(gnmFile$Genome)
  for(i in 11: length(GenomeNm))
  {
    
    tmpExcl=readRDS(annExcl[grep(pattern = GenomeNm[i],x=annExcl)])
    
    
    tmpOperonMapping=operonMapper_files[grep(pattern = GenomeNm[i],x = operonMapper_files)]
    tmpeggNOG=eggNOGDB_files[grep(pattern = str_c(as.character(gnmFile$Genbank_accession[which(gnmFile$Genome==GenomeNm[i])]),collapse = "|"),x = eggNOGDB_files)]
    tmpNovFam=NovelFamDB_files[grep(pattern = str_c(as.character(gnmFile$Genbank_accession[which(gnmFile$Genome==GenomeNm[i])]),collapse = "|"),x = NovelFamDB_files)]
    tmpdbCAN2=dbCAN2_files[grep(pattern = GenomeNm[i],x = dbCAN2_files)]
    tmpIntegrated=integrated_proteins_files[grep(pattern = GenomeNm[i],x = integrated_proteins_files)]
    
    
    op_proteins=lapply(tmpOperonMapping,function(x){
      read.fasta(file = x,seqtype = "AA",whole.header = T,set.attributes = F,as.string = T)
    })
    op_proteinsOverall=do.call(c,op_proteins)
    
    
    eggNOG_proteins=lapply(tmpeggNOG,function(x){
      read.fasta(file = x,seqtype = "AA",whole.header = T,set.attributes = F,as.string = T)
    })
    eggNOG_proteinsOverall=do.call(c,eggNOG_proteins)
    
    novFam_proteins=lapply(tmpNovFam,function(x){
      read.fasta(file = x,seqtype = "AA",whole.header = T,set.attributes = F,as.string = T)
    })
    novFam_proteinsOverall=do.call(c,novFam_proteins)
    
    dbCAN2_proteins=lapply(tmpdbCAN2,function(x){
      read.fasta(file = x,seqtype = "AA",whole.header = T,set.attributes = F,as.string = T)
    })
    dbCAN2_proteinsOverall=do.call(c,dbCAN2_proteins)
    
    integrated_proteins=lapply(tmpIntegrated,function(x){
      read.fasta(file = x,seqtype = "AA",whole.header = T,set.attributes = F,as.string = T)
    })
    integrated_proteinsOverall=do.call(c,integrated_proteins)
    
    Allproteins=c(op_proteinsOverall,eggNOG_proteinsOverall,novFam_proteinsOverall,
                  dbCAN2_proteinsOverall,integrated_proteinsOverall)
    names(Allproteins)=str_replace_all(string = names(Allproteins),pattern = " .*",replacement = "")
    proteinsList=rep("",nrow(tmpExcl))
    nms_proteins=rep("",nrow(tmpExcl))
  str_crdVec=vector()
    for(pr1 in 1: nrow(tmpExcl))
    {
      
      x=tmpExcl[pr1,c("str_Crd","NCBIStart","eggNOGStart","operonStart","dbCAN2Start","RGIStart")];
      x=as.character(as.vector(unname(x))); 
      str_Crd_Idx=which(tmpExcl$str_Crd==x[1]);
      locus_tags_gbk=tmpExcl$locus_tag[str_Crd_Idx];
      locus_tags_eggNOG=str_replace_all(tmpExcl$query.x[str_Crd_Idx],pattern = " .*",replacement = "")
      locus_tags_operon=str_replace_all(tmpExcl$GeneID.y[str_Crd_Idx],pattern = " .*",replacement = "")
      locus_tags_dbCAN2=tmpExcl$dbCAN2_GeneID[str_Crd_Idx]
      str_crdVec[pr1]=tmpExcl$str_Crd[str_Crd_Idx]    
      
        
      if(locus_tags_gbk!="" && !is.na(locus_tags_gbk))
      {
        sqnces=as.character(Allproteins[locus_tags_gbk]);
        nms_proteins[pr1]=locus_tags_gbk;
        proteinsList[pr1]=sqnces   
        next;
        
      } else {
      
        
        if(locus_tags_eggNOG!="" && !is.na(locus_tags_eggNOG))
        {
          
          sqnces=as.character(Allproteins[locus_tags_eggNOG]);
          nms_proteins[pr1]=locus_tags_eggNOG;
          proteinsList[pr1]=sqnces 
          next;
          
        } else {
          
          if(locus_tags_operon!="" && !is.na(locus_tags_operon))
          {
            
            sqnces=as.character(Allproteins[locus_tags_operon]);
            nms_proteins[pr1]=locus_tags_operon;
            proteinsList[pr1]=sqnces;
            next;
            
          } else {
            if(locus_tags_dbCAN2!="" && !is.na(locus_tags_dbCAN2))
            {
              
              sqnces=as.character(Allproteins[locus_tags_dbCAN2]);
              nms_proteins[pr1]=locus_tags_dbCAN2;
              proteinsList[pr1]=sqnces 
              next;
              
            } else {
              if(x[6]!="" && !is.na(x[6]))
              {
                
                sqnces=as.character(tmpExcl$Predicted_Protein[which(tmpExcl$str_Crd==x[1])]);
                nms_proteins[pr1]=locus_tags;
                proteinsList[pr1]=sqnces
                next;  
                
              }
              
            }
            
            
          }
          
          
        }
        
        
      }
      
      
      
      
      
    }
    names(proteinsList)=nms_proteins
    
    tmpExcl$ProteinSequences=proteinsList[match(tmpExcl$str_Crd,str_crdVec)]
     
    ### Read gbk files:
    
    gbkFile=tryCatch(readGenBank(file = str_c(gbkFolder,as.character(gnmFile$GenbankID[which(gnmFile$Genome==GenomeNm[i])[1]]),".gb")),error=function(x){NA})
    rfsqFile=tryCatch(readGenBank(file = str_c(gbkFolder,str_c(as.character(gnmFile$GenbankID[which(gnmFile$Genome==GenomeNm[i])[1]])),".gb")),error=function(x){NA})
    
    # saveRDS(gbkFile,str_c("H:/",gnmFile$Organism[i],"_ori_gbkobj.rds"))
    # saveRDS(rfsqFile,str_c("H:/",gnmFile$Organism[i],"_ori_rfsqobj.rds"))
    
    ### Read annotation files:
    
    geneDf=data.frame(chromosome=as.character(tmpExcl$Overall_GenomeID))
    
    geneDf$GeneRange=as.character(unlist(lapply(tmpExcl$str_Crd,function(x){
      
      xx=unlist(str_split(string = x,pattern ="_" ));
      xx=str_c(xx[1:2],collapse = "-");
      return(xx)
    })))
      
    geneDf$WhichStrand=as.character(unlist(lapply(tmpExcl$str_Crd,function(x){
      
      xx=unlist(str_split(string = x,pattern ="_" ));
      xx=xx[3];
      return(xx)
    })))
    
    strandMissing=which(is.na(geneDf$WhichStrand))
    if(length(strandMissing)> 0)
    {
    geneDf$WhichStrand[strandMissing]=as.character(tmpExcl$Strand.y[strandMissing]  )
      
    }
    geneDf$GeneStart= as.numeric( unlist(lapply(tmpExcl$str_Crd,function(x){
      
      xx=unlist(str_split(string = x,pattern ="_" ));
      xx=xx[1];
      return(xx)
    })))
    geneDf$GeneEnd= as.numeric(unlist(lapply(tmpExcl$str_Crd,function(x){
      
      xx=unlist(str_split(string = x,pattern ="_" ));
      xx=xx[2];
      return(xx)
    })))
    
    geneDf$type="gene"
    geneDf$locus_tag=tmpExcl$locus_tag
    geneDf$locus_tag[which(is.na(geneDf$locus_tag))]=tmpExcl$str_Crd[which(is.na(geneDf$locus_tag))]
    geneDf$gene_id=geneDf$locus_tag
    
    GRanges_obj=makeGRangesFromDataFrame(geneDf,
                             keep.extra.columns=TRUE,
                             ignore.strand=FALSE,
                             seqinfo=NULL,
                             seqnames.field=c("chromosome"),
                             start.field="GeneStart",
                             end.field="GeneEnd",
                             strand.field="WhichStrand",
                             starts.in.df.are.0based=FALSE)

    
    if(!is.na(gbkFile))
    {
      
      gbkFile@genes=GRanges_obj
      gbkFile@cds=GRanges_obj
      gbkFile@exons=GRanges_obj
      gbkFile@transcripts=GRanges_obj
      saveRDS(object = gbkFile,file = str_c("S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/",GenomeNm[i],"_Overallannotations_grange_obj.rds"))
      
    } else {
      
      rfsqFile@genes=GRanges_obj
      rfsqFile@cds=GRanges_obj
      rfsqFile@exons=GRanges_obj
      rfsqFile@transcripts=GRanges_obj
      saveRDS(object = rfsqFile,file = str_c("S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/",GenomeNm[i],"_Overallannotations_grange_obj.rds"))
      
      
    }
    
    
write.fasta(sequences = proteinsList,names = names(proteinsList),file.out = str_c("S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/",GenomeNm[i],"_FinalCompiled_proteins.faa"))
    
     
  }
  
  
  
}