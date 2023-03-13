coordinate_based_non_redundant_genes=function(dC_gff_eN_O_N_P_R_KO_Sht)
{
  dC_gff_eN_O_N_P_R_KO_Sht[c('GeneStart','GeneEnd','GeneStrand')]=str_split_fixed(string = dC_gff_eN_O_N_P_R_KO_Sht$str_Crd,pattern = "_",n = 3)
  startTable=as.data.frame(table(dC_gff_eN_O_N_P_R_KO_Sht$GeneStart))
  stpTable=as.data.frame(table(dC_gff_eN_O_N_P_R_KO_Sht$GeneEnd))
  
  UC_startTable=startTable[which(startTable$Freq>1),]
  UC_stpTable=stpTable[which(stpTable$Freq>1),]
  
  UC_startTableIdx=unlist(lapply(UC_startTable$Var1,function(x){which(dC_gff_eN_O_N_P_R_KO_Sht$GeneStart==x)}))
  UC_stpTableIdx=unlist(lapply(UC_stpTable$Var1,function(x){which(dC_gff_eN_O_N_P_R_KO_Sht$GeneEnd==x)}))
  
  tbkIdx=setdiff(1:nrow(dC_gff_eN_O_N_P_R_KO_Sht),unique(c(UC_startTableIdx,UC_stpTableIdx)))
  to_be_keptDf=dC_gff_eN_O_N_P_R_KO_Sht[tbkIdx,]
  
  
  if(dim(UC_startTable)[1]>0)
  {
    UC_startTableDf=dC_gff_eN_O_N_P_R_KO_Sht[UC_startTableIdx,]
    UC_startTableDf$GeneStart=as.numeric(  UC_startTableDf$GeneStart)
    UC_startTableDf$GeneEnd=as.numeric(  UC_startTableDf$GeneEnd)
    
    ### Genome check:
    UC_startTableDf$Overall_GenomeID=apply(UC_startTableDf[,intersect(GenomeID_columns,colnames(UC_startTableDf))],1,function(x){
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
    
    #### For those rows with same start coordinates but different stop coordinates:
    
    unq_start=unique(UC_startTableDf$GeneStart)
    to_be_joined_rows=list()
    kk=1
    
    for( j2 in 1: length(unq_start))
    {
      tmpStartIdx=which(UC_startTableDf$GeneStart==unq_start[j2])
      cmbMat=as.matrix(combinat::combn(tmpStartIdx,2))
      
      for(k2 in 1: ncol(cmbMat))
      {
        if(length(unique(UC_startTableDf$Overall_GenomeID[cmbMat[,k2]]))==1)
        {
          
          
          if(diff(UC_startTableDf$GeneEnd[cmbMat[,k2]])<=50)
          {
            to_be_joined_rows[[kk]]=cmbMat[,k2]
            kk=kk+1   
          } else {
            
            line1KO=unique(c(UC_startTableDf$KO[cmbMat[1,k2]],
                             UC_startTableDf$ko4eggNOG_KO[cmbMat[1,k2]],
                             UC_startTableDf$ko4Operon_KO[cmbMat[1,k2]]))
            
            line1KO=line1KO[!is.na(line1KO)]
            line1KO=line1KO[line1KO!=""]
            line1KO=unique(line1KO)
            
            
            line2KO=unique(c(UC_startTableDf$KO[cmbMat[2,k2]],
                             UC_startTableDf$ko4eggNOG_KO[cmbMat[2,k2]],
                             UC_startTableDf$ko4Operon_KO[cmbMat[2,k2]]))
            
            line2KO=line2KO[!is.na(line2KO)]
            line2KO=line2KO[line2KO!=""]
            line2KO=unique(line2KO)
            
            
            if(any(line1KO %in% line2KO) || any(line2KO %in% line1KO))
            {
              
              to_be_joined_rows[[kk]]=cmbMat[,k2]
              kk=kk+1   
              
            } else {
              
              line1COG=unique(c(UC_startTableDf$COG4Operon[cmbMat[1,k2]],
                                UC_startTableDf$COG4gbk_rfsq[cmbMat[1,k2]],
                                UC_startTableDf$COG4eggNOG[cmbMat[1,k2]]))
              
              line1COG=line1COG[!is.na(line1COG)]
              line1COG=line1COG[line1COG!=""]
              line1COG=unique(line1COG)
              
              
              
              line2COG=unique(c(UC_startTableDf$COG4Operon[cmbMat[2,k2]],
                                UC_startTableDf$COG4gbk_rfsq[cmbMat[2,k2]],
                                UC_startTableDf$COG4eggNOG[cmbMat[2,k2]]))
              
              line2COG=line2COG[!is.na(line2COG)]
              line2COG=line2COG[line2COG!=""]
              line2COG=unique(line2COG)
              
              
              if(any(line1COG %in% line2COG) || any(line2COG %in% line1COG))
              {
                
                to_be_joined_rows[[kk]]=cmbMat[,k2]
                kk=kk+1   
                
              }
              
              
              
            }
            
          }
        }
      }
      
    }
    
    to_be_joined_rows_overall=unique(unlist(to_be_joined_rows))
    
    to_be_joined_rows=unique(sapply(to_be_joined_rows, function(x) 
      sort(unique(unlist(to_be_joined_rows[sapply(to_be_joined_rows, function(y) 
        any(x %in% y))])))))
    res.list <- to_be_joined_rows
    # take set difference between contents of list elements and accumulated elements
    res.list[-1] <- mapply("setdiff", res.list[-1],
                           head(Reduce(c, to_be_joined_rows, accumulate=TRUE), -1))
    
    to_be_joined_rows=res.list
    retainedIdx=setdiff(1:nrow(UC_startTableDf),to_be_joined_rows_overall)
    if(length(retainedIdx) >0)
    {
      UC_startTableRetained=UC_startTableDf[retainedIdx,]
      ncbistart=ifelse(test = !is.na(UC_startTableRetained$locus_tag),UC_startTableRetained$GeneStart,"")
      ncbiend=ifelse(test = !is.na(UC_startTableRetained$locus_tag),UC_startTableRetained$GeneEnd,"")
      
      dbCAN2start=ifelse(test = !is.na(UC_startTableRetained$dbCAN2_Cluster),UC_startTableRetained$GeneStart,"")
      dbCAN2end=ifelse(test = !is.na(UC_startTableRetained$dbCAN2_Cluster),UC_startTableRetained$GeneEnd,"")
      
      eggNOGstart=ifelse(test = !is.na(UC_startTableRetained$eggNOG_OGs.x),UC_startTableRetained$GeneStart,"")
      eggNOGend=ifelse(test = !is.na(UC_startTableRetained$eggNOG_OGs.x),UC_startTableRetained$GeneEnd,"")
      
      ## added later..
      operonstart=ifelse(test = !is.na(UC_startTableRetained$Operon),UC_startTableRetained$PosLeft,"")
      
      if(length(grep(pattern = "PosRight",x=colnames(UC_startTableRetained))) > 0)
      {
        operonend=ifelse(test = !is.na(UC_startTableRetained$Operon),UC_startTableRetained$PosRight,"")
        
      }
      
      if(length(grep(pattern = "postRight",x=colnames(UC_startTableRetained))) > 0)
      {
        operonend=ifelse(test = !is.na(UC_startTableRetained$Operon),UC_startTableRetained$postRight,"")
        
      }
      
      rgistart=ifelse(test = !is.na(UC_startTableRetained$ARO),UC_startTableRetained$Start.y,"")
      rgiend=ifelse(test = !is.na(UC_startTableRetained$ARO),UC_startTableRetained$Stop.y,"")
      
      UC_startTableRetained$PULStart=""
      UC_startTableRetained$PULEnd=""
      
      if(any(colnames(UC_startTableRetained) %in% "pulid"))
      {
        pulstart=ifelse(test = !is.na(UC_startTableRetained$pulid),UC_startTableRetained$start,"")
        pulend=ifelse(test = !is.na(UC_startTableRetained$pulid),UC_startTableRetained$end,"")
        
        UC_startTableRetained$PULStart=pulstart
        UC_startTableRetained$PULEnd=pulend
      }
      
      
      UC_startTableRetained$NCBIStart=ncbistart
      UC_startTableRetained$NCBIEnd=ncbiend
      
      UC_startTableRetained$dbCAN2Start=dbCAN2start
      UC_startTableRetained$dbCAN2End=dbCAN2end
      
      UC_startTableRetained$eggNOGStart=eggNOGstart
      UC_startTableRetained$eggNOGEnd=eggNOGend
      
      UC_startTableRetained$operonStart=operonstart
      UC_startTableRetained$operonEnd=operonend
      
      ## added later..
      UC_startTableRetained$RGIStart=rgistart
      UC_startTableRetained$RGIEnd=rgiend
      
      
    } else {
      UC_startTableRetained=data.frame()
    }
    
    if(length(to_be_joined_rows) > 0)
    {
      for(jn1 in 1: length(to_be_joined_rows))
      {
        
        tmpIdx=to_be_joined_rows[[jn1]]
        tmpDf=droplevels(as.data.frame(UC_startTableDf[tmpIdx,]))
        gbkIdx=which(tmpDf$locus_tag!="")
        
        if(length(gbkIdx)!=0)
        {
          
          tmpAnn=data.frame(
            NCBIStart=min(unique(tmpDf$GeneStart[gbkIdx])),
            NCBIEnd=max(unique(tmpDf$GeneEnd[gbkIdx]))
          )
          
          
          tmpAnn$Strand=unique(apply(tmpDf[,grep(pattern = "Strand",x=colnames(tmpDf))],1,function(x){
            xx=x[!is.na(x)];
            xx=xx[xx!=""];
            return(unique(xx))
          }))
          
          tmpAnn$str_Crd=apply(tmpAnn,1,function(x){
            str_c(x,collapse = "_")
            
          })
          
          tmpAnn=tmpAnn[,c("str_Crd","NCBIStart","NCBIEnd")]
          
          dbCAN2_idx=which(!is.na(tmpDf$dbCAN2_Cluster))
          eggNOG_idx=which(!is.na(tmpDf$eggNOG_OGs.x))
          operon_idx=which(!is.na(tmpDf$Operon))
          RGI_idx=which(!is.na(tmpDf$ARO))
          
          
          if(length(dbCAN2_idx)>0)
          {
            tmpAnn$dbCAN2Start=min(tmpDf$GeneStart[dbCAN2_idx])
            tmpAnn$dbCAN2End=max(tmpDf$GeneEnd[dbCAN2_idx])
            
          } else {
            tmpAnn$dbCAN2Start=""
            tmpAnn$dbCAN2End=""
            
          }
          
          
          if(length(eggNOG_idx)>0)
          {
            tmpAnn$eggNOGStart=min(unique(tmpDf$GeneStart[eggNOG_idx]))
            tmpAnn$eggNOGEnd=max(unique(tmpDf$GeneEnd[eggNOG_idx]))
            
          } else {
            tmpAnn$eggNOGStart=""
            tmpAnn$eggNOGEnd=""
            
          }
          
          
          if(length(operon_idx)>0)
          {
            tmpAnn$operonStart=min(unique(tmpDf$GeneStart[operon_idx]))
            tmpAnn$operonEnd=max(unique(tmpDf$GeneEnd[operon_idx]))
            
          } else {
            tmpAnn$operonStart=""
            tmpAnn$operonEnd=""
            
          }
          
          
          if(length(RGI_idx)>0)
          {
            tmpAnn$RGIStart=min(unique(tmpDf$Start.y[RGI_idx]))
            tmpAnn$RGIEnd=max(unique(tmpDf$Stop.y[RGI_idx]))
            
          } else {
            tmpAnn$RGIStart=""
            tmpAnn$RGIEnd=""
            
          }
          
          tmpAnn$PULStart=""
          tmpAnn$PULEnd=""
          
          if(any(colnames(tmpDf) %in% "pulid"))
          {
            PUL_idx=which(!is.na(tmpDf$pulid))
            if(length(PUL_idx)>0)
            {
              tmpAnn$PULStart=min(unique(tmpDf$start[PUL_idx]))
              tmpAnn$PULEnd=max(unique(tmpDf$end[PUL_idx]))
              
            }
            
          }
          
          
          
          
          tmpAnn2=(apply(tmpDf[,-grep(pattern = as.character(str_c(colnames(tmpAnn),collapse = "|")),x=colnames(tmpDf))],2,function(x){
            x=x[!is.na(x)];
            str_c(unique(x),collapse = "///")
            
            
          }))
          
          
          tmpAnn2df=as.matrix(t(unname(tmpAnn2)),ncol=length(tmpAnn2))
          colnames(tmpAnn2df)=names(tmpAnn2)
          tmpAnn2df=as.data.frame(tmpAnn2df)
          
          
        } else {
          if(length(unique(tmpDf$str_Crd))==1)
          {
            tmpAnn=data.frame(str_Crd=tmpDf$str_Crd[1],
                              NCBIStart=tmpDf$GeneStart[1],
                              NCBIEnd=tmpDf$GeneEnd[1]) 
            
            
            dbCAN2_idx=which(!is.na(tmpDf$dbCAN2_Cluster))
            eggNOG_idx=which(!is.na(tmpDf$eggNOG_OGs.x))
            operon_idx=which(!is.na(tmpDf$Operon))
            
            operon_idx=which(!is.na(tmpDf$Operon))
            RGI_idx=which(!is.na(tmpDf$ARO))
            
            
            
            if(length(dbCAN2_idx)>0)
            {
              tmpAnn$dbCAN2Start=min(unique(tmpDf$GeneStart[dbCAN2_idx]))
              tmpAnn$dbCAN2End=max(unique(tmpDf$GeneEnd[dbCAN2_idx]))
              
            } else {
              tmpAnn$dbCAN2Start=""
              tmpAnn$dbCAN2End=""
              
            }
            
            
            if(length(eggNOG_idx)>0)
            {
              tmpAnn$eggNOGStart=min(unique(tmpDf$GeneStart[eggNOG_idx]))
              tmpAnn$eggNOGEnd=max(unique(tmpDf$GeneEnd[eggNOG_idx]))
              
            } else {
              tmpAnn$eggNOGStart=""
              tmpAnn$eggNOGEnd=""
              
            }
            
            
            if(length(operon_idx)>0)
            {
              tmpAnn$operonStart=min(unique(tmpDf$GeneStart[operon_idx]))
              tmpAnn$operonEnd=max(unique(tmpDf$GeneEnd[operon_idx]))
              
            } else {
              tmpAnn$operonStart=""
              tmpAnn$operonEnd=""
              
            }
            
            
            if(length(RGI_idx)>0)
            {
              tmpAnn$RGIStart=min(unique(tmpDf$Start.y[RGI_idx]))
              tmpAnn$RGIEnd=max(unique(tmpDf$Stop.y[RGI_idx]))
              
            } else {
              tmpAnn$RGIStart=""
              tmpAnn$RGIEnd=""
              
            }
            
            
            tmpAnn$PULStart=""
            tmpAnn$PULEnd=""
            
            if(any(colnames(UC_startTableRetained) %in% "pulid"))
            {
              PUL_idx=which(!is.na(tmpDf$pulid))
              
              if(length(PUL_idx)>0)
              {
                tmpAnn$PULStart=min(unique(tmpDf$start[PUL_idx]))
                tmpAnn$PULEnd=max(unique(tmpDf$end[PUL_idx]))
                
              }
            }
            
            
            tmpAnn2=(apply(tmpDf[,-grep(pattern = as.character(str_c(colnames(tmpAnn),collapse = "|")),x=colnames(tmpDf))],2,function(x){
              x=x[!is.na(x)];
              str_c(unique(x),collapse = "///")
              
              
            }))
            
            
            tmpAnn2df=as.matrix(t(unname(tmpAnn2)),ncol=length(tmpAnn2))
            colnames(tmpAnn2df)=names(tmpAnn2)
            tmpAnn2df=as.data.frame(tmpAnn2df)
            
            
            
          } else {
            
            ### This scenario is a bit exceptional:
            ### There are two genes from genbank/refseq with overlapping regions and same annotations.
            ### However, we retain both of them since they are present in genbank/refseq.
            
            if(length(tmpDf$locus_tag)==2)
            {
              tmpAnn=data.frame(str_Crd=tmpDf$str_Crd,
                                NCBIStart=tmpDf$GeneStart,
                                NCBIEnd=tmpDf$GeneEnd)
              
              dbCAN2_idx=which(!is.na(tmpDf$dbCAN2_Cluster))
              eggNOG_idx=which(!is.na(tmpDf$eggNOG_OGs.x))
              operon_idx=which(!is.na(tmpDf$Operon))
              RGI_idx=which(!is.na(tmpDf$ARO))
              
              
              if(length(dbCAN2_idx)>0)
              {
                tmpAnn$dbCAN2Start=min(unique(tmpDf$GeneStart[dbCAN2_idx]))
                tmpAnn$dbCAN2End=max(unique(tmpDf$GeneEnd[dbCAN2_idx]))
                
              } else {
                tmpAnn$dbCAN2Start=""
                tmpAnn$dbCAN2End=""
                
              }
              
              
              if(length(eggNOG_idx)>0)
              {
                tmpAnn$eggNOGStart=min(unique(tmpDf$GeneStart[eggNOG_idx]))
                tmpAnn$eggNOGEnd=max(unique(tmpDf$GeneEnd[eggNOG_idx]))
                
              } else {
                tmpAnn$eggNOGStart=""
                tmpAnn$eggNOGEnd=""
                
              }
              
              
              if(length(operon_idx)>0)
              {
                tmpAnn$operonStart=min(unique(tmpDf$GeneStart[operon_idx]))
                tmpAnn$operonEnd=max(unique(tmpDf$GeneEnd[operon_idx]))
                
              } else {
                tmpAnn$operonStart=""
                tmpAnn$operonEnd=""
                
              }
              
              
              if(length(RGI_idx)>0)
              {
                tmpAnn$RGIStart=min(unique(tmpDf$Start.y[RGI_idx]))
                tmpAnn$RGIEnd=max(unique(tmpDf$Stop.y[RGI_idx]))
                
              } else {
                tmpAnn$RGIStart=""
                tmpAnn$RGIEnd=""
                
              }
              
              
              tmpAnn$PULStart=""
              tmpAnn$PULEnd=""
              
              if(any(colnames(UC_startTableRetained) %in% "pulid"))
              {
                PUL_idx=which(!is.na(tmpDf$pulid))
                
                if(length(PUL_idx)>0)
                {
                  tmpAnn$PULStart=min(unique(tmpDf$start[PUL_idx]))
                  tmpAnn$PULEnd=max(unique(tmpDf$end[PUL_idx]))
                  
                }
              }
              
              
              tmpAnn2df=(tmpDf[,-grep(pattern = as.character(str_c(colnames(tmpAnn),collapse = "|")),x=colnames(tmpDf))])
              
              
              
              
            }
          }
          
          
        }
        
        
        # tmpAnn2df=tmpAnn2df[!duplicated(tmpAnn2df),]
        
        tmpAnnFinal=cbind(tmpAnn,tmpAnn2df)
        
        for(tmp123 in 1: ncol(tmpAnnFinal))
        {
          tmpAnnFinal[,tmp123]=as.character(tmpAnnFinal[,tmp123])
          
        }
        
        
        if(jn1 ==1)
        {
          newAnn2=tmpAnnFinal
          
        } else {
          
          newAnn2=rbind(newAnn2,tmpAnnFinal[,colnames(newAnn2)])
          
        }
      }
      
      
      newAnn2=newAnn2[!duplicated(newAnn2),]
      
      if(class(UC_startTableRetained)!="numeric")
      {
        UC_startOverall=rbind(UC_startTableRetained,
                              newAnn2[,colnames(UC_startTableRetained)])
        
      } else {
        UC_startOverall=newAnn2
      }
      
    } else {
      UC_startOverall=UC_startTableRetained
    }
    
    
    ###################################################
  }
  
  if(dim(UC_stpTable)[1]>0)
  {
    UC_stpTableDf=dC_gff_eN_O_N_P_R_KO_Sht[UC_stpTableIdx,]  
    UC_stpTableDf$GeneStart=as.numeric(  UC_stpTableDf$GeneStart)
    UC_stpTableDf$GeneEnd=as.numeric(  UC_stpTableDf$GeneEnd)
    
    ### Genome check:
    UC_stpTableDf$Overall_GenomeID=apply(UC_stpTableDf[,intersect(GenomeID_columns,colnames(UC_stpTableDf))],1,function(x){
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
    
    
    unq_stop=unique(UC_stpTableDf$GeneEnd)
    to_be_joined_rows=list()
    kk=1
    
    for( j2 in 1: length(unq_stop))
    {
      tmpStpIdx=which(UC_stpTableDf$GeneEnd==unq_stop[j2])
      cmbMat=as.matrix(combinat::combn(tmpStpIdx,2))
      
      for(k2 in 1: ncol(cmbMat))
      {
        if(length(unique(UC_stpTableDf$Overall_GenomeID[cmbMat[,k2]]))==1)
        {
          
          if(diff(UC_stpTableDf$GeneStart[cmbMat[,k2]])<=50)
          {
            to_be_joined_rows[[kk]]=cmbMat[,k2]
            kk=kk+1   
          } else {
            
            line1KO=unique(c(UC_stpTableDf$KO[cmbMat[1,k2]],
                             UC_stpTableDf$ko4eggNOG_KO[cmbMat[1,k2]],
                             UC_stpTableDf$ko4Operon_KO[cmbMat[1,k2]]))
            
            line1KO=line1KO[!is.na(line1KO)]
            line1KO=line1KO[line1KO!=""]
            line1KO=unique(line1KO)
            
            
            line2KO=unique(c(UC_stpTableDf$KO[cmbMat[2,k2]],
                             UC_stpTableDf$ko4eggNOG_KO[cmbMat[2,k2]],
                             UC_stpTableDf$ko4Operon_KO[cmbMat[2,k2]]))
            
            line2KO=line2KO[!is.na(line2KO)]
            line2KO=line2KO[line2KO!=""]
            line2KO=unique(line2KO)
            
            
            if(any(line1KO %in% line2KO) || any(line2KO %in% line1KO))
            {
              
              to_be_joined_rows[[kk]]=cmbMat[,k2]
              kk=kk+1   
              
            } else {
              
              line1COG=unique(c(UC_stpTableDf$COG4Operon[cmbMat[1,k2]],
                                UC_stpTableDf$COG4gbk_rfsq[cmbMat[1,k2]],
                                UC_stpTableDf$COG4eggNOG[cmbMat[1,k2]]))
              
              line1COG=line1COG[!is.na(line1COG)]
              line1COG=line1COG[line1COG!=""]
              line1COG=unique(line1COG)
              
              
              
              line2COG=unique(c(UC_stpTableDf$COG4Operon[cmbMat[2,k2]],
                                UC_stpTableDf$COG4gbk_rfsq[cmbMat[2,k2]],
                                UC_stpTableDf$COG4eggNOG[cmbMat[2,k2]]))
              
              line2COG=line2COG[!is.na(line2COG)]
              line2COG=line2COG[line2COG!=""]
              line2COG=unique(line2COG)
              
              
              if(any(line1COG %in% line2COG) || any(line2COG %in% line1COG))
              {
                
                to_be_joined_rows[[kk]]=cmbMat[,k2]
                kk=kk+1   
                
              }
              
              
              
            }
            
          }
          
        }
      }
      
    }
    
    to_be_joined_rows_overall=unique(unlist(to_be_joined_rows))
    
    to_be_joined_rows=unique(sapply(to_be_joined_rows, function(x) 
      sort(unique(unlist(to_be_joined_rows[sapply(to_be_joined_rows, function(y) 
        any(x %in% y))])))))
    res.list <- to_be_joined_rows
    # take set difference between contents of list elements and accumulated elements
     res.list[-1] <- mapply("setdiff", res.list[-1],
                           head(Reduce(c, to_be_joined_rows, accumulate=TRUE), -1))
    
     to_be_joined_rows=res.list
    retainedIdx=setdiff(1:nrow(UC_stpTableDf),to_be_joined_rows_overall)
    if(length(retainedIdx)>0)
    {
      UC_stpTableRetained=UC_stpTableDf[retainedIdx,]
      ncbistart=ifelse(test = !is.na(UC_stpTableRetained$locus_tag),UC_stpTableRetained$GeneStart,"")
      ncbiend=ifelse(test = !is.na(UC_stpTableRetained$locus_tag),UC_stpTableRetained$GeneEnd,"")
      
      dbCAN2start=ifelse(test = !is.na(UC_stpTableRetained$dbCAN2_Cluster),UC_stpTableRetained$GeneStart,"")
      dbCAN2end=ifelse(test = !is.na(UC_stpTableRetained$dbCAN2_Cluster),UC_stpTableRetained$GeneEnd,"")
      
      eggNOGstart=ifelse(test = !is.na(UC_stpTableRetained$eggNOG_OGs.x),UC_stpTableRetained$GeneStart,"")
      eggNOGend=ifelse(test = !is.na(UC_stpTableRetained$eggNOG_OGs.x),UC_stpTableRetained$GeneEnd,"")
      
      operonstart=ifelse(test = !is.na(UC_stpTableRetained$Operon),UC_stpTableRetained$PosLeft,"")
      
      if(length(grep(pattern = "PosRight",x=colnames(UC_stpTableRetained))) > 0)
      {
        operonend=ifelse(test = !is.na(UC_stpTableRetained$Operon),UC_stpTableRetained$PosRight,"")
        
      }
      
      if(length(grep(pattern = "postRight",x=colnames(UC_stpTableRetained))) > 0)
      {
        operonend=ifelse(test = !is.na(UC_stpTableRetained$Operon),UC_stpTableRetained$postRight,"")
        
      }
      
      rgistart=ifelse(test = !is.na(UC_stpTableRetained$ARO),UC_stpTableRetained$Start.y,"")
      rgiend=ifelse(test = !is.na(UC_stpTableRetained$ARO),UC_stpTableRetained$Stop.y,"")
      
      UC_stpTableRetained$PULStart=""
      UC_stpTableRetained$PULEnd=""
      
      
      if(any(colnames(UC_stpTableRetained) %in% "pulid"))
      {
        pulstart=ifelse(test = !is.na(UC_stpTableRetained$pulid),UC_stpTableRetained$start,"")
        pulend=ifelse(test = !is.na(UC_stpTableRetained$pulid),UC_stpTableRetained$end,"")
        
        UC_stpTableRetained$PULStart=pulstart
        UC_stpTableRetained$PULEnd=pulend
      }
      
      
      
      UC_stpTableRetained$NCBIStart=ncbistart
      UC_stpTableRetained$NCBIEnd=ncbiend
      
      UC_stpTableRetained$dbCAN2Start=dbCAN2start
      UC_stpTableRetained$dbCAN2End=dbCAN2end
      
      UC_stpTableRetained$eggNOGStart=eggNOGstart
      UC_stpTableRetained$eggNOGEnd=eggNOGend
      
      UC_stpTableRetained$operonStart=operonstart
      UC_stpTableRetained$operonEnd=operonend
      
      UC_stpTableRetained$RGIStart=rgistart
      UC_stpTableRetained$RGIEnd=rgiend
      
    } else {
      
      UC_stpTableRetained=data.frame()
    }
    
    if(length(to_be_joined_rows) > 0)
    {
      for(jn1 in 1: length(to_be_joined_rows))
      {
        
        tmpIdx=to_be_joined_rows[[jn1]]
        tmpDf=droplevels(as.data.frame(UC_stpTableDf[tmpIdx,]))
        gbkIdx=which(tmpDf$locus_tag!="")
        
        if(length(gbkIdx)!=0)
        {
          
          tmpAnn=data.frame(
                            NCBIStart=min(unique(tmpDf$GeneStart[gbkIdx])),
                            NCBIEnd=max(unique(tmpDf$GeneEnd[gbkIdx]))
                            )
          
          
          tmpAnn$Strand=unique(apply(tmpDf[,grep(pattern = "Strand",x=colnames(tmpDf))],1,function(x){
            xx=x[!is.na(x)];
            xx=xx[xx!=""];
            return(unique(xx))
          }))
          
          tmpAnn$str_Crd=apply(tmpAnn,1,function(x){
            str_c(x,collapse = "_")
            
          })
          
          tmpAnn=tmpAnn[,c("str_Crd","NCBIStart","NCBIEnd")]
          dbCAN2_idx=which(!is.na(tmpDf$dbCAN2_Cluster))
          eggNOG_idx=which(!is.na(tmpDf$eggNOG_OGs.x))
          operon_idx=which(!is.na(tmpDf$Operon))
          RGI_idx=which(!is.na(tmpDf$ARO))
          
          
          if(length(dbCAN2_idx)>0)
          {
            tmpAnn$dbCAN2Start=min(unique(tmpDf$GeneStart[dbCAN2_idx]))
            tmpAnn$dbCAN2End=max(unique(tmpDf$GeneEnd[dbCAN2_idx]))
            
          } else {
            tmpAnn$dbCAN2Start=""
            tmpAnn$dbCAN2End=""
            
          }
          
          
          if(length(eggNOG_idx)>0)
          {
            tmpAnn$eggNOGStart=min(unique(tmpDf$GeneStart[eggNOG_idx]))
            tmpAnn$eggNOGEnd=max(unique(tmpDf$GeneEnd[eggNOG_idx]))
            
          } else {
            tmpAnn$eggNOGStart=""
            tmpAnn$eggNOGEnd=""
            
          }
          
          
          if(length(operon_idx)>0)
          {
            tmpAnn$operonStart=min(unique(tmpDf$GeneStart[operon_idx]))
            tmpAnn$operonEnd=max(unique(tmpDf$GeneEnd[operon_idx]))
            
          } else {
            tmpAnn$operonStart=""
            tmpAnn$operonEnd=""
            
          }
          
          
          if(length(RGI_idx)>0)
          {
            tmpAnn$RGIStart=min(unique(tmpDf$Start.y[RGI_idx]))
            tmpAnn$RGIEnd=max(unique(tmpDf$Stop.y[RGI_idx]))
            
          } else {
            tmpAnn$RGIStart=""
            tmpAnn$RGIEnd=""
            
          }
          
          tmpAnn$PULStart=""
          tmpAnn$PULEnd=""
          
          if(any(colnames(UC_stpTableRetained) %in% "pulid"))
          {
            PUL_idx=which(!is.na(tmpDf$pulid))
            
            if(length(PUL_idx)>0)
            {
              tmpAnn$PULStart=min(unique(tmpDf$start[PUL_idx]))
              tmpAnn$PULEnd=max(unique(tmpDf$end[PUL_idx]))
              
            }
          }
          
          
          
          tmpAnn2=(apply(tmpDf[,-grep(pattern = as.character(str_c(colnames(tmpAnn),collapse = "|")),x=colnames(tmpDf))],2,function(x){
            x=x[!is.na(x)];
            str_c(unique(x),collapse = "///")
            
            
          }))
          
          
          tmpAnn2df=as.matrix(t(unname(tmpAnn2)),ncol=length(tmpAnn2))
          colnames(tmpAnn2df)=names(tmpAnn2)
          tmpAnn2df=as.data.frame(tmpAnn2df)
          
          
        } else {
          if(length(unique(tmpDf$str_Crd))==1)
          {
            tmpAnn=data.frame(str_Crd=tmpDf$str_Crd[1],
                              NCBIStart=tmpDf$GeneStart[1],
                              NCBIEnd=tmpDf$GeneEnd[1]) 
            
            
            dbCAN2_idx=which(!is.na(tmpDf$dbCAN2_Cluster))
            eggNOG_idx=which(!is.na(tmpDf$eggNOG_OGs.x))
            operon_idx=which(!is.na(tmpDf$Operon))
            operon_idx=which(!is.na(tmpDf$Operon))
            RGI_idx=which(!is.na(tmpDf$ARO))
            
            
            
            
            if(length(dbCAN2_idx)>0)
            {
              tmpAnn$dbCAN2Start=min(unique(tmpDf$GeneStart[dbCAN2_idx]))
              tmpAnn$dbCAN2End=max(unique(tmpDf$GeneEnd[dbCAN2_idx]))
              
            } else {
              tmpAnn$dbCAN2Start=""
              tmpAnn$dbCAN2End=""
              
            }
            
            
            if(length(eggNOG_idx)>0)
            {
              tmpAnn$eggNOGStart=min(unique(tmpDf$GeneStart[eggNOG_idx]))
              tmpAnn$eggNOGEnd=max(unique(tmpDf$GeneEnd[eggNOG_idx]))
              
            } else {
              tmpAnn$eggNOGStart=""
              tmpAnn$eggNOGEnd=""
              
            }
            
            
            if(length(operon_idx)>0)
            {
              tmpAnn$operonStart=min(unique(tmpDf$GeneStart[operon_idx]))
              tmpAnn$operonEnd=max(unique(tmpDf$GeneEnd[operon_idx]))
              
            } else {
              tmpAnn$operonStart=""
              tmpAnn$operonEnd=""
              
            }
            
            
            if(length(RGI_idx)>0)
            {
              tmpAnn$RGIStart=min(unique(tmpDf$Start.y[RGI_idx]))
              tmpAnn$RGIEnd=max(unique(tmpDf$Stop.y[RGI_idx]))
              
            } else {
              tmpAnn$RGIStart=""
              tmpAnn$RGIEnd=""
              
            }
            
            tmpAnn$PULStart=""
            tmpAnn$PULEnd=""
            
            if(any(colnames(UC_stpTableRetained) %in% "pulid"))
            {
              PUL_idx=which(!is.na(tmpDf$pulid))
              
              if(length(PUL_idx)>0)
              {
                tmpAnn$PULStart=min(unique(tmpDf$start[PUL_idx]))
                tmpAnn$PULEnd=max(unique(tmpDf$end[PUL_idx]))
                
              }
            }
            
            
            
            
            tmpAnn2=(apply(tmpDf[,-grep(pattern = as.character(str_c(colnames(tmpAnn),collapse = "|")),x=colnames(tmpDf))],2,function(x){
              x=x[!is.na(x)];
              str_c(unique(x),collapse = "///")
              
              
            }))
            
            
            tmpAnn2df=as.matrix(t(unname(tmpAnn2)),ncol=length(tmpAnn2))
            colnames(tmpAnn2df)=names(tmpAnn2)
            tmpAnn2df=as.data.frame(tmpAnn2df)
            
            
            
          } else {
            
            ### This scenario is a bit exceptional:
            ### There are two genes from genbank/refseq with overlapping regions and same annotations.
            ### However, we retain both of them since they are present in genbank/refseq.
            
            if(length(tmpDf$locus_tag)==2)
            {
              tmpAnn=data.frame(str_Crd=tmpDf$str_Crd,
                                NCBIStart=tmpDf$GeneStart,
                                NCBIEnd=tmpDf$GeneEnd)
              
              dbCAN2_idx=which(!is.na(tmpDf$dbCAN2_Cluster))
              eggNOG_idx=which(!is.na(tmpDf$eggNOG_OGs.x))
              operon_idx=which(!is.na(tmpDf$Operon))
              
              operon_idx=which(!is.na(tmpDf$Operon))
              RGI_idx=which(!is.na(tmpDf$ARO))
              
              
              
              
              if(length(dbCAN2_idx)>0)
              {
                tmpAnn$dbCAN2Start=min(unique(tmpDf$GeneStart[dbCAN2_idx]))
                tmpAnn$dbCAN2End=max(unique(tmpDf$GeneEnd[dbCAN2_idx]))
                
              } else {
                tmpAnn$dbCAN2Start=""
                tmpAnn$dbCAN2End=""
                
              }
              
              
              if(length(eggNOG_idx)>0)
              {
                tmpAnn$eggNOGStart=min(unique(tmpDf$GeneStart[eggNOG_idx]))
                tmpAnn$eggNOGEnd=max(unique(tmpDf$GeneEnd[eggNOG_idx]))
                
              } else {
                tmpAnn$eggNOGStart=""
                tmpAnn$eggNOGEnd=""
                
              }
              
              
              if(length(operon_idx)>0)
              {
                tmpAnn$operonStart=min(unique(tmpDf$GeneStart[operon_idx]))
                tmpAnn$operonEnd=max(unique(tmpDf$GeneEnd[operon_idx]))
                
              } else {
                tmpAnn$operonStart=""
                tmpAnn$operonEnd=""
                
              }
              
              
              if(length(RGI_idx)>0)
              {
                tmpAnn$RGIStart=min(unique(tmpDf$Start.y[RGI_idx]))
                tmpAnn$RGIEnd=max(unique(tmpDf$Stop.y[RGI_idx]))
                
              } else {
                tmpAnn$RGIStart=""
                tmpAnn$RGIEnd=""
                
              }
              
              tmpAnn$PULStart=""
              tmpAnn$PULEnd=""
              
              if(any(colnames(UC_stpTableRetained) %in% "pulid"))
              {
                PUL_idx=which(!is.na(tmpDf$pulid))
                
                if(length(PUL_idx)>0)
                {
                  tmpAnn$PULStart=min(unique(tmpDf$start[PUL_idx]))
                  tmpAnn$PULEnd=max(unique(tmpDf$end[PUL_idx]))
                  
                }
              }
              
              
              tmpAnn2df=(tmpDf[,-grep(pattern = as.character(str_c(colnames(tmpAnn),collapse = "|")),x=colnames(tmpDf))])
              
              
              
              
            }
          }
          
          
        }
        
        
        # tmpAnn2df=tmpAnn2df[!duplicated(tmpAnn2df),]
        
        tmpAnnFinal=cbind(tmpAnn,tmpAnn2df)
        
        for(tmp123 in 1: ncol(tmpAnnFinal))
        {
          tmpAnnFinal[,tmp123]=as.character(tmpAnnFinal[,tmp123])
          
        }
        if(jn1 ==1)
        {
          newAnn2=tmpAnnFinal
          
        } else {
          
          newAnn2=rbind(newAnn2,tmpAnnFinal[,colnames(newAnn2)])
          
        }
      }
      
      
      newAnn2=newAnn2[!duplicated(newAnn2),]
      
      if(class(UC_stpTableRetained)!="numeric")
      {
        UC_stpOverall=rbind(UC_stpTableRetained,
                            newAnn2[,colnames(UC_stpTableRetained)])
        
      } else {
        UC_stpOverall=newAnn2
      }
      
    } else {
      UC_stpOverall=UC_stpTableRetained
    }
    
  }
  
  if(length(tbkIdx)>0)
  {
    
    tbkIdx=setdiff(1:nrow(dC_gff_eN_O_N_P_R_KO_Sht),unique(c(UC_startTableIdx,UC_stpTableIdx)))
    to_be_keptDf=dC_gff_eN_O_N_P_R_KO_Sht[tbkIdx,]
    
    
    to_be_keptDf$Overall_GenomeID=apply(to_be_keptDf[,intersect(GenomeID_columns,colnames(to_be_keptDf))],1,function(x){
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
    
    if(length(tbkIdx) >0)
    {
      
      ncbistart=ifelse(test = !is.na(to_be_keptDf$locus_tag),to_be_keptDf$GeneStart,"")
      ncbiend=ifelse(test = !is.na(to_be_keptDf$locus_tag),to_be_keptDf$GeneEnd,"")
      
      dbCAN2start=ifelse(test = !is.na(to_be_keptDf$dbCAN2_Cluster),to_be_keptDf$GeneStart,"")
      dbCAN2end=ifelse(test = !is.na(to_be_keptDf$dbCAN2_Cluster),to_be_keptDf$GeneEnd,"")
      
      eggNOGstart=ifelse(test = !is.na(to_be_keptDf$eggNOG_OGs.x),to_be_keptDf$GeneStart,"")
      eggNOGend=ifelse(test = !is.na(to_be_keptDf$eggNOG_OGs.x),to_be_keptDf$GeneEnd,"")
      
      operonstart=ifelse(test = !is.na(to_be_keptDf$Operon),to_be_keptDf$PosLeft,"")
      
      if(length(grep(pattern = "PosRight",x=colnames(to_be_keptDf))) > 0)
      {
        operonend=ifelse(test = !is.na(to_be_keptDf$Operon),to_be_keptDf$PosRight,"")
        
      }
      
      if(length(grep(pattern = "postRight",x=colnames(to_be_keptDf))) > 0)
      {
        operonend=ifelse(test = !is.na(to_be_keptDf$Operon),to_be_keptDf$postRight,"")
        
      }
      
      
      RGIstart=ifelse(test = !is.na(to_be_keptDf$ARO),to_be_keptDf$Start.y,"")
      RGIend=ifelse(test = !is.na(to_be_keptDf$ARO),to_be_keptDf$Stop.y,"")
      
      to_be_keptDf$NCBIStart=ncbistart
      to_be_keptDf$NCBIEnd=ncbiend
      
      to_be_keptDf$dbCAN2Start=dbCAN2start
      to_be_keptDf$dbCAN2End=dbCAN2end
      
      to_be_keptDf$eggNOGStart=eggNOGstart
      to_be_keptDf$eggNOGEnd=eggNOGend
      
      to_be_keptDf$operonStart=operonstart
      to_be_keptDf$operonEnd=operonend
      
      to_be_keptDf$RGIStart=RGIstart
      to_be_keptDf$RGIEnd=RGIend
      
      to_be_keptDf$PULStart=""
      to_be_keptDf$PULEnd=""
      
      if(any(colnames(to_be_keptDf) %in% "pulid"))
      {
        pulstart=ifelse(test = !is.na(to_be_keptDf$pulid),to_be_keptDf$start,"")
        pulend=ifelse(test = !is.na(to_be_keptDf$pulid),to_be_keptDf$end,"")
        
        to_be_keptDf$PULStart=pulstart
        to_be_keptDf$PULEnd=pulend
        
        
      }
      
    } else {
      to_be_keptDf=data.frame()
    }
    
  }
  
  return(new(Class = "non_redundant_genes",to_be_keptDf=to_be_keptDf,
             UC_stpOverall=UC_stpOverall,
             UC_startOverall=UC_startOverall))
  
}