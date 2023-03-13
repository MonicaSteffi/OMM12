remove_duplicate_genes=function(tmpOverallFinGff,GenbankFin,refseqFin)
{
  
unqGenbank_genomes=unique(GenbankFin$Locus_tag)  
  
k3=1

for(g1 in 1: length(unqGenbank_genomes))
{
  tmpIdx=grep(x = tmpOverallFinGff$attributes,pattern = unqGenbank_genomes[g1])
  
  if(length(tmpIdx) > 0)
  {
    tmpDf=tmpOverallFinGff[tmpIdx,]
    
    if(length(tmpIdx) !=1)
    {
      
      
      if(length(unique(tmpDf$source)) != nrow(tmpDf))
      {
        unqSources=unique(tmpDf$source)    
        u1_df=tmpDf
        
        for(u1 in 1: length(unqSources))
        {
          tmp_u1=which(u1_df$source==unqSources[u1])
          
          if(length(tmp_u1) ==1)
          {
            unqSc_df=u1_df[tmp_u1,]
            
          } else {
            
            tmp_u1_df=u1_df[tmp_u1,]
            unqSc_df=data.frame(seqid=str_c(unique(tmp_u1_df$seqid),collapse = "&"),
                                source=str_c(unique(tmp_u1_df$source),collapse = "&"),
                                type=str_c(unique(tmp_u1_df$type),collapse = "&"),
                                start=min(unlist(c(tmp_u1_df$start,tmp_u1_df$end))),
                                end=max(unlist(c(tmp_u1_df$start,tmp_u1_df$end)))
                                )
            unqSc_df$score=""
            unqSc_df$strand=unique(tmp_u1_df$strand)
            unqSc_df$phase=""
            attrList=unlist(lapply(tmp_u1_df$attributes,function(x){
              
              unlist(str_split(string = x,pattern = ";"));
              
            }))
            attrList=str_c(attrList,collapse = ";")
            unqSc_df$attributes=attrList
            unqSc_df$var=str_c(unique(tmp_u1_df$var),collapse = "&")
            
          }
        
          if(u1==1)
          {
            Fin_u1_df=unqSc_df
          } else {
            
            Fin_u1_df=rbind(Fin_u1_df,unqSc_df)
          }
          
        }
        tmpDf=Fin_u1_df
      } 
    }
   
    if(k3==1)
    {
      NoDupGff=tmpDf
      k3=k3+1
    } else {
      
      NoDupGff=rbind(NoDupGff,tmpDf)
    }
  }
 
  
}

unqRefSeq_genomes=unique(refseqFin$Locus_tag)  

k3=1

for(g1 in 1: length(unqRefSeq_genomes))
{
  tmpIdx=grep(x = tmpOverallFinGff$attributes,pattern = unqRefSeq_genomes[g1])
  
  if(length(tmpIdx) > 0)
  {
    tmpDf=tmpOverallFinGff[tmpIdx,]
    
    if(length(tmpIdx) !=1)
    {
      
      
      if(length(unique(tmpDf$source)) != nrow(tmpDf))
      {
        unqSources=unique(tmpDf$source)    
        u1_df=tmpDf
        
        for(u1 in 1: length(unqSources))
        {
          tmp_u1=which(u1_df$source==unqSources[u1])
          
          if(length(tmp_u1) ==1)
          {
            unqSc_df=u1_df[tmp_u1,]
            
          } else {
            
            tmp_u1_df=u1_df[tmp_u1,]
            unqSc_df=data.frame(seqid=str_c(unique(tmp_u1_df$seqid),collapse = "&"),
                                source=str_c(unique(tmp_u1_df$source),collapse = "&"),
                                type=str_c(unique(tmp_u1_df$type),collapse = "&"),
                                start=min(unlist(c(tmp_u1_df$start,tmp_u1_df$end))),
                                end=max(unlist(c(tmp_u1_df$start,tmp_u1_df$end)))
            )
            unqSc_df$score=""
            unqSc_df$strand=unique(tmp_u1_df$strand)
            unqSc_df$phase=""
            attrList=unlist(lapply(tmp_u1_df$attributes,function(x){
              
              unlist(str_split(string = x,pattern = ";"));
              
            }))
            attrList=str_c(attrList,collapse = ";")
            unqSc_df$attributes=attrList
            unqSc_df$var=str_c(unique(tmp_u1_df$var),collapse = "&")
            
          }
          
          if(u1==1)
          {
            Fin_u1_df=unqSc_df
          } else {
            
            Fin_u1_df=rbind(Fin_u1_df,unqSc_df)
          }
          
        }
        tmpDf=Fin_u1_df
      } 
    }
    if(k3==1)
    {
      NoDupGff_rfsq=tmpDf
      k3=k3+1
    } else {
      
      NoDupGff_rfsq=rbind(NoDupGff_rfsq,tmpDf)
    }
  }
  
  
}

FinNoDupGff=rbind(NoDupGff,NoDupGff_rfsq)
FinNoDupGff=FinNoDupGff[!duplicated(FinNoDupGff),]

return(FinNoDupGff)

}