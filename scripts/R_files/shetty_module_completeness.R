shettyModulesFile=readLines("S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Extra/Shetty_et_al_MDbMM16-master/rna_seq/data_files/ModulesCustom/curatedGMMs.v1.07_new.txt")
shettyModules_slash_idx=grep(pattern = "///",x = shettyModulesFile)
shettyModules=shettyModulesFile[shettyModules_slash_idx+1]
shettyModules=shettyModules[!is.na(shettyModules)]
FinKOList=list()
FinMod=vector()
for(i in 1: length(shettyModules_slash_idx))
{
  if((i+1)>length(shettyModules_slash_idx))
  {
    ll=length(shettyModulesFile)
    
  } else{
    ll=shettyModules_slash_idx[i+1]
    
  }
  
  ovList=shettyModulesFile[shettyModules_slash_idx[i]:ll]
  ovList=ovList[ovList!="///"]
  modIdx=grep(pattern = "^M",x = ovList)
  ovModule=ovList[modIdx]
  ovKO=as.list(ovList[-modIdx])
  
  print(str_c(i," : ", length(ovModule)))
  FinKOList[[i]]=str_replace_all(unlist(ovKO),pattern = "\t$",replacement = "")
  FinMod[i]=ovModule
  
}

names(FinKOList)=FinMod
newFinMatshtv2Files=list.files(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/",pattern = "_new_FinMatShtv2.rds",full.names = T)
newFinMatshtv2_genome=str_replace_all(string = newFinMatshtv2Files,pattern = ".*/",replacement = "")
newFinMatshtv2_genome=str_replace_all(string = newFinMatshtv2_genome,pattern = "_new.*",replacement = "")

for( gnm1 in 1: length(newFinMatshtv2Files))
{
  newFinMatshtv2=readRDS(file = newFinMatshtv2Files[gnm1])
  koVec=unique(unlist(str_split(string = newFinMatshtv2$KO_overall,pattern = "///")))
  modList=lapply(FinKOList,function(x){
    
    unique(unlist(str_split(string = x,pattern = "\t|,| ")))
  })
  
  Genome_FinKOList=FinKOList# [grep(pattern = str_c(koVec,collapse = "|"),x = modList)]
  
  tot_blocks=rep(0,length(Genome_FinKOList))
  pr_blocks=rep(0,length(Genome_FinKOList))
  for(i in 1: length(Genome_FinKOList))
  {
    tot_blocks[i]=length(Genome_FinKOList[[i]])
    tmpPr=0
    
    for(j in 1: length(Genome_FinKOList[[i]]))
    {
      tmpKO=str_trim(Genome_FinKOList[[i]][j],side = "both")
      tabPr=unique(unlist(str_locate_all(string = tmpKO,pattern = "\t")))
      cmPr=unique(unlist(str_locate_all(string = tmpKO,pattern = ",")))
      
      if(length(tabPr) > 0 && length(cmPr) == 0)
      {
        ko_list=unlist(str_split(string = tmpKO,pattern = "\t"))
        
        if(any(ko_list %in% koVec))
        {
          pr_blocks[i]=pr_blocks[i]+1    
        }
      } else {
        
        if( length(tabPr) ==0 && length(cmPr) > 0)
        {
          ko_list=unlist(str_split(string = tmpKO,pattern = ","))
          
          if(all(ko_list %in% koVec))
          {
            pr_blocks[i]=pr_blocks[i]+1    
            
          }
          
        } else {
          
          if(length(tabPr) > 0 && length(cmPr) > 0)
          {
            
            tmpKO_expr=str_replace_all(str_replace_all(string = tmpKO,pattern = "\t",replacement = "|"),pattern = ",",replacement = "&")
            tmpKO_elements=unlist(str_split(string = tmpKO_expr,pattern = "\\||&"))
            tmpKO_sym=unlist(str_split(pattern = "",string = tmpKO_expr))[unique(unlist(str_locate_all(string = tmpKO_expr,pattern = "\\||&")))]
            
            # res=""
            # for(kk in 1: length(tmpKO_sym))
            # {
            #   if(tmpKO_sym[kk]=="|" ) # && kk <length(tmpKO_sym)
            #   {
            #     if(any(c(tmpKO_elements[kk],tmpKO_elements[kk+1]) %in% koVec))
            #     {
            #       res="yes"
            #      
            #       next;
            #     }
            #    
            #    
            #   }
            #  
            #   if(tmpKO_sym[kk] == "&")
            #   {
            #     if(res=="yes")
            #     {
            #       if(tmpKO_elements[kk+1] %in% koVec)
            #       {
            #         res="yes"
            #       } else {
            #         res="no"
            #       }
            #     }
            #    
            #   }
            #  
            # }
            #
            # if(res =="yes")
            # {
            #   pr_blocks[i]=pr_blocks[i]+1    
            #  
            # }
            #
            
            tmp_and_split=unlist(str_split(string = tmpKO_expr,pattern = "&"))
            tmp_and_split_dec=rep("",length(tmp_and_split))
            
            for(t1 in 1: length(tmp_and_split))
            {
              
              tmp_t2=unlist(str_split(string = tmp_and_split[t1],pattern = "\\|"))
              
              if(any(tmp_t2 %in% koVec))
              {
                
                tmp_and_split_dec[t1]="yes"
              } else {
                tmp_and_split_dec[t1]="no"
              }
              
              
              
              
            }
            
            if(all(tmp_and_split_dec =="yes"))
            {
              pr_blocks[i]=pr_blocks[i]+1
            }
            
            
          } else {
            
            if(tmpKO %in% koVec)
            {
              
              pr_blocks[i]=pr_blocks[i]+1    
              
            }
          }
        }
        
      }
      
      
      
    }
    
    
  }
  blocks_df=data.frame(total_blocks=tot_blocks,present_blocks=pr_blocks)
  blocks_df$fr_blocks=blocks_df$present_blocks/blocks_df$total_blocks
  blocks_df$Modules=names(Genome_FinKOList)
  blocks_df[c('Module_id','Module_name')]=str_split_fixed(string = blocks_df$Modules,pattern = "\t",n=2)
  saveRDS(object = blocks_df,file = str_c("S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/",newFinMatshtv2_genome[gnm1],"_shetty_module_completeness.rds"))
}
