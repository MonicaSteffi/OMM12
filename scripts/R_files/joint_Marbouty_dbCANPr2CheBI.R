joint_Marbouty_dbCANPr2CheBI=function()
{
  
  library("stringr")
  
  library("RCurl")
  
  library("data.table")
  
  library(foreach)
  
  
  kofamPath="S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/KOFAMSCAN_results_new//"
  gffPath="S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/gffFiles/"
  ProteinPath="S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/"
  OperonPath="S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/OperonMapper_new/"
  cgcPath="S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/dbCAN2_new/"
  
  
  list_koFiles=list.files(path = kofamPath,full.names = T,pattern = "_result_all_tab_edited")
  gnmNms=str_replace_all(string = list_koFiles,pattern = ".*/",replacement = "")
  gnmNms=str_replace_all(string = gnmNms,pattern = "_result.*",replacement = "")
  
  
  
  
  
  
  rm_onlyDigits=function(digCnt,TotCnt,x)
    
  {
    
    rtString=vector()
    
    k=1
    
    for(i in 1:length(digCnt))
      
    {
      if((TotCnt[i]) != (digCnt[i]))
        
      {
        rtString[k]=x[i]
        k=k+1
      }
      
    }
    
    return(rtString)
  }
  
  
  
  for(gnmNmsThis in gnmNms)
  {
    
    FinalAnn=readRDS(str_c("S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/",gnmNmsThis,"_allAnnJoint.rds"))
    
    ### Annotating UniProt
    
    uniprot_IP=(lapply(FinalAnn$Product_CGC,function(x){str_split(string = x,pattern = "\\|")[[1]]}))
    
    uniprot_IP_digCnt=lapply(uniprot_IP,function(x){str_count(x, "[0-9]")})
    
    uniprot_IP_TotCnt=lapply(uniprot_IP,function(x){str_count(x)})
    
    
    uniprot_names=rep("",length(uniprot_IP))
    
    uniprot_compounds=rep("",length(uniprot_IP))
    
    ov_uniprot_names=list()
    
    ov_uniprot_compounds=list()
    
    
    library("stringr")
    
    library(RCurl)
    
    kk=1
    
    
    for(i in  1: length(uniprot_IP)) # %do%
      
    {
      
      
      if(any(uniprot_IP[[i]]!="") && !is.na(unlist(uniprot_IP[[i]])))
        
      {
        
        message(str_c("++",i ,"--->",uniprot_IP[[i]], " is running!"))
        
        ind_up=unlist(str_split(string = uniprot_IP[[i]],pattern = "\\|"))
        
        ind_up=rm_onlyDigits(digCnt = uniprot_IP_digCnt[[i]],TotCnt = uniprot_IP_TotCnt[[i]],x = uniprot_IP[[i]])  
        
        ind_up=ind_up[ind_up!=""]
        
        
        ## remove only numbered-strings
        
        if(length(ind_up)>0)
          
        {
          
          listIdx=grep(ind_up,names(ov_uniprot_compounds))
          
          if(length(listIdx)>0)
            
          {
            
            uniprot_names[i]=str_c(unique(unlist(ov_uniprot_names[listIdx])),collapse = ";")
            
            uniprot_compounds[i]=str_c(unique(unlist(ov_uniprot_compounds[match(ind_up,names(ov_uniprot_compounds))])),collapse = ";")
            
            break
            
          } else {
            
            for(jj in 1:length(ind_up)) 
              
            {
              
              message("++++",str_c(jj ,"--->",ind_up[jj], " is running!"))
              
              tmpind_up=which(ind_up==ind_up[jj])
              
              
              URL1=str_c("https://www.uniprot.org/uniprot/?query=",ind_up[jj],"+AND+taxonomy%3ABacteria&sort=score&format=tab&columns=protein+names,chebi")# AND+reviewed%3Ayes+ # id,entry+name,
              # x <- getURL(URL1)  
              
              # tmpFlNm=str_replace_all(ind_up[jj],pattern = "[[:punct:]]",replacement = "_")
              # download.file(destfile = str_c("/dss/dssfs02/lwp-dss-0001/u7b03/u7b03-dss-0000/ra72seh/",tmpFlNm),url = URL1,method = "curl",headers = T,mode = "wb",) # auto
              # out=tryCatch(read.csv(str_c("/dss/dssfs02/lwp-dss-0001/u7b03/u7b03-dss-0000/ra72seh/",tmpFlNm),sep="\t",header = T),error=function(x){"###"})
              # file.remove(str_c("/dss/dssfs02/lwp-dss-0001/u7b03/u7b03-dss-0000/ra72seh/",tmpFlNm))
              # if(out!="###")
              out=tryCatch(fread(URL1,sep = "\t",quote = F,header = T,stringsAsFactors = FALSE,colClasses = c("character","character"),blank.lines.skip = TRUE),error=function(x){"###"})
              
              if(any(class(out)=="data.frame")==TRUE)
              {
                unqCmpds=unique(str_trim(unlist(str_split(names(sort(table(out[,2]),decreasing = T)),pattern = ";")),side = "both"))
                unqCmpds=unqCmpds[unqCmpds!=""]
                
                uniprot_names[tmpind_up]=    names(sort(table(out[,1]),decreasing = T)[1])
                uniprot_compounds[tmpind_up]=    str_c(unqCmpds,collapse = ";")
                
                tmpNames=list(rep(uniprot_names[tmpind_up],length(ind_up)))
                tmpCompounds=list(rep(uniprot_compounds[tmpind_up],length(ind_up)))
                names(tmpNames)=ind_up
                names(tmpCompounds)=ind_up
                
                ov_uniprot_compounds=c(ov_uniprot_compounds,tmpCompounds)
                ov_uniprot_names=c(ov_uniprot_names,tmpNames)
                
                break
              } 
              
            }
            
            
          }
          
          
          
          
          
          
          
        }
        
      }
      
    }
    FinalAnn$Product_Name=uniprot_names
    FinalAnn$Compound=uniprot_compounds    
    saveRDS(FinalAnn,str_c("S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/",gnmNmsThis,"_FinalAnnWithCheBIUniprot.rds"))
  }
  
  
  # return(FinalAnn)
}


joint_Marbouty_dbCANPr2CheBI()