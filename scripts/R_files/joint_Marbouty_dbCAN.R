joint_Marbouty_dbCAN=function(list_cgc_files,gnmNmsThis)
{
  library(foreach)
  
  list_cgc_files=sort(list_cgc_files)
  cgc_out=read.csv(file = list_cgc_files[3],sep = "\t",header = F)
  cgc_tc_tf_out=tryCatch(read.csv(file = list_cgc_files[1],sep = "\t",header = F) ,error=function(x){NA} )
  cgc_tf_out=tryCatch ( read.csv(file = list_cgc_files[2],sep = "\t",header = F)   ,error=function(x){NA} )
  cgc_families=read.csv(file = "S:/Abilash/Marbouty_annotation_refseq_genbank 11_2021/Data/Marbouty_annotations_for_integrated_proteins/dbCAN2/CAZyDB.07292021.fam.subfam.ec.txt",header =F,sep = "\t" )
  colnames(cgc_families)=c("CAZy_family","Uniprot","EC")
  overall_cgc=cgc_out
  overall_cgc$V5=str_c("CGC_",overall_cgc$V5)
  
  if(!is.na(cgc_tc_tf_out))
  {
    cgc_tc_tf_out$V5=str_c("CGC_TC_TF_",cgc_tc_tf_out$V5)
    overall_cgc=rbind(overall_cgc,cgc_tc_tf_out)
    
  }
  
  if(!is.na(cgc_tf_out))
  {
    cgc_tf_out$V5=str_c("CGC_TF_",cgc_tf_out$V5)
    overall_cgc=rbind(overall_cgc,cgc_tf_out)
    
  }
  
  overall_cgc=overall_cgc[-which(overall_cgc$V1=="+++++"),]
  overall_cgc$V13=str_c(gnmNmsThis,"_",overall_cgc$V7,"_",overall_cgc$V8)
  overall_cgc=overall_cgc[,c(2,5:13)]
  colnames(overall_cgc)=c("DB","CGCID","GenomeID","Start","End","GeneID","Strand","Description","Product_CGC","ID")
  overall_cgc$Product_CGC=as.character(overall_cgc$Product_CGC)
  
  cazyIdx=which(overall_cgc$DB=="CAZyme")
  TCIdx=which(overall_cgc$DB=="TC")
  TFIdx=which(overall_cgc$DB=="TF")
  
  ## For CAZyme
  
  overall_cgc$Product_CGC[cazyIdx]=str_replace_all(string = overall_cgc$Product_CGC[cazyIdx],pattern = ";ID=.*",replacement = "")
  overall_cgc$Product_CGC[cazyIdx]=str_replace_all(string = overall_cgc$Product_CGC[cazyIdx],pattern = "DB=",replacement = "")
  
  lapply_ec=lapply(overall_cgc$Product_CGC[cazyIdx],function(x){table(unlist(str_split(as.character(cgc_families$EC)[match(unlist(str_split(x,pattern = "\\|")),as.character(cgc_families$CAZy_family))],pattern = "\\|")))})
  overall_cgc$Product_CGC[cazyIdx]=unlist(lapply(lapply_ec,function(x){names(which(x==max(x))[1])}))
  
  ### For TC
  
  overall_cgc$Product_CGC[TCIdx]=str_replace_all(string = overall_cgc$Product_CGC[TCIdx],pattern = ";ID=.*",replacement = "")
  overall_cgc$Product_CGC[TCIdx]=str_replace_all(string = overall_cgc$Product_CGC[TCIdx],pattern = "DB=gnl\\|TC-DB\\|",replacement = "")
  
  
  ### For TF
  
  overall_cgc$Product_CGC[TFIdx]=str_replace_all(string = overall_cgc$Product_CGC[TFIdx],pattern = ";ID=.*",replacement = "")
  overall_cgc$Product_CGC[TFIdx]=str_replace_all(string = overall_cgc$Product_CGC[TFIdx],pattern = "DB=",replacement = "")
  overall_cgc$Product_CGC[TFIdx]=str_replace_all(string = overall_cgc$Product_CGC[TFIdx],pattern = "tr\\|",replacement = "")
  overall_cgc$Product_CGC=str_replace_all(string = overall_cgc$Product_CGC,pattern = "sp\\|",replacement = "")
  
  ### Annotating UniProt
  uniprot_IP=(lapply(overall_cgc$Product_CGC,function(x){str_split(string = x,pattern = "\\|")[[1]]}))
  uniprot_IP_digCnt=lapply(uniprot_IP,function(x){str_count(x, "[0-9]")})
  uniprot_IP_TotCnt=lapply(uniprot_IP,function(x){str_count(x)})
  
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
  uniprot_names=rep("",length(uniprot_IP))
  uniprot_compounds=rep("",length(uniprot_IP))
  ov_uniprot_names=list()
  ov_uniprot_compounds=list()
  
  library("stringr")
  library(RCurl)
  kk=1
  
  foreach (i = 1: length(uniprot_IP)) %do%
  {
      
    if(any(uniprot_IP[[i]]!=""))
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
  overall_cgc$Product_Name=uniprot_names
  overall_cgc$Compound=uniprot_compounds
  return(overall_cgc)
}

# Databases:
# http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/
# http://bcb.unl.edu/dbcan_pul/Webserver/static/DBCAN-PUL/dbCAN-PUL/PUL0003.out/

# Script-assistance:
# https://github.com/WatsonLab/PULpy/blob/master/Snakefile
# https://bi.snu.ac.kr/Courses/bio02/HMMER_tutorial.pdf
# https://github.com/gpertea/gsrc/blob/master/scripts/pfam_scan.pl
# https://github.com/WatsonLab/PULpy/blob/master/scripts/parse.py


