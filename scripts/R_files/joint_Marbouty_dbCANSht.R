joint_Marbouty_dbCANSht=function(list_cgc_filesThis,gnmNmsThis)
{
  ## Downloaded the background database from https://bcb.unl.edu/dbCAN2/download/Databases/V11/
  
  cgc_families=read.csv(file = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/dbCAN2_new/CAZyDB.08062022.fam.subfam.ec.txt",header =F,sep = "\t" )
  colnames(cgc_families)=c("CAZy_family","Uniprot","EC")
  
  
  library(foreach)
  library(vroom)

  list_cgc_filesThis=sort(list_cgc_filesThis)
  xx=grep(pattern = "cgc_standard.out.txt",x=list_cgc_filesThis)
  
  
  if(gnmNmsThis != "escherichia_coli_strain_Mt1B1")
  {
    
  
  
  
  cgc_out=read.csv(file = list_cgc_filesThis[xx],sep = "\t",header = F)
  
  
  
  colnames(cgc_out)=c("dbCAN2_Cluster","GeneType","GenomeID","dbCAN2_GeneID","Cluster_Start","Cluster_End","Strand","dbCAN2_annotation")
  cgc_out$dbCAN2_Cluster=str_c("CGC_TC_",cgc_out$dbCAN2_Cluster)
  
  cgc_tc_tf_outFiles=list_cgc_filesThis[grep(pattern = "/cgc_tc_tf/",x=list_cgc_filesThis)]
  cgc_tf_outFiles=list_cgc_filesThis[grep(pattern = "/cgc_tf/",x=list_cgc_filesThis)]
  
  if(length(cgc_tc_tf_outFiles)> 0)
  {
    
  cgc_tc_tf_outList=lapply(cgc_tc_tf_outFiles,function(x){read.csv(file = x,header = F,sep = "\t")})
  cgc_tc_tf_outList=lapply(cgc_tc_tf_outList,function(x){x=x[,-1]; return(x)})
  cgc_tc_tf_outList=lapply(cgc_tc_tf_outList,function(x)
    {
    
    if(ncol(x)==11)
    {
      x[,10]=apply(x[,c(10,11)],1,function(x){str_c(unique(x),collapse = ";")})
      x=x[,-11]
    } 
    
    return(x)
    
    })
  
  cgc_tc_tf_outDf=(do.call(what = rbind,args = cgc_tc_tf_outList))
  colnames(cgc_tc_tf_outDf)=c("GeneType","Downstream_distance","Upstream_distance","dbCAN2_Cluster","GenomeID","Cluster_Start","Cluster_End","dbCAN2_GeneID","Strand","dbCAN2_annotation")
  cgc_tc_tf_outDf$dbCAN2_Cluster=str_c("CGC_TC_TF_",cgc_tc_tf_outDf$dbCAN2_Cluster)
  
  } else {
  
    cgc_tc_tf_outDf=as.data.frame(GeneType= character(), Downstream_distance=numeric(),Upstream_distance=numeric(),
                               Cluster=character(), GenomeID=character(), Cluster_Start=numeric(), Cluster_End=numeric(),
                               dbCAN2_GeneID=character(),Strand=character(),dbCAN2_annotation=character())  
    
}
  
  
  
  if(length(cgc_tf_outFiles)> 0)
  {
    cgc_tf_outList=lapply(cgc_tf_outFiles,function(x){read.csv(file = x,header = F,sep = "\t")})
    cgc_tf_outList=lapply(cgc_tf_outList,function(x){x=x[,-1]; return(x)})
    cgc_tf_outList=lapply(cgc_tf_outList,function(x)
    {
      
      if(ncol(x)==11)
      {
        x[,10]=apply(x[,c(10,11)],1,function(x){str_c(unique(x),collapse = ";")})
        x=x[,-11]
      } 
      
      return(x)
      
    })
    
    cgc_tf_outDf=(do.call(what = rbind,args = cgc_tf_outList))
    colnames(cgc_tf_outDf)=c("GeneType","Downstream_distance","Upstream_distance","dbCAN2_Cluster","GenomeID","Cluster_Start","Cluster_End","dbCAN2_GeneID","Strand","dbCAN2_annotation")
    cgc_tf_outDf$dbCAN2_Cluster=str_c("CGC_TF_",cgc_tf_outDf$dbCAN2_Cluster)
    
  } else {
    cgc_tf_outDf=as.data.frame(GeneType= character(), Downstream_distance=numeric(),Upstream_distance=numeric(),
                               Cluster=character(), GenomeID=character(), Cluster_Start=numeric(), Cluster_End=numeric(),
                               dbCAN2_GeneID=character(),Strand=character(),dbCAN2_annotation=character())
  }
  
  
  
  overall_cgc=cgc_out
  overall_cgc=rbind(overall_cgc,cgc_tc_tf_outDf[,colnames(overall_cgc)],cgc_tf_outDf[,colnames(overall_cgc)])
  overall_cgc$CGCID=apply(overall_cgc[,c("GenomeID" , "Cluster_Start","Cluster_End","Strand")],1,function(x){str_c(x,collapse = "_")})
  
  overall_cgc$dbCAN2_annotation=as.character(overall_cgc$dbCAN2_annotation)
  overall_cgc$dbCAN2_annotation_Ec=""
  overall_cgc$dbCAN2_annotation_Uniprot=""
  
  
  cazyIdx=which(overall_cgc$GeneType=="CAZyme")
  TCIdx=which(overall_cgc$GeneType=="TC")
  TFIdx=which(overall_cgc$GeneType=="TF")
  STPIdx=which(overall_cgc$GeneType=="STP")
  
  
  ## For CAZyme
  
  overall_cgc$dbCAN2_annotation[cazyIdx]=str_replace_all(string = overall_cgc$dbCAN2_annotation[cazyIdx],pattern = ";ID=.*",replacement = "")
  overall_cgc$dbCAN2_annotation[cazyIdx]=str_replace_all(string = overall_cgc$dbCAN2_annotation[cazyIdx],pattern = "DB=",replacement = "")
  
  lapply_ec=lapply(overall_cgc$dbCAN2_annotation[cazyIdx],function(x){table(unlist(str_split(as.character(cgc_families$EC)[match(unlist(str_split(x,pattern = "\\|")),as.character(cgc_families$CAZy_family))],pattern = "\\|")))})
  
  lapply_ec=lapply(overall_cgc$dbCAN2_annotation[cazyIdx],function(x) {
    unqXX=unique(unlist(str_split(string = x,pattern = "\\|")));
    yy=unique(unlist(str_split(cgc_families$EC[match(unqXX,cgc_families$CAZy_family )],pattern = "\\|")))
  str_c(yy[!is.na(yy)],collapse = "|")
  })
  
  lapply_uniprot=lapply(overall_cgc$dbCAN2_annotation[cazyIdx],function(x) {
    unqXX=unique(unlist(str_split(string = x,pattern = "\\|")));
    yy=unique(unlist(str_split(cgc_families$Uniprot[match(unqXX,cgc_families$CAZy_family )],pattern = "\\|")))
    str_c(yy[!is.na(yy)],collapse = "|")
  })
  
  overall_cgc$dbCAN2_annotation_Ec[cazyIdx]=unlist(lapply_ec)
  overall_cgc$dbCAN2_annotation_Uniprot[cazyIdx]=unlist(lapply_uniprot)
  
  ### For TC
  
  overall_cgc$dbCAN2_annotation[TCIdx]=str_replace_all(string = overall_cgc$dbCAN2_annotation[TCIdx],pattern = ";ID=.*",replacement = "")
  overall_cgc$dbCAN2_annotation[TCIdx]=str_replace_all(string = overall_cgc$dbCAN2_annotation[TCIdx],pattern = "DB=gnl\\|TC-DB\\|",replacement = "")
  
  
  ### For TF
  
  overall_cgc$dbCAN2_annotation[TFIdx]=str_replace_all(string = overall_cgc$dbCAN2_annotation[TFIdx],pattern = ";ID=.*",replacement = "")
  overall_cgc$dbCAN2_annotation[TFIdx]=str_replace_all(string = overall_cgc$dbCAN2_annotation[TFIdx],pattern = "DB=",replacement = "")
  overall_cgc$dbCAN2_annotation[TFIdx]=str_replace_all(string = overall_cgc$dbCAN2_annotation[TFIdx],pattern = "tr\\|",replacement = "")
  overall_cgc$dbCAN2_annotation=str_replace_all(string = overall_cgc$dbCAN2_annotation,pattern = "sp\\|",replacement = "")
  overall_cgcFin=overall_cgc

  } else {
    
    
    
    for(xC in 1: length(xx))
    {
      
  ###########################
      
      cgc_out=read.csv(file = list_cgc_filesThis[xx[xC]],sep = "\t",header = F)
      
      
      
      colnames(cgc_out)=c("dbCAN2_Cluster","GeneType","GenomeID","dbCAN2_GeneID","Cluster_Start","Cluster_End","Strand","dbCAN2_annotation")
      cgc_out$dbCAN2_Cluster=str_c("CGC_TC_",cgc_out$dbCAN2_Cluster)
      
      cgc_tc_tf_outFiles=list_cgc_filesThis[grep(pattern = "/cgc_tc_tf/",x=list_cgc_filesThis)]
      cgc_tf_outFiles=list_cgc_filesThis[grep(pattern = "/cgc_tf/",x=list_cgc_filesThis)]
      
      if(length(cgc_tc_tf_outFiles)> 0)
      {
        
        cgc_tc_tf_outList=lapply(cgc_tc_tf_outFiles,function(x){read.csv(file = x,header = F,sep = "\t")})
        cgc_tc_tf_outList=lapply(cgc_tc_tf_outList,function(x){x=x[,-1]; return(x)})
        cgc_tc_tf_outList=lapply(cgc_tc_tf_outList,function(x)
        {
          
          if(ncol(x)==11)
          {
            x[,10]=apply(x[,c(10,11)],1,function(x){str_c(unique(x),collapse = ";")})
            x=x[,-11]
          } 
          
          return(x)
          
        })
        
        cgc_tc_tf_outDf=(do.call(what = rbind,args = cgc_tc_tf_outList))
        colnames(cgc_tc_tf_outDf)=c("GeneType","Downstream_distance","Upstream_distance","dbCAN2_Cluster","GenomeID","Cluster_Start","Cluster_End","dbCAN2_GeneID","Strand","dbCAN2_annotation")
        cgc_tc_tf_outDf$dbCAN2_Cluster=str_c("CGC_TC_TF_",cgc_tc_tf_outDf$dbCAN2_Cluster)
        
      } else {
        
        cgc_tc_tf_outDf=as.data.frame(GeneType= character(), Downstream_distance=numeric(),Upstream_distance=numeric(),
                                      Cluster=character(), GenomeID=character(), Cluster_Start=numeric(), Cluster_End=numeric(),
                                      dbCAN2_GeneID=character(),Strand=character(),dbCAN2_annotation=character())  
        
      }
      
      
      
      if(length(cgc_tf_outFiles)> 0)
      {
        cgc_tf_outList=lapply(cgc_tf_outFiles,function(x){read.csv(file = x,header = F,sep = "\t")})
        cgc_tf_outList=lapply(cgc_tf_outList,function(x){x=x[,-1]; return(x)})
        cgc_tf_outList=lapply(cgc_tf_outList,function(x)
        {
          
          if(ncol(x)==11)
          {
            x[,10]=apply(x[,c(10,11)],1,function(x){str_c(unique(x),collapse = ";")})
            x=x[,-11]
          } 
          
          return(x)
          
        })
        
        cgc_tf_outDf=(do.call(what = rbind,args = cgc_tf_outList))
        colnames(cgc_tf_outDf)=c("GeneType","Downstream_distance","Upstream_distance","dbCAN2_Cluster","GenomeID","Cluster_Start","Cluster_End","dbCAN2_GeneID","Strand","dbCAN2_annotation")
        cgc_tf_outDf$dbCAN2_Cluster=str_c("CGC_TF_",cgc_tf_outDf$dbCAN2_Cluster)
        
      } else {
        cgc_tf_outDf=as.data.frame(GeneType= character(), Downstream_distance=numeric(),Upstream_distance=numeric(),
                                   Cluster=character(), GenomeID=character(), Cluster_Start=numeric(), Cluster_End=numeric(),
                                   dbCAN2_GeneID=character(),Strand=character(),dbCAN2_annotation=character())
      }
      
      
      
      overall_cgc=cgc_out
      overall_cgc=rbind(overall_cgc,cgc_tc_tf_outDf[,colnames(overall_cgc)],cgc_tf_outDf[,colnames(overall_cgc)])
      overall_cgc$CGCID=apply(overall_cgc[,c("GenomeID" , "Cluster_Start","Cluster_End","Strand")],1,function(x){str_c(x,collapse = "_")})
      
      overall_cgc$dbCAN2_annotation=as.character(overall_cgc$dbCAN2_annotation)
      overall_cgc$dbCAN2_annotation_Ec=""
      overall_cgc$dbCAN2_annotation_Uniprot=""
      
      
      cazyIdx=which(overall_cgc$GeneType=="CAZyme")
      TCIdx=which(overall_cgc$GeneType=="TC")
      TFIdx=which(overall_cgc$GeneType=="TF")
      STPIdx=which(overall_cgc$GeneType=="STP")
      
      
      ## For CAZyme
      
      overall_cgc$dbCAN2_annotation[cazyIdx]=str_replace_all(string = overall_cgc$dbCAN2_annotation[cazyIdx],pattern = ";ID=.*",replacement = "")
      overall_cgc$dbCAN2_annotation[cazyIdx]=str_replace_all(string = overall_cgc$dbCAN2_annotation[cazyIdx],pattern = "DB=",replacement = "")
      
      lapply_ec=lapply(overall_cgc$dbCAN2_annotation[cazyIdx],function(x){table(unlist(str_split(as.character(cgc_families$EC)[match(unlist(str_split(x,pattern = "\\|")),as.character(cgc_families$CAZy_family))],pattern = "\\|")))})
      
      lapply_ec=lapply(overall_cgc$dbCAN2_annotation[cazyIdx],function(x) {
        unqXX=unique(unlist(str_split(string = x,pattern = "\\|")));
        yy=unique(unlist(str_split(cgc_families$EC[match(unqXX,cgc_families$CAZy_family )],pattern = "\\|")))
        str_c(yy[!is.na(yy)],collapse = "|")
      })
      
      lapply_uniprot=lapply(overall_cgc$dbCAN2_annotation[cazyIdx],function(x) {
        unqXX=unique(unlist(str_split(string = x,pattern = "\\|")));
        yy=unique(unlist(str_split(cgc_families$Uniprot[match(unqXX,cgc_families$CAZy_family )],pattern = "\\|")))
        str_c(yy[!is.na(yy)],collapse = "|")
      })
      
      overall_cgc$dbCAN2_annotation_Ec[cazyIdx]=unlist(lapply_ec)
      overall_cgc$dbCAN2_annotation_Uniprot[cazyIdx]=unlist(lapply_uniprot)
      
      ### For TC
      
      overall_cgc$dbCAN2_annotation[TCIdx]=str_replace_all(string = overall_cgc$dbCAN2_annotation[TCIdx],pattern = ";ID=.*",replacement = "")
      overall_cgc$dbCAN2_annotation[TCIdx]=str_replace_all(string = overall_cgc$dbCAN2_annotation[TCIdx],pattern = "DB=gnl\\|TC-DB\\|",replacement = "")
      
      
      ### For TF
      
      overall_cgc$dbCAN2_annotation[TFIdx]=str_replace_all(string = overall_cgc$dbCAN2_annotation[TFIdx],pattern = ";ID=.*",replacement = "")
      overall_cgc$dbCAN2_annotation[TFIdx]=str_replace_all(string = overall_cgc$dbCAN2_annotation[TFIdx],pattern = "DB=",replacement = "")
      overall_cgc$dbCAN2_annotation[TFIdx]=str_replace_all(string = overall_cgc$dbCAN2_annotation[TFIdx],pattern = "tr\\|",replacement = "")
      overall_cgc$dbCAN2_annotation=str_replace_all(string = overall_cgc$dbCAN2_annotation,pattern = "sp\\|",replacement = "")
      
      
      
  ###########################    
      if(xC ==1) 
      {
        
        overall_cgcFin=overall_cgc
      } else {
        
        overall_cgcFin=rbind(overall_cgcFin,overall_cgc)
      }
    }
    
  }
  
  overall_cgcFin$str_Crd=apply(overall_cgcFin[,c("Cluster_Start", "Cluster_End","Strand")],1,function(x){str_c(x,collapse = "_")})
  return(overall_cgcFin)
  
  }


