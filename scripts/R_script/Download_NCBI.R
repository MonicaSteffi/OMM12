genbnk_dwnld=function()
  {

########################################################################################################################
##### --- New ----------------------------------------------------------------------------------------------------------
########################################################################################################################

library("dplyr")
library("data.table")
library("RCurl")


## Download the assembly summary genbank file from ncbi
## wget https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt


colnames_assembly_summary=c(
  "assembly_accession","bioproject","biosample","wgs_master","refseq_category", "taxid", "species_taxid", "organism_name", "infraspecific_name","isolate", "version_status","assembly_level","release_type","genome_rep",
  "seq_rel_date","asm_name","submitter","gbrs_paired_asm","paired_asm_comp", "ftp_path","excluded_from_refseq","relation_to_type_material"
)
assembly_summary=fread("/dss/dssfs02/lwp-dss-0001/u7b03/u7b03-dss-0000/ra72seh/DB/assembly_summary_genbank.txt",skip = 2,sep = "\t",quote = "")
colnames(assembly_summary)=colnames_assembly_summary

assembly_summary$strain_name=str_replace_all(string = assembly_summary$infraspecific_name,pattern =  "strain=",replacement = "")

req_speciesName=c("Akkermansia muciniphila","Acutalibacter muris","Bifidobacterium animalis","Bacteroides caecimuris","Blautia coccoides","Enterocloster clostridioformis","Clostridium innocuum","Enterococcus faecalis","Flavonifractor plautii","Limosilactobacillus reuteri","Muribaculum intestinale","Turicimonas muris")
req_strainName=c("YL44","KB18","YL2","I48","YL58","YL32","I46","KB1","YL31","I49","YL27","YL45")



req_Org=str_c(req_speciesName,req_strainName,sep = "_")
req_Org=str_replace_all(string = req_Org,pattern = " ",replacement = "_")

baseFolder="/dss/dssfs02/lwp-dss-0001/u7b03/u7b03-dss-0000/ra72seh/DB/Marbouty/"

if (!file.exists(baseFolder)){
  dir.create(baseFolder)
}

for(spNr in 1:length(req_speciesName))
{
  rwSpecies=grep(pattern = str_replace_all(string = req_speciesName[spNr],pattern = " ",replacement = "|"),x = assembly_summary$organism_name)
  
  rwStrains=grep(pattern = str_c(req_strainName[spNr],"$"),x=assembly_summary$strain_name)
  
  reqRows=intersect(rwSpecies,rwStrains)
  LtRw=order(as.Date(assembly_summary$seq_rel_date[reqRows]),decreasing = T)
  LtRw=reqRows[LtRw[1]]
  ftp_lnk=assembly_summary$ftp_path[LtRw]
  http_lnk=str_replace_all(string = ftp_lnk,pattern = "ftp://",replacement = "https://")
  lnk=str_c(http_lnk,"/")

  
  filenames <- getURL( lnk , dirlistonly=T )
  filenames <- ( strsplit(filenames, "\r*\n")[[1]] )
  filenames=str_replace_all(string = filenames,pattern = "</a>.*",replacement = "")
  filenames=str_replace_all(string = filenames,pattern = ".*>",replacement = "")
  filenames=filenames[filenames!=""]
  fileLinks=str_c(lnk,filenames)
  fileLinks=fileLinks[-grep(pattern = "Parent Directory",x=fileLinks)]
  
  
  #loop through all of those files and save them to your working directory
  for ( i in 1: length(fileLinks )){ 
    
    #determine the year directory
    pthName=str_c(baseFolder,req_Org[spNr])
    #if the directory doesn't exist, make it!
    if (!file.exists(pthName)){
      dir.create(pthName)
    }
    fName=str_replace_all(string = fileLinks[i],pattern = ".*/",replacement = "")
    download.file( fileLinks[i] , destfile =  paste(  pthName , "/" , fName,sep = "" ) )
  }	  

  
    
  }

}

genbnk_dwnld()