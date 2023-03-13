library("stringr")
library("data.table")
gnkList=list.files(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Downloaded_from_genbank_and_refseq/Marbouty_Genbank/",full.names = T,recursive = T,pattern = "_feature_table.txt")
rfsqList=list.files(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Downloaded_from_genbank_and_refseq/Marbouty_refseq/",full.names = T,recursive = T,pattern = "_feature_table.txt")
ftPttrn=str_replace_all(string = gnkList,pattern = ".*Marbouty_Genbank//",replacement = "")
ftPttrn=str_replace_all(string = ftPttrn,pattern = "/GC.*",replacement = "")
ftPttrn=str_replace_all(string = ftPttrn,pattern = "/",replacement = "")
ovList=c(gnkList,rfsqList)

for(i in 1: length(ftPttrn))
{
  
FtFiles=ovList[grep(pattern = ftPttrn[i],x=ovList)]  
dbIdx=grep(pattern = "Marbouty_refseq",x=FtFiles)
dbNm=c("Genbank","Genbank")
dbNm[dbIdx]="RefSeq"

for(j in 1:length(FtFiles))
{
  ff=fread(file = FtFiles[j],sep = "\t",skip = 1)
  colnames(ff)=c("feature",	"class",	"assembly"	,"assembly_unit",	"seq_type"	,"chromosome",	"genomic_accession",
                 "start",	"end",	"strand",	"product_accession",	"non-redundant_refseq",	"related_accession",
                 "name",	"symbol",	"GeneID","locus_tag",	"feature_interval_length",	"product_length"	,"attributes")
  ff$str_crd=apply(ff[,c(8,9,10)],1,function(x){str_c(str_trim(as.character(x),side = "both"),collapse="_")})
  ff$genomic_accession=str_replace_all(string = ff$genomic_accession,pattern = "NZ_",replacement = "")
 
  unqstr_crd=unique(ff$str_crd)
 for(k in 1: length(unqstr_crd))
 {
   tmpCrd=which(ff$str_crd==unqstr_crd[k])
   fv=data.frame(chromosome=ff$genomic_accession[tmpCrd[1]],start=ff$start[tmpCrd[1]],end=ff$end[tmpCrd[1]],
                 strand=ff$strand[tmpCrd[1]],Gene=ff$attributes[tmpCrd[1]],Protein=ff$locus_tag[tmpCrd[1]],Crd=ff$str_crd[tmpCrd[1]])
   fv$Gene=str_replace_all(string = fv$Gene,pattern = "old_locus_tag=",replacement = "")
   if(k==1)
   {
     ovFv=fv
     
   } else {
     
     ovFv=rbind(ovFv,fv)
   }
   
   
 }
 ovFv$db=dbNm[j]  
 
 if(j==1)
 {
   Fv=ovFv
   
 } else {
   Fv=rbind(Fv,ovFv)
   
 }
 
 }


FvGnbk=Fv[Fv$db=="Genbank",]
FvRfsq=Fv[Fv$db=="RefSeq",]
overallFv=FvGnbk
overallFv=rbind(overallFv,FvRfsq[match(setdiff(FvRfsq$Crd,FvGnbk$Crd),FvRfsq$Crd),])
write.table(overallFv[,c(1:4,6,5)],file = str_c("S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/PULs_new/",ftPttrn[i],"_joint_feature_table.txt"),quote = F,sep = "\t",row.names = F,col.names = F)
message(str_c(i, "x is done!!"))
}



########!/usr/bin/env -S Rscript --vanilla
#!/usr/bin/env Rscript

#################################
#
# Functions
#
#################################

process_pairs <- function(idx, td) {
  
  mindbc <- dbcan[dbcan$protein_id %in% td$protein_id,]
  
  mindbc <- unique(mindbc[,c("protein_id","hmm")])
  
  gh <- data.frame(protein_id="",hmm="", stringsAsFactors=FALSE)
  if (nrow(mindbc) > 0) {
    gh <- aggregate(mindbc$hmm, by=list(protein_id=mindbc$protein_id), paste, collapse=";")
  }
  colnames(gh) <- c("protein_id","hmm")
  
  cmin = idx -1
  while(cmin > 0) {
    
    contains = 0
    for (j in cmin:(cmin-4)) {
      
      if (j>0) {
        
        dist <- td$start[j+1] - td$end[j]
        if (dist > 500) {
          break
        }
        
        # check if it is susC/D
        if (td$sus[j] == "susC" || td$sus[j] == "susD") {
          contains = 1
        }
        
        # check if it is a GH
        prot_id = td$protein_id[j]
        ghvec <- gh$hmm[gh$protein_id==prot_id]
        if (length(ghvec) > 0) {
          contains = 1
        }
        
      }
      
    }
    
    if (contains > 0) {
      cmin <- cmin -1
    } else {
      break
    }
  }
  
  cmin <- cmin + 1
  
  cmax = idx + 1
  while(1) {
    
    contains = 0
    for (j in cmax:(cmax+4)) {
      
      if (j <= nrow(td)) {
        
        dist <- td$start[j] - td$end[j-1]
        if (dist > 500) {
          break
        }
        
        # check if it is susC/D
        if (td$sus[j] == "susC" || td$sus[j] == "susD") {
          contains = 1
        }
        
        # check if it is a GH
        prot_id = td$protein_id[j]
        ghvec <- gh$hmm[gh$protein_id==prot_id]
        if (length(ghvec) > 0) {
          contains = 1
        }
      }
      
    }
    
    if (contains > 0) {
      cmax <- cmax + 1
    } else {
      break
    }
  }
  
  cmax <- cmax - 1
  
  if (cmax > cmin) {
    
    pultdf <- merge(td[cmin:cmax,], gh, all.x=TRUE, sort=FALSE)
    pultdf <- pultdf[order(pultdf$order),]
    pultdf[is.na(pultdf)] <- ""
    pultdf$pulid <- rep(paste("PUL",PULCOUNTER, sep=""), nrow(pultdf))
    
    print(paste("created PUL", PULCOUNTER))
    print(paste(cmin,"-",cmax))
    
    ALLPULS <<- rbind(ALLPULS, pultdf)
    
    PULCOUNTER <<- PULCOUNTER + 1
  }
  return(cmax)
  
}

#################################
#
# code taken and edited from https://github.com/WatsonLab/PULpy
#
#################################
hmmFile="S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/dbCAN2_new/dbCAN-HMMdb-V11.txt"
PULdb="S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/PULs_new/"


db_names=list.files(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/PULs_new/",pattern = "V11_dbCANhmmscan_output.filt",full.names = T)
pfam_names=list.files(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/PULs_new/",pattern = ".pfam",full.names = T)
ft_names=list.files(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/Data/Marbouty_annotations_for_integrated_proteins/PULs_new/",pattern = "_joint_feature_table.txt",recursive = T,full.names = T)

genomes_list=str_replace_all(string = pfam_names,pattern = ".*PULs_new/",replacement = "")
genomes_list=str_replace_all(string = genomes_list,pattern = "\\.pfam",replacement = "")
genomes_list=str_c(genomes_list,"V11")
pul_tsv_names=str_c(PULdb,genomes_list,"_pul.tsv")
pul_tsv_sum_names=str_c(PULdb,genomes_list,"_pul.sum.tsv")

# get command line arguments as an array
# args <- commandArgs(trailingOnly = TRUE)
for (i5 in 1:length(db_names))
{
  # arguments
  ft.name <- ft_names[i5] # args[1]
  pf.name <- pfam_names[i5] # args[2]
  db.name <- db_names[i5] # args[3]
  
  # file base
  file.out <- pul_tsv_names[i5] # args[4]
  file.sum <- pul_tsv_sum_names[i5] # args[5]
  
  # genome
  genome <- genomes_list[i5] # args[6]
  
  cn <- c("contig","start","end","strand","protein_id","protein_name")
  ft <- read.table(ft.name, 
                   header=FALSE, 
                   sep="\t", 
                   stringsAsFactors=FALSE, 
                   col.names=cn)
  
  ft$order = 1:nrow(ft)
  
  cn <- c("protein_id","as","ae","es","ee","ha","hn","type","hs","he","hl","bit","e","sig","clan","active")
  pfam <- read.table(pf.name, 
                     skip=29, header=FALSE, 
                     stringsAsFactors=FALSE, 
                     fill=TRUE, 
                     col.names=cn)
  
  cn <- c("hmm","hit_len","protein_id","query_len","evalue","hit_start","hit_end","query_start","query_end","cov")
  dbcan <- read.table(db.name, 
                      sep="\t", 
                      header=FALSE, 
                      stringsAsFactors=FALSE, 
                      col.names=cn)
  
  dbcan$hmm <- gsub(".hmm","",dbcan$hmm)
  
  # GLOBAL VARIABLES
  PULCOUNTER <- 1
  ALLPULS <- data.frame(protein_id=character(),
                        contig=character(),
                        start=numeric(),
                        end=numeric(),
                        strand=character(),
                        protein_name=character(),
                        order=numeric(), 
                        active=character(), 
                        sus=character(), 
                        hmm=character(), 
                        pulid=character(),
                        stringsAsFactors=FALSE)
  
  
  # filter pfam
  halign_prop = (pfam$he - pfam$hs + 1) / pfam$hl
  pfam <- pfam[halign_prop>=0.6,]
  
  # find and annotate susC/susD
  pfam$sus <- rep("none", nrow(pfam))
  pfam$sus[grep("PF00593",pfam$ha)] <- "susC"
  pfam$sus[grep("PF13715",pfam$ha)] <- "susC"
  pfam$sus[grep("PF07715",pfam$ha)] <- "susC"
  pfam$sus[grep("PF07980",pfam$ha)] <- "susD"
  pfam$sus[grep("PF12741",pfam$ha)] <- "susD"
  pfam$sus[grep("PF12771",pfam$ha)] <- "susD"
  pfam$sus[grep("PF14322",pfam$ha)] <- "susD"
  
  # limit pfam to relevant rows
  pfam <- pfam[pfam$sus != "none",]
  
  # merge with feature table
  ftp <- merge(ft, pfam, by.x="protein_id", by.y="protein_id", all.x=TRUE, sort=FALSE)[,c(1:7,22,23)]
  ftp <- unique(ftp)
  ftp <- ftp[order(ftp$order),]
  
  # unique list of contigs
  cons <- unique(ftp$contig)
  
  # go through each contig
  for (c in cons[order(cons)]) {
    
    tdf <- ftp[ftp$contig==c,]
    tdf <- tdf
    tdf[is.na(tdf)] <- ""
    
    i <- 1
    while(i <= nrow(tdf)) {
      #for (i in 1:nrow(tdf)) {
      
      #print(c)
      #print(nrow(tdf))
      #print(i)
      if (i>=nrow(tdf)) {
        break
      }
      if (tdf$sus[i] == "susC" && tdf$sus[i+1] == "susD") {
        # we have a pair - do something
        print(paste("pair at",c,i))
        i <- process_pairs(i,tdf)
        print(i)
      }
      
      if (i>=nrow(tdf)) {
        break
      } 
      
      if (tdf$sus[i] == "susD" && tdf$sus[i+1] == "susC") {
        # we have a pair - do something
        print(paste("pair at",c, i))
        i <- process_pairs(i,tdf)
        print(i)
      }
      
      i <- i+1
    }
    
  }
  
  
  if (nrow(ALLPULS) >= 2) {
    two <- ALLPULS$end[1:(nrow(ALLPULS)-1)]
    one <- ALLPULS$start[2:nrow(ALLPULS)]
    
    ALLPULS$dist <- c(0,one - two)
    
    ALLPULS$genome <- rep(genome, nrow(ALLPULS))
    
    ALLPULS <- ALLPULS[, c("genome","pulid","protein_id","contig","start","end","strand","dist","protein_name","sus","hmm","active")]
    
    ALLPULS <- ALLPULS[order(ALLPULS$contig, ALLPULS$start),]
    
    write.table(ALLPULS, file.out, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
    
  } else {
    ALLPULS$dist <- numeric()
    file.create(file.out)
  }
  
  
  if (nrow(ALLPULS) >= 2) {
    
    ALLPULS$GENE <- ALLPULS$hmm
    ALLPULS$GENE[ALLPULS$sus!=""] <- ALLPULS$sus[ALLPULS$sus!=""]
    ALLPULS$GENE[ALLPULS$GENE==""] <- "unk"
    
    agg <- aggregate(ALLPULS$GENE, by=list(pulid=ALLPULS$pulid), paste, sep="-", collapse="-")
    uni <- unique(ALLPULS[,c("genome","pulid","contig")])
    start <- aggregate(ALLPULS$start,  by=list(pulid=ALLPULS$pulid), function(x) return(x[1]))
    end   <- aggregate(ALLPULS$end,  by=list(pulid=ALLPULS$pulid), function(x) return(x[length(x)]))
    
    pos <- merge(start, end, by="pulid")
    deets <- merge(uni, pos, by="pulid")
    
    out <- merge(deets,agg,by="pulid")
    
    colnames(out) <- c("pulid","genome","contigid","start","end","pattern")
    
    out <- out[,c("genome","pulid","contigid","start","end","pattern")]
    
    out <- out[order(out$contigid, out$start),]
    
    write.table(out, file.sum, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
  } else {
    file.create(file.sum)
  }
  
}

######################

list_dirs=list.dirs(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/",full.names = T,recursive = T)
muri_dirs=list_dirs[grep(pattern = "muribaculum_intestinale_",x=list_dirs)]
fl2=gsub("muribaculum_intestinale_", "muribaculum_intestinales_", muri_dirs)
file.rename(muri_dirs,fl2)

list_files=list.files(path = "S:/Abilash/Marbouty_annotation_refseq_genbank_2022/",pattern = "muribaculum_intestinale_",full.names = T,recursive = T)
fl2=gsub("muribaculum_intestinale_", "muribaculum_intestinales_", list_files)
file.rename(list_files,fl2)
