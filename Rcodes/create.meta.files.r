rm(list=ls())
debug <- FALSE
#debug <- TRUE
if (debug){
    path <- '~/data/accord/data/analysis/icaps/icaps.thiazide.all_races/'
    setwd('~/data/accord/data/analysis/')
    pheno <- 'ICAPS_po'
    modeltype <- 'logistic'
}else{
    args<-(commandArgs(TRUE))
    path<-args[1]
    pheno <- args[2]
    modeltype <- args[3]
}

cat('Creating meta analysis input files for phenotype ',pheno,' and modeltype ',modeltype,'...\n')

library(data.table)
##ADD A2 to Genotyped file
geno.file.name <- list.files(file.path(path,"association_cv"),pattern = "chr0")
geno.file.name <- geno.file.name[grep(pheno,geno.file.name)]
geno.file.name <- geno.file.name[grep("assoc",geno.file.name)]
cat("Loading file: ",geno.file.name,"\n")
geno.file <- fread(file.path(path,"association_cv",geno.file.name),header=TRUE,nrows = 3)
#head(geno.file)

if(!any(grepl("A2",colnames(geno.file))) || !any(grepl("SE",colnames(geno.file)))){
  geno.file <- fread(file.path(path,"association_cv",geno.file.name),header=TRUE)
  #print('hello1')
  if(!any(grepl("A2",colnames(geno.file)))){
    bim.file <- read.table("/home/accord/data/geno_data/post_qc.unc.uva.merged.bim",header=FALSE)
    #bim.file <- read.table("~/NCSU_Projects/post_qc.unc.uva.merged.bim",header=FALSE)
    A2 <- bim.file[match(geno.file$SNP,bim.file[,2]),6]
    geno.file <- cbind(geno.file[,1:4],A2,geno.file[,5:ncol(geno.file)])
  }
  #print('hello2')
  if(!any(grepl("SE",colnames(geno.file)))){
      #print('SE1')
      if (modeltype == 'linear'){
          SE <- geno.file$BETA/geno.file$STAT
      }else{
	  cat('Error: Double check your SE. This should not run for ICAPS. ',file.path(path,"association_cv",geno.file.name),'\n')##FIX
          #SE <- exp(log(geno.file$OR)/geno.file$STAT)
      }
    #length(SE)
    #print('SE2')
    geno.file <- cbind(geno.file,SE)
  }
  #print('hello3')
  write.table(geno.file,file.path(path,"association_cv",geno.file.name),row.names=FALSE,quote=FALSE)
  cat("Done\n")
}else{
  cat("A2 column and SE column already exist in file...skipping\n")
}

if (debug){
    head(geno.file)
}
## CREATE ALL IMPUTED FILE FOR METAL


imp.files <- list.files(file.path(path,"association_cv","imputed_chunks","imputed_chunks_forMeta"),pattern = "assoc")
#imp.files <- imp.files[grep("chr0",imp.files,invert=TRUE)]

out.name <- paste0("allChrImputed_forMetaAnalysis.",pheno,".assoc")
cat('Merging imputed_forMeta association results to create: ',out.name,'\n')
#if(file.exists(file.path(path,"association_cv",out.name))){
#  file.remove(file.path(path,"association_cv",out.name))
#}
pheno.imp.files <- imp.files[grep(pheno,imp.files)]
#pheno.imp.files <- pheno.imp.files[grep("allChr",pheno.imp.files,invert=TRUE)]
output.file <- NULL
for(f in pheno.imp.files){
  cat("Merging:",f," file number ",which(pheno.imp.files==f)," of ",length(pheno.imp.files),"\n")
  temp.file <- fread(file.path(path,"association_cv","imputed_chunks","imputed_chunks_forMeta",f),header=TRUE)
  output.file <- rbind(output.file,temp.file)
}
write.table(output.file,file.path(path,"association_cv",out.name),row.names=FALSE,quote=FALSE,sep="\t")

cat("Done\n")

