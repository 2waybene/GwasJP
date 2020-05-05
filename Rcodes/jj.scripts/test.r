## create file of all gwas results
rm(list=ls())
debug <- FALSE
debug <- TRUE
#oneoff <- FALSE ## This is temporarily set to TRUE for the ICAPS study.
oneoff <- TRUE

if (debug){
    setwd('~/data/accord/data/analysis/')
    path <- '~/data/accord/data/analysis/icaps/icaps.thiazide.all_races/'
    pheno <- 'ICAPS_po'
    modeltype <- 'logistic'
}else{
    args<-commandArgs(T)
    path<-args[1]
    pheno <- args[2]
    modeltype <- args[3]
}

#chr<-args[2]
#chunk<-args[3]
library(data.table)##install.packages('data.tables',type="source",repos = "http://Rdatatable.github.io/data.table")

meta.file <- fread(file.path(path,"association_cv",paste0("plink_meta_",pheno,".meta")),sep=" ",header=TRUE)
#meta.file <- read.table(file.path(path,"association_cv",paste0("plink_meta_",pheno,".meta")),sep="\t",header=FALSE)
setkey(meta.file,SNP)
geno.file <- fread(file.path(path,"association_cv",paste0("chr0.",pheno,".assoc.",modeltype)),sep=" ",header=TRUE)
setkey(geno.file,SNP)
imp.file <- fread(file.path(path,"association_cv",paste0("allChrImputed.",pheno,".assoc")),sep=" ",header=TRUE)
setkey(imp.file,SNP)
imp.sub <- fread(file.path(path,"association_cv",paste0("allChrImputed_forMetaAnalysis.",pheno,".assoc")),sep="\t",header=TRUE)
setkey(imp.sub,SNP)



head(geno.file)
summary(geno.file$P)

head(imp.sub)
summary(imp.sub$P)

hist(imp.sub$P)

hist(imp.sub$P[which(imp.sub$P<.00001)])


length(imp.sub$P[which(imp.sub$P>.000000001)])

