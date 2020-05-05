##########################
# Get GWAS Peaks
# Author: Jon Leirer
# Date: 10/11/2019
# Description: Copied from ./bin/get_peaks.r
#   Made changes to ensure adjusted is interpreted as logical


rm(list=ls())

library(BiocGenerics)
library(IRanges)
library(XVector)
library(GenomicRanges)
library(data.table)

debug <- FALSE
#debug <- TRUE
if (debug){
    path <- '~/data/accord/data/analysis/Metformin/met_90d_allraces4'
    maf <- .03
    adjusted <- FALSE
}else{
    args<-commandArgs(T)
    path<-args[1]
    maf<-as.numeric(args[2])
    adjusted<-as.logical(args[3])
}

print(args)

##maf <- .3

p_value_cutoff<-1e-5
delta<-75000




# source("https://bioconductor.org/biocLite.R")
# biocLite()
# biocLite("S4Vectors")
# biocLite("IRanges")
# biocLite("XVector")
# biocLite("GenomicRanges")
# install.packages("/home/drotroff/R/BiocGenerics_0.10.0.tar.gz",repos=NULL,type="source")
# install.packages('S4Vectors')
# install.packages("/home/drotroff/R/IRanges_2.0.1.tar.gz",repos=NULL,type="source")
# suppressMessages(library(XVector,lib="/home/drotroff/R/library/"))
# suppressMessages(library(GenomicRanges,lib="/home/drotroff/R/library/"))



if(adjusted){
    cat('Loading adjusted p-vals...\n')
}else{
   cat('Loading unadjusted p-vals...\n')
}


#length(is.na(data[,P_unadj]))

phenotypes<-unlist(read.table(paste0(path,"/phenotypes.txt"),stringsAsFactors=F))
#pheno_i <- 1
for (pheno_i in 1:length(phenotypes)) {
    ## Load data
    data <- fread(file.path(path,"association_cv",paste0("allChrom_",phenotypes[pheno_i],"_gc_adj.csv")))
    if (adjusted){
        loc<-data$BP[data$P_adj<=p_value_cutoff]
    }else{
        loc<-data$BP[data$P_unadj<=p_value_cutoff]        
    }
    #summary(loc)
    #head(loc)
    #head(data[data$P_unadj<=p_value_cutoff,])
    if (length(loc)>0){
        loc_delta<-cbind(loc-delta,loc+delta)
        gr<-GRanges(seqnames="*",IRanges(loc_delta[,1],loc_delta[,2]),strand="*") 
        rr<-as.data.frame(reduce(gr))
        out<-c()
        for (rr_i in 1:dim(rr)[1]) {
            out<-rbind(out,data[data$BP>=rr$start[rr_i] & data$BP<=rr$end[rr_i],])
        }
        out<-out[order(out$BP),]        
        write.table(out,paste0(path,"/peak_data/",phenotypes[pheno_i],"_gc_adj.assoc.linear"),row.names=F,quote=F,sep="\t")
    }
}
