rm(list=ls())
library(BiocGenerics)
library(IRanges)
library(XVector)
library(GenomicRanges)
library(data.table)

debug <- FALSE
debug <- TRUE
if (debug){
    path <- '~/data/accord/data/analysis/Sulfonylurea/sulf_90d_black2'
    maf <- .03
    adjusted <- TRUE
}else{
    args<-commandArgs(T)
    path<-args[1]
    maf<-as.numeric(args[2])
    adjusted<-args[3]
}


##maf <- .3

p_value_cutoff<-1e-5
delta<-100000




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
peaksnps <- c()
#pheno_i <- 1
for (pheno_i in 1:length(phenotypes)) {
    ## Load data
    data <- fread(file.path(path,"association_cv",paste0("allChrom_",phenotypes[pheno_i],"_gc_adj.csv")))
    ## Temporarily rename all '.' SNPs
    data$SNP[data$SNP=='.'] <- paste0('Unknown',seq(1,length(data$SNP[data$SNP=='.'])))
    if (adjusted){
        loc<-data$BP[data$P_adj<=p_value_cutoff]
        rs.loc <- data$SNP[data$P_adj<=p_value_cutoff]
    }else{
        loc<-data$BP[data$P_unadj<=p_value_cutoff]        
        rs.loc <- data$SNP[data$P_unadj<=p_value_cutoff]
    }
    #head(data[data$SNP=='.',])
    #dim(data[data$SNP=='.',])
    #head(paste0('Unknown',seq(1,length(data$SNP[data$SNP=='.']))))
    

    ## DEBUG
    #summary(loc)
    #head(loc)
    #head(rs.loc)
    #head(data[data$P_unadj<=p_value_cutoff,])
    #data[data$P_adj==min(data$P_adj),]
    #which(rs.loc %in% data$SNP[data$P_adj==min(data$P_adj)])
    #loc<-loc[c(1,360,361)]
    #rs.loc<-rs.loc[c(1,360,361)]
    #data[data$SNP %in% rs.loc]
    ## DEBUG
    if (length(loc)>0){
        loc_delta<-cbind(loc-delta,loc+delta)
        ## DEBUG        
        #loc_delta<-cbind(loc[1:3]-delta,loc[1:3]+delta)
        ## DEBUG
        gr<-GRanges(seqnames="*",IRanges(loc_delta[,1],loc_delta[,2]),strand="*") 
        rr<-as.data.frame(reduce(gr))
        #dim(rr)
        #head(rr)        
        out<-c()
        for (rr_i in 1:dim(rr)[1]) {
            out<-rbind(out,data[data$BP>=rr$start[rr_i] & data$BP<=rr$end[rr_i],])
        }
        out<-out[order(out$BP),]        
        #dim(out)
        #head(out)        
        #summary(out$BP)
        wout <- out
        wout$SNP[grep('Unknown',out$SNP)] <- '.'
        write.table(wout,paste0(path,"/peak_data/",phenotypes[pheno_i],"_gc_adj.assoc.linear"),row.names=F,quote=F,sep="\t")
        if (adjusted){
            wout <- data.frame(MarkerName=wout$SNP,'P-value'=wout$P_adj)
            #head(wout)
            colnames(wout)[2] <- "P-value"
            write.table(wout,paste0(path,"/peak_data/",phenotypes[pheno_i],"_lz.batchmode_metal"),row.names=F,quote=F,sep="\t")
        }else{
            wout <- data.frame(MarkerName=wout$SNP,'P-value'=wout$P_unadj)
            #head(wout)
            colnames(wout)[2] <- "P-value"            
            write.table(wout,paste0(path,"/peak_data/",phenotypes[pheno_i],"_lz.batchmode_metal"),row.names=F,quote=F,sep="\t")
        }
        options(scipen = 999)
        ## Make file for # of plots (basically, contains the center for each plot (peak))
        peak.range <- delta # what's the range of a skyscraper
        peak.rsids <- data.frame(snp='rsNA',chr=99,start=0,end=0,flank=delta,run='yes',m2zargs='showAnnot=T',stringsAsFactors=FALSE)
        peak.rsids <- peak.rsids[-1,]
        #rstmp <- rs.loc[2]
        #for (tmpi in rstmp)
        counter <- 1
        for(rstmp in rs.loc){
            tmp.data <- out[#intersect(
                                intersect(which(out$CHR == out$CHR[which(out$SNP==rstmp)]),
                                    which(out$BP<=out$BP[out$SNP==rstmp]+peak.range)),]
                            #    which(out$BP>=out$BP[out$SNP==rstmp]-peak.range)),]
            dim(tmp.data)
            head(tmp.data)
            dim(out)
            #tmp.data
            ## DEBUG ADD CHECK FO UNADJ vs ADJ)
            ## Check if we are currently on min SNP in peak. Save if we are ##
            tmp.data[which(tmp.data$P_adj==min(tmp.data$P_adj)),]
            if (rstmp == tmp.data$SNP[which(tmp.data$P_adj==min(tmp.data$P_adj))][1]){
                peak.rsids[counter,] <- c(rstmp,
                                            tmp.data$CHR[which(tmp.data$SNP %in% rstmp)],
                                            NA,
                                            NA,
                                            delta,
                                            'yes',
                                            'showAnnot=T')
                counter <- counter + 1
            }
        }
        peak.rsids$snp[grep('Unknown',peak.rsids$snp)] <- '.'
        write.table(peak.rsids,paste0(path,"/peak_data/",phenotypes[pheno_i],"_lz.batchmode_hitspec"),row.names=F,quote=F,sep="\t")  
    }    
}


