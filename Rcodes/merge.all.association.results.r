## create file of all gwas results
rm(list=ls())
debug <-FALSE
#debug <- TRUE
#oneoff <- FALSE ## This is temporarily set to TRUE for the ICAPS study.
#oneoff <- FALSE

if (debug){
    #setwd('~/NCSU_Projects/Accord/accord_data/analysis/')
    #path <- '~/NCSU_Projects/Accord/accord_data/analysis/Sulfonylurea/sulf_90d_black2/'
    #pheno <- 'sulf_hba1c'
    #modeltype <- 'linear'
    setwd('~/data/accord/data/analysis/')
    path <- '~/data/accord/data/analysis/mindcog/mindcog_mri_gm/'
    pheno <- 'gm'
    modeltype <- 'linear'
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

head(imp.sub)
head(geno.file)
head(meta.file)
tail(imp.file[imp.file$CHR==4])
out <- merge(imp.sub,geno.file,all=TRUE,by=c("SNP","CHR","A1","A2"))
if (nrow(meta.file)>1){
out <- merge(out,meta.file,all=TRUE,by=c("SNP","CHR","A1","A2"))
}
out <- merge(imp.file,out,all=TRUE,by=c("SNP","CHR"))

# cols.rem <- c(-12,-13,-22,-23,-32,-33,-34,-35)#c("CHR.y","BP.y","CHR","BP")
# out <- out[,cols.rem,with=FALSE]
colnames(out)



#remove SNPs with duplicated SNP IDs
dup.snp <- unique(out$SNP[duplicated(out$SNP)])
dup.snp <- dup.snp[which(dup.snp!='.')]
out <- out[!out$SNP%in%dup.snp,]

#head(out)

#head(out$P.y.)

#unique(colnames(out))

##set column headers
#if (oneoff){
    ## Need to calculate OR from BETA (because I do not have time to recalculate the associations)
#    out$BETA <- exp(out$BETA)
#    out$SE <- exp(out$SE)
#} ## The BETA is now an OR and it will be remedied in the next lines of code:

if (modeltype == 'logistic'){
  if(length(which(colnames(imp.file)=="N"))>0){
    new.cols <- c("SNP","CHR","BP","A1_allImp","A2_allImp","FRQ_allImp","HWE_CHI2_allImp","OR_allImp",
                  "SE_allImp","P_allImp", "N_allImp",
                  "A1","A2","BP_subImp","FRQ_subImp","HWE_CHI2_subImp", "OR_subImp","SE_subImp",
                  "P_subImp", "N_subImp","BP_geno","TEST_geno","NMISS_geno","OR_geno","SE_geno",
                  "L95_geno","U95_geno","STAT_geno","P_geno","BP_meta","N_meta","P_meta","P_rndm_meta",
                  "OR_meta","OR_rndm_meta","Q_meta","I_meta")
  }else{
    new.cols <- c("SNP","CHR","BP","A1_allImp","A2_allImp","FRQ_allImp","HWE_CHI2_allImp","OR_allImp",
                  "SE_allImp","P_allImp", #"N_allImp",
                  "A1","A2","BP_subImp","FRQ_subImp","HWE_CHI2_subImp", "OR_subImp","SE_subImp",
                  "P_subImp", "N_subImp","BP_geno","TEST_geno","NMISS_geno","OR_geno","SE_geno",
                  "L95_geno","U95_geno","STAT_geno","P_geno","BP_meta","N_meta","P_meta","P_rndm_meta",
                  "OR_meta","OR_rndm_meta","Q_meta","I_meta")
  }  
  
}else{
  if(length(which(colnames(imp.file)=="N"))>0){
    new.cols <- c("SNP","CHR","BP","A1_allImp","A2_allImp","FRQ_allImp","HWE_CHI2_allImp","BETA_allImp",
                  "SE_allImp","P_allImp", "N_allImp",
                  "A1","A2","BP_subImp","FRQ_subImp","HWE_CHI2_subImp", "BETA_subImp","SE_subImp",
                  "P_subImp", "N_subImp","BP_geno","TEST_geno","NMISS_geno","BETA_geno",
                  "STAT_geno","P_geno","SE_geno","BP_meta","N_meta","P_meta","P_rndm_meta",
                  "OR_meta","OR_rndm_meta","Q_meta","I_meta")
    }else{
      new.cols <- c("SNP","CHR","BP","A1_allImp","A2_allImp","FRQ_allImp","HWE_CHI2_allImp","BETA_allImp",
                    "SE_allImp","P_allImp", #"N_allImp",
                    "A1","A2","BP_subImp","FRQ_subImp","HWE_CHI2_subImp", "BETA_subImp","SE_subImp",
                    "P_subImp", "N_subImp","BP_geno","TEST_geno","NMISS_geno","BETA_geno",
                    "STAT_geno","P_geno","SE_geno","BP_meta","N_meta","P_meta","P_rndm_meta",
                    "OR_meta","OR_rndm_meta","Q_meta","I_meta")
  }

}


colnames(out)<- new.cols

#backup <- out

## Create columns for ultimate p-value and direction
out$P_unadj <- 99
out$Type <- 'UNKNOWN'
out$Direction <- '??' 
#which(colnames(out) %in% 'Direction')
gc()

summary(out$P_rndm_meta)
summary(out$OR_meta)

## If a meta analysis P-value exists, store it
#dim(out[!is.na(P_rndm_meta)])
#head(out[!is.na(P_rndm_meta)])
#head(out[!is.na(P_rndm_meta),P_unadj])
out[!is.na(P_rndm_meta),"P_unadj"] <- out[!is.na(P_rndm_meta),"P_rndm_meta"]
out[!is.na(P_rndm_meta),"Type"] <- rep('META',nrow(out[!is.na(P_rndm_meta),"P_unadj"]))
#out[!is.na(P_rndm_meta),"Direction"] <- rep('META',nrow(out[!is.na(P_rndm_meta),"P_unadj"]))
## FIX ME ADD DIRECTION HERE ##
#head(out[!is.na(P_rndm_meta)])
gc()

## If a geno analysis P-value exists (and meta analysis does not exist), store it
#dim(out[!is.na(P_rndm_meta)])
#head(out[is.na(P_rndm_meta)&!is.na(P_geno)])
#head(out[!is.na(P_)])
#head(out[is.na(P_rndm_meta)&!is.na(P_geno),P_unadj])
out[is.na(P_rndm_meta)&!is.na(P_geno),"P_unadj"] <- out[is.na(P_rndm_meta)&!is.na(P_geno),"P_geno"]
out[is.na(P_rndm_meta)&!is.na(P_geno),"Type"] <- rep('GENO',nrow(out[is.na(P_rndm_meta)&!is.na(P_geno),"P_unadj"]))
#head(out[is.na(P_rndm_meta)&!is.na(P_geno)])
gc()

## Finally, if a allImp analysis P-value exists (and meta analysis AND P_geno do not exist), store it
#dim(out[!is.na(P_rndm_meta)])
#head(out[!is.na(P_allImp)&is.na(P_rndm_meta)&is.na(P_geno)])
#head(out[!is.na(P_)])
#head(out[!is.na(P_allImp)&is.na(P_rndm_meta)&is.na(P_geno),P_adj])
out[!is.na(P_allImp)&is.na(P_rndm_meta)&is.na(P_geno),"P_unadj"] <- out[!is.na(P_allImp)&is.na(P_rndm_meta)&is.na(P_geno),"P_allImp"]
out[!is.na(P_allImp)&is.na(P_rndm_meta)&is.na(P_geno),"Type"] <- rep('IMPU',nrow(out[!is.na(P_allImp)&is.na(P_rndm_meta)&is.na(P_geno),"P_unadj"]))
#head(out[!is.na(P_allImp)&is.na(P_rndm_meta)&is.na(P_geno)])
gc()

## Overwrite the allImp BP whenever there is a genotype BP value
out[!is.na(BP_geno),"BP"] <- out[!is.na(BP_geno),"BP_geno"]
#length(which(is.na(out[,BP])))
gc()

#################################################################


## Drop probes without valid P-values in any case (geno, allImp, meta)
#dim(out[is.na(P_allImp)&is.na(P_rndm_meta)&is.na(P_geno)])
## 160168 x 40
#dim(out[P_unadj==99,"P_unadj"])
#head(out[P_unadj==99,"P_unadj"])
out[P_unadj==99,"P_unadj"] <- NA
#dim(out[is.na(P_unadj)])
## There are 160168 probes with no valid pval
#head(out[is.na(P_unadj)])
#tmp <- dim(out)
## Drop NAs
out <- out[!is.na(P_unadj),]
#dim(out)-tmp
gc()

# #########################
# # Compute genomic control
# #########################
# gc_probes<-unlist(read.table("/home/accord/data/geno_data/LD.prune.in",header=F,stringsAsFactors=F))
# #gc_probes<-unlist(read.table("~/NCSU_Projects/Accord/Data/LD.prune.in",header=F,stringsAsFactors=F))
# #head(gc_probes)
# #length(gc_probes)
# 
# #gc_idx<-chip$SNP%in%gc_probes
# #length(data[Type=='GENO',SNP])
# #length(data[,SNP])
# gc_idx<-which(out[,SNP]%in%gc_probes)
# #length(gc_idx)
# #gc_idx<-which(out[,SNP]%in%gc_probes)
# gc_p<-out[gc_idx,P_geno]
# #length(gc_p)
# #head(gc_p)
# #summary(gc_p)
# chi2_stats<-qchisq(na.omit(gc_p),df=1,lower.tail=F)
# #?qchisq
# #head(chi2_stats)
# gc <- median(chi2_stats)/qchisq(0.5,df=1)
# #gc
# 
# ## Output info
# output<-file(paste0(path,"/outputs/",pheno,"_info.txt"),open="w")
# cat("Phenotype:         ",pheno,"\n",
#     #"MAF cutoff:        ",maf_cutoff,"\n",
#     #"Chromosomes:       ",paste(chr,collapse=" "),"\n",
#     "Genomic inflation: ",format(gc,digits=4),"\n",
#     "Genotyped:         ",nrow(out[Type=='GENO',]),"\n",
#     "Imputed (all):           ",nrow(out[Type=='IMPU',]),"\n",
#     "Meta:           ",nrow(out[Type=='META',]),"\n",
#     sep="",file=output)
# close(output)
# gc()
#########################
# Adjust p-values
#########################
# chi2_stats<-qchisq(out[,P_unadj],df=1,lower.tail=F)
# chi2_stats<-chi2_stats/gc
# adjP<-pchisq(chi2_stats,df=1,lower.tail=F)
# out$P_adj <- adjP
# gc()

##summary(data$P_adj)




#add minor allele frequency column
frq<-read.table(file.path(path,"association_cv","plink.frq"),header=T,stringsAsFactors=F)
MAF <- out$FRQ_allImp
idx <-match(frq$SNP,out$SNP)
MAF[idx[which(!is.na(idx))]] <- frq[which(!is.na(idx)),"MAF"]
out <- cbind(out,MAF)

write.table(out,file=file.path(path,"association_cv",paste0("all_commonVariantAnalysis_",pheno,"_results.tsv")),row.names=FALSE)
cat("***FINISHED***/n")
#incorporate into pipeline before plink meta
# test13 <- imp.sub[,c("CHR","BP","SNP","BETA","SE","P","A1","A2"),with=FALSE]
# write.table(test13,file=file.path(path,"association_cv",paste0("test.assoc")),row.names=FALSE,sep=" ",quote=FALSE)
# 
# SE <- geno.file$BETA/geno.file$STAT
# geno.test <- cbind(geno.file,SE)
# write.table(geno.test,file=file.path(path,"association_cv",paste0("geno.assoc")),row.names=FALSE,sep=" ",quote=FALSE)
# 
# 


