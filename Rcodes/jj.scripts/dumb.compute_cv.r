rm(list=ls())
debug <- FALSE
debug <- TRUE
if (debug){
    path <- '~/NCSU_Projects/Accord/Data/analysis/Subgroups/prelim_clustering/'#~/data/accord/data/analysis/GV/GV_allraces/'
    path <- '~/data/accord/data/analysis/GV/GV_allraces/'
    setwd('~/data/accord/data/analysis/')
    chr<-3
    chunk<-35
}else{
    args<-commandArgs(T)
    path<-args[1]
    chr<-args[2]
    chunk<-args[3]   
}


# Check if imputed genotype file exists
f_gz<-paste0("../imputation/outputs/chr",chr,".imputed.",chunk,".gz")
if (!file.exists(f_gz)) {
  cat("File",f_gz,"does not exist\n")
  quit()
}

# Parameters
info_threshold<-0.5
read_size<-1000

# Models to fit
models<-unlist(read.table(paste0(path,"/pheno_data/models.txt"),stringsAsFactors=F))
phenos<-sapply(strsplit(models,"~"),"[",1)

# Load phenotype
source("bin/load_pheno_data.r")
pheno_data<-load_pheno_data(path)

# read in model types
model.types <- unlist(read.table(paste0(path,"/modeltypes.txt"),stringsAsFactors=F))

# Get list of samples with no NA
retained_samples<-complete.cases(pheno_data)

# Initialize variables for HWE computation
n_obs<-sum(retained_samples)
O<-c(0,0,0)
E<-c(0,0,0)

#check for snps in list
snp.index <- NULL
if(file.exists(file.path(path,"snp_list.txt"))){
    print("You have a snp_list.txt file!!!")
  snp.list <- read.table(file.path(path,"snp_list.txt"),stringsAsFactors=F)
  system.time(file.snps<-read.table(pipe(paste0("cut -d ' ' -f2 ../imputation/outputs/chr",chr,".imputed.",chunk,"_info")),header=T,stringsAsFactors=F))
  snp.index <- which(file.snps[,1]%in%as.character(snp.list[,1]))
  if(length(snp.index)==0){
    stop(paste0("No SNP's in snp_list.txt were found in ../imputation/outputs/chr",chr,".imputed.",chunk,"_info"))
  } else{
      print("Found it!!!")
  }
    
}

# Load genotype info metric
geno_data_info<-read.table(pipe(paste0("cut -d ' ' -f5 ../imputation/outputs/chr",chr,".imputed.",chunk,"_info")),header=T)
geno_data_info<-geno_data_info>info_threshold
geno_data_i<-1

# Get index values for the 3 genotypes
geno_2_idx<-(1:length(retained_samples))*3
geno_1_idx<-geno_2_idx-1
geno_0_idx<-geno_2_idx-2

####################
# Initialize output files
####################
file_list<-list()
for (i in 1:length(phenos)) {
  file_list[[i]]<-file(paste0(path,"/oddsratio/chr",chr,".",chunk,".",phenos[i],".assoc"),open="w")
  cat("CHR\tSNP\tBP\tA1\tA2\tFRQ\tHWE_CHI2\tBETA\tSE\tP\n",file=file_list[[i]])
}

if(!is.null(snp.index)){
  read.range <- snp.index-1
  read_size <- 1
}else{
  read.range <- 0:ceiling(1e9/read_size)
}

read_i <- read.range[1]



for (read_i in read.range) {

  ####################
  # Load genotype data
  ####################
  obj<-try(geno_data<-read.table(gzfile(f_gz),header=F,stringsAsFactors=F,nrow=read_size,skip=read_i*read_size),silent=T)
  if (is(obj, "try-error")) break
  # Strip SNP info
  snp_data<-geno_data[,1:5]
  geno_data<-as.matrix(geno_data[,-c(1:5)])
  
  ####################
  # Fit model
####################
    if (debug){
        snp_i <- 1
    }
  for (snp_i in 1:dim(snp_data)[1]) {
    # Fit model if snp passed info threshold
    if (geno_data_info[geno_data_i]) {
      dosage<-geno_data[snp_i,geno_1_idx]+2*geno_data[snp_i,geno_2_idx]
      # Filter out samples with total posterior genotype probability < 0.1
      sum_geno<-geno_data[snp_i,geno_0_idx]+geno_data[snp_i,geno_1_idx]+geno_data[snp_i,geno_2_idx]
      dosage[sum_geno<.1]<-NA
      # Compute MAF
      maf<-mean(dosage[retained_samples],na.rm=T)/2
      # Make sure dosage reflects minor allele (e.g. dosage of 2 = 2 minor alleles)
      flip<-0
      if (maf>.5) {
        dosage<-2*geno_data[snp_i,geno_0_idx]+geno_data[snp_i,geno_1_idx]
        dosage[sum_geno<.1]<-NA
        maf<-mean(dosage[retained_samples],na.rm=T)/2
      } else {
        temp<-snp_data[snp_i,4]
        snp_data[snp_i,4]<-snp_data[snp_i,5]
        snp_data[snp_i,5]<-temp
      }
      # Compute HWE p-value
      d<-geno_data[snp_i,geno_2_idx]
      O[3]<-sum(d[retained_samples])
      d<-geno_data[snp_i,geno_1_idx]
      O[2]<-sum(d[retained_samples])
      d<-geno_data[snp_i,geno_0_idx]
      O[1]<-sum(d[retained_samples])
      E[3]<-n_obs*maf*maf
      E[2]<-n_obs*maf*(1-maf)*2
      E[1]<-n_obs*(1-maf)*(1-maf)
      if (flip) {
        O<-rev(O)
      }
      hwe_chi2<-round(sum((O-E)^2/E),digits=4)
      maf<-round(maf,digits=6)
      # Add dosage data to dataset
      test_data<-cbind(pheno_data,dosage)
      # Fit model(s)
      for (i in 1:length(phenos)) {
        f<-as.formula(paste0(models[i],"+dosage"))
        model.type <- model.types[i]
        if(model.type=="linear"){
          fit<-lm(f,test_data,na.action=na.omit)
          # Get model stats
          coef<-summary(fit)$coef
          if ("dosage"%in%row.names(coef)) {
            beta<-format(coef["dosage","Estimate"],digits=6)
            se<-format(coef["dosage","Std. Error"],digits=6)
            p<-format(coef["dosage","Pr(>|t|)"],digits=6)
          } else {
            beta<-NA
            se<-NA
            p<-NA
          }
        }else if(model.type=="logistic"){
          fit <- glm(f,test_data,family="binomial")
          # Get model stats
          coef<-summary(fit)$coef
          if ("dosage"%in%row.names(coef)) {
            beta<-format(coef["dosage","Estimate"],digits=6)
            se<-format(coef["dosage","Std. Error"],digits=6)
            p<-format(coef["dosage","Pr(>|z|)"],digits=6)
          } else {
            beta<-NA
            se<-NA
            p<-NA
          }
        }

        # Print output
        cat(paste0(chr,"\t",snp_data[snp_i,2],"\t",snp_data[snp_i,3],"\t",
                   snp_data[snp_i,4],"\t",snp_data[snp_i,5],"\t",maf,"\t",hwe_chi2,"\t",beta,"\t",se,"\t",p,"\n"),file=file_list[[i]])
        #print(summary(fit))
      }
    }
    # Increment info counter
    geno_data_i<-geno_data_i+1
  }
}

####################
# Close files
####################
for (i in 1:length(phenos)) {
  close(file_list[[i]])
}

warnings()
