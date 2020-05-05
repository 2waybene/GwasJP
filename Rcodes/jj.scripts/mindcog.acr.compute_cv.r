## This file creates the two interaction models (APOE and ACR) requested
## by Jasmin Divers for the MINDCOG consortium 

rm(list=ls())

require(data.table)
debug <- FALSE
#debug <- TRUE
if (debug){
    path <- '~/data/accord/data/analysis/mindcog/mindcog_cog_dsst'
    setwd('~/data/accord/data/analysis')
    phename = 'dsst'
    chr<-21
    chunk<-4
    f_gz <- '~/data/accord/data/geno_data/mindcog.apoe.chunk/chr21.imputed.4'
    #read_i = 26969
    #snp_i = 121
    #geno_data_i = 921
}else{
  ## DO NOT ADD DEBUG INFO IN HERE
    args <- commandArgs(T)
     path <- args[1]
     chr <- args[2]
     chunk <- args[3]
     f_gz <- args[4]    
}

load_pheno_data<-function(path,phenopath) {
  data<-read.table(phenopath,header=T,stringsAsFactors=F)
  names(data)[1]="FID"
  names(data)[2]="IID"
  #head(data)  
  #pull phenotypes so they don't get centered/scaled #this is what is new from old version
  phenos <- unlist(read.table(file.path(path,"phenotypes.txt"),stringsAsFactors=FALSE))

  # Scale continuous variables
  data_class<-sapply(data,class)
  data_levels<-sapply(data,function(x) length(levels(factor(x))))
  scale_idx<-data_class=="numeric" | (data_class=="integer" & data_levels>8) & c(0,0,rep(1,dim(data)[2]-2))
  scale_idx[which(names(data_class)%in%phenos)] <- FALSE
  std<-function(x){(x-mean(x,na.rm=T))/sd(x,na.rm=T)}
  data[,scale_idx]<-data.frame(sapply(data[,scale_idx],std))

  # Define categorical variables
  vars<-names(data)[data_class=="integer" & data_levels<=8]
  if(length(grep("prt_hba1c_timegap",vars))>0) vars <- vars[which(vars!="prt_hba1c_timegap")]
  for (i in 1:length(vars)) {data[,vars[i]]<-factor(data[,vars[i]])}

  # Create full list for all samples with imputed genotypes
  samples<-read.table(pipe(paste0("cut -d ' ' -f1,2 ../imputation/data/chr1.unphased.fam")),header=F,stringsAsFactors=F)
  samples$order<-seq(len=nrow(samples))
  names(samples)<-c("FID","IID","order")
  pheno_data<-merge(samples,data,sort=F,all.x=T)
  pheno_data<-pheno_data[sort.list(pheno_data$order), -3]

  if(any(grepl("AX",colnames(pheno_data)))){
      class(pheno_data[,grep("AX",colnames(pheno_data))]) <- "numeric"
  }
  #head(pheno_data)
  #head(samples)
  return(pheno_data)
}


## debug
cat(path)

## This info comes from the preceding bash script now CAN DELETE
# Check if imputed genotype file exists
#f_gz<-paste0("../imputation/outputs/chr",chr,".imputed.",chunk,".gz")
#if (!file.exists(f_gz)) {
#  cat("File",f_gz,"does not exist\n")
#  quit()
#}

# if (debug){
#     f_gz<-paste0("../imputation/outputs/test/chr",chr,".imputed.",chunk,"")
#     }

# Parameters
info_threshold<-0.5
read_size<-1000

#cat('hello')

# Models to fit
#models<-unlist(read.table(paste0(path,"/pheno_data/models.txt"),stringsAsFactors=F))
models<-unlist(read.table(paste0(path,"/acr/models.txt"),stringsAsFactors=F))
phenos<-sapply(strsplit(models,"~"),"[",1)

# Load phenotype
#source("bin/load_pheno_data.r")
pheno_data<-load_pheno_data(path,paste0(path,"/acr/pheno_data.acr.txt"))
head(pheno_data)

# read in model types
model.types <- unlist(read.table(paste0(path,"/modeltypes.txt"),stringsAsFactors=F))

# Get list of samples with no NA
retained_samples<-complete.cases(pheno_data)

head(pheno_data[retained_samples,])

# Initialize variables for HWE computation
n_obs<-sum(retained_samples)
O<-c(0,0,0)
E<-c(0,0,0)

#cat('hello2')
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

#cat('hello3')

####################
# Initialize output files
####################
file_list<-list()
for (i in 1:length(phenos)) {
  file_list[[i]]<-file(paste0(path,"/acr/association_cv/imputed_chunks/chr",chr,".",chunk,".",phenos[i],".assoc"),open="w")
  if(model.types[i]=="linear"){
      cat("CHR\tSNP\tBP\tA1\tA2\tTEST\tFRQ\tHWE_CHI2\tBETA\tSE\tP\tN\n",file=file_list[[i]])
  }else{ #### THIS IS NOT AN OR
     cat("CHR\tSNP\tBP\tA1\tA2\tTEST\tFRQ\tHWE_CHI2\tOR\tSE\tP\tN\n",file=file_list[[i]])
 }

}

if(!is.null(snp.index)){
  read.range <- snp.index-1
  read_size <- 1
}else{
  read.range <- 0:ceiling(1e9/read_size)
}

#head(read.range)
#length(read.range)
#read_i <- read.range[length(read.range)]
#read_i <- read.range[1]

#cat('hello4')

cat(f_gz)
beta = c()
se = c()
p = c()
n = c()
coeffs = c("dosage","log10ACR","log10ACR:dosage")
#tmp = 1
#read.range <- 1
for (read_i in read.range) {

  ####################
  # Load genotype data
  ####################
  #obj<-try(geno_data<-read.table(gzfile(f_gz),header=F,stringsAsFactors=F,nrow=read_size,skip=read_i*read_size),silent=T)
    ##if (is(obj, "try-error")) break
    #cat('Do you even fread, bro?\n')    
    geno_data<-try(fread(f_gz,header=F,verbose=F,showProgress=F,nrow=read_size,skip=read_i*read_size))
    if (any(class(geno_data)=="try-error")){
	    cat(paste0('try-error: ',f_gz,'\n'))
	    break
    }
    #cat('Do you?\n')
    #geno_data<-try(fread(?gzfile(f_gz),header=F,nrows=read_size,skip=read_i*read_size),silent=T)
    #geno_data<-fread(f_gz,header=F,showProgress=T)
    #geno_data<-fread(gzfile(f_gz),header=F,showProgress=T)

    #geno_data<-fread(f_gz,header=F,nrows=read_size,skip=read_i*read_size)
    ##geno_data<-try(fread(f_gz,header=F),silent=T)
  head(geno_data)[,1:10]
  #cat('Hello world')
  
  #convert to data frame so that you can use the same code as before
  geno_data<- as.data.frame(geno_data)
  # Strip SNP info
  snp_data<-geno_data[,1:5]
  geno_data<-as.matrix(geno_data[,-c(1:5)])
  #cat(dim(geno_data),'\n')
  ####################
  # Fit model
  ####################
  #if (debug){
  #    snp_i <- 1
  #    geno_data[snp_i,1:20]
  #    snp_data[snp_i,]
  #    summary(geno_data[snp_i,])
  #    length(geno_data[snp_i,])/3
  #    summary(dosage)
  #    summary(sum_geno)
  #    head(geno_1_idx)
  #    head(geno_0_idx)
  #    head(geno_2_idx)
  #}
  #snp_i = 1

  for (snp_i in 1:dim(snp_data)[1]) {
  #for (snp_i in 1:2) {
 #cat('am I here?1')    
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
      ## A1 is now the minor allele, A2 is the major allele
      
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
      i=1
      for (i in 1:length(phenos)) {
	#cat('am I here?1')
        f<-as.formula(paste0(models[i],"+dosage+log10ACR*dosage"))
        model.type <- model.types[i]
        if(model.type=="linear"){
          # sub.data <- test_data[sample(1:nrow(test_data),nrow(test_data)/3),]
          # summary(lm(f,sub.data,na.action=na.omit))
          fit<-lm(f,test_data,na.action=na.omit)
          # Get model stats
          coef<-summary(fit)$coef
          for (coeffy in coeffs){
            if (coeffy%in%row.names(coef)) {
              beta<-c(beta,format(coef[coeffy,"Estimate"],digits=6))
              se<-c(se,format(coef[coeffy,"Std. Error"],digits=6))
              p<-c(p,format(coef[coeffy,"Pr(>|t|)"],digits=6))
              n<-c(n,length(fit$fitted.values))
            } else {
              beta<-c(beta,NA)
              se<-c(se,NA)
              p<-c(p,NA)
              n<-c(n,NA)
            }          
          }
        }else if(model.type=="logistic"){
	  #cat('am I here?')
          fit <- glm(f,test_data,family="binomial")
          # Get model stats
          coef<-summary(fit)$coef
          if ("dosage"%in%row.names(coef)) {
            beta<-format(exp(coef["dosage","Estimate"]),digits=6)
            se<-format(exp(coef["dosage","Std. Error"]),digits=6)
            p<-format(coef["dosage","Pr(>|z|)"],digits=6)
            n<- length(fit$fitted.values)
          } else {
            beta<-NA
            se<-NA
            p<-NA
            n<-NA
          }
        }
        
        for (c in 1:length(coeffs)){
          #print read_i
          # Print output
          cat(paste0(chr,"\t",snp_data[snp_i,2],"\t",snp_data[snp_i,3],"\t",
                   snp_data[snp_i,4],"\t",snp_data[snp_i,5],"\t",coeffs[c],"\t",maf[c],"\t",hwe_chi2[c],"\t",beta[c],"\t",se[c],"\t",p[c],"\t",n[c],"\n"),file=file_list[[i]])
        }
        beta = c()
        se = c()
        p = c()
        n = c()        
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
