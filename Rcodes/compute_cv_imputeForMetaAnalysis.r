rm(list=ls())
debug <- FALSE
#debug <- TRUE
if (debug){
    path <- '~/NCSU_Projects/Accord/accord_data/analysis/Metformin/met_90d_allraces4_noGlyArmCorrect'#'~/NCSU_Projects/Accord/accord_data/analysis/Metformin/met_90d_allraces4_intArm'
    setwd('~/NCSU_Projects/Accord/accord_data/analysis')#setwd('~/NCSU_Projects/Accord/accord_data/analysis')
    #path <- '~/data/accord/data/analysis/icaps/icaps.thiazide.all_races/'
    #setwd('~/data/accord/data/analysis/')
    chr <- 8
    chunk <- 14
}else{
    args <- commandArgs(T)
    path <- args[1]
    chr <- args[2]
    chunk <- args[3]
}



library(data.table)
# Check if imputed genotype file exists
f_gz<-paste0("../imputation_forMETA/chr",chr,".imputed.",chunk,".txt")
if (!file.exists(f_gz)) {
  cat("File",f_gz,"does not exist\n")
  quit()
}

#check for snps in list
if(file.exists(file.path(path,"snp_list.txt"))){
  print("You have a snp_list.txt file!!!")
  snp.list <- read.table(file.path(path,"snp_list.txt"),stringsAsFactors=F)
  system.time(file.snps<-read.table(pipe(paste0("cut -d ' ' -f2 ../imputation/outputs/chr",chr,".imputed.",chunk,"_info")),header=T,stringsAsFactors=F))
  snp.exists <- which(file.snps[,1]%in%as.character(snp.list[,1]))
  if(length(snp.exists)==0){
    stop(paste0("No SNP's in snp_list.txt were found in ../imputation/outputs/chr",chr,".imputed.",chunk,"_info"))
  } else{
    print("Found it!!!")
  }
}

# Models to fit
models<-unlist(read.table(paste0(path,"/pheno_data/models.txt"),stringsAsFactors=F))
phenos<-sapply(strsplit(models,"~"),"[",1)

# Read in bim file to overwrite A1 and A2 for META matching...
bim.file <- fread("/home/accord/data/geno_data/post_qc.unc.uva.merged.bim",header=FALSE)

#Read subjects only in UNC set for imputed input to META
#note: FID and IID are switched
unc_only_subj <- read.table("/home/accord/data/geno_data/files.for.meta.analysis/IIDs.unc.only.no.uva.data.txt",header=FALSE)

# Load phenotype
source("bin/load_pheno_data.r")
pheno_data<-load_pheno_data(path)

exclude_subjects <- which(!pheno_data[,"FID"]%in%as.character(unc_only_subj[,1]))

# read in model types
model.types <- unlist(read.table(paste0(path,"/modeltypes.txt"),stringsAsFactors=F))

## Get list of samples with no NA
retained_samples<-complete.cases(pheno_data)
retained_samples[exclude_subjects] <- FALSE

# Initialize variables for HWE computation
n_obs<-sum(retained_samples)
O<-c(0,0,0)
E<-c(0,0,0)




# Get index values for the 3 genotypes
geno_2_idx<-(1:length(retained_samples))*3
geno_1_idx<-geno_2_idx-1
geno_0_idx<-geno_2_idx-2


  ####################
  # Load genotype data
  ####################
  #obj<-try(geno_data<-read.table(gzfile(f_gz),header=F,stringsAsFactors=F,nrow=read_size,skip=read_i*read_size),silent=T)
  geno_data<-try(fread(f_gz,header=F),silent=T)
  if (any(class(geno_data)%in%"try-error")) stop("error on fread...most likely empty")
  
  #convert to data frame so that you can use the same code as before
  geno_data<- as.data.frame(geno_data)
  
  #check for snps in list
  #snp.index <- NULL
  if(file.exists(file.path(path,"snp_list.txt"))){
    print("You have a snp_list.txt file!!!")
    snp.list <- read.table(file.path(path,"snp_list.txt"),stringsAsFactors=F)
    snp.index <- which(geno_data[,2]%in%as.character(snp.list[,1]))
    if(length(snp.index)==0){
      stop(paste0("No SNP's in snp_list.txt were found in ",f_gz))
    } else{
      print("Found it!!!")
      geno_data <- geno_data[snp.index,]
    }
    
  }

  
  ####################
  # Initialize output files
  ####################
  file_list<-list()
  for (i in 1:length(phenos)){
    file_list[[i]]<-file(paste0(path,"/association_cv/imputed_chunks/imputed_chunks_forMeta/chr",chr,".",chunk,".",phenos[i],".assoc"),open="w")
    if(model.types[i]=="linear"){
      cat("CHR\tSNP\tBP\tA1\tA2\tFRQ\tHWE_CHI2\tBETA\tSE\tP\tN\n",file=file_list[[i]])
    }else{
      cat("CHR\tSNP\tBP\tA1\tA2\tFRQ\tHWE_CHI2\tOR\tSE\tP\tN\n",file=file_list[[i]])
    }
  }
  
  # Strip SNP info
  snp_data<-geno_data[,1:5]
  geno_data<-as.matrix(geno_data[,-c(1:5)])

  
    
  ####################
  # Fit model
  ####################
  #snp_i <- 1
  for (snp_i in 1:nrow(snp_data)) {
    # Fit model if snp passed info threshold
    rs.id <- snp_data[snp_i,2]
      dosage<-geno_data[snp_i,geno_1_idx]+2*geno_data[snp_i,geno_2_idx]
      # Filter out samples with total posterior genotype probability < 0.1
      sum_geno<-geno_data[snp_i,geno_0_idx]+geno_data[snp_i,geno_1_idx]+geno_data[snp_i,geno_2_idx]
      dosage[sum_geno<.1]<-NA
      # Compute MAF
      maf<-mean(dosage[retained_samples],na.rm=T)/2
      # Make sure dosage reflects minor allele (e.g. dosage of 2 = 2 minor alleles)
      
      ## REASSIGN ALLELES FOR META ANALYSIS-MUST MATCH UVA STRAND
      bim.vals <- bim.file[V2==rs.id]
      A1 <- bim.vals[,V5]
      A2 <- bim.vals[,V6]
      
      flip<-0
      if (maf>.5) {
        dosage<-2*geno_data[snp_i,geno_0_idx]+geno_data[snp_i,geno_1_idx]
        dosage[sum_geno<.1]<-NA
        maf<-mean(dosage[retained_samples],na.rm=T)/2
        flip<-1
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
      if (flip==1) {
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
          fit<-lm(f,test_data[retained_samples,],na.action=na.omit)
        # Get model stats
        coef<-summary(fit)$coef
        if ("dosage"%in%row.names(coef)) {
          beta<-format(coef["dosage","Estimate"],digits=6)
          se<-format(coef["dosage","Std. Error"],digits=6)
          p<-format(coef["dosage","Pr(>|t|)"],digits=6)
          n<- length(fit$fitted.values)
        } else {
          beta<-NA
          se<-NA
          p<-NA
          n<- NA
        }
        }else if(model.type=="logistic"){
          fit <- glm(f,test_data[retained_samples,],family="binomial")
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
            n<- NA
          }
        }
        # Print output
        cat(paste0(chr,"\t",snp_data[snp_i,2],"\t",snp_data[snp_i,3],"\t",
                   A1,"\t", #A1
                   A2,"\t", #A2
                   maf,"\t",hwe_chi2,"\t",beta,"\t",se,"\t",p,"\t",n,"\n"),file=file_list[[i]])
        #print(summary(fit))
      }
  }

####################
# Close files
####################
for (i in 1:length(phenos)) {
  close(file_list[[i]])
}
