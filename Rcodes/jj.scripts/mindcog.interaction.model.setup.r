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
    chr<-19
    chunk<-10
     f_gz <- '~/data/accord/data/geno_data/mindcog.apoe.chunk/chr19.imputed.10'
}else{
  ## DO NOT ADD DEBUG INFO IN HERE
    args <- commandArgs(T)
     path <- args[1]
     phename = args[2]
     chr <- 19
     chunk <- 10
     f_gz <- '~/data/accord/data/geno_data/mindcog.apoe.chunk/chr19.imputed.10'
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
head(pheno_data) ## The top ones are empty
dim(pheno_data)
# read in model types
model.types <- unlist(read.table(paste0(path,"/modeltypes.txt"),stringsAsFactors=F))

# Get list of samples with no NA
retained_samples<-complete.cases(pheno_data)
#head(retained_samples)
which(retained_samples)

# Initialize variables for HWE computation
#n_obs<-sum(retained_sampless)
#O<-c(0,0,0)
#E<-c(0,0,0)

#cat('hello2')
#check for snps in list
snp.index <- NULL
if(file.exists(file.path("~/data/accord/data/geno_data/mindcog.apoe.chunk/apoe.txt"))){
  #print("You have a snp_list.txt file!!!")
  snp.list <- read.table(file.path("~/data/accord/data/geno_data/mindcog.apoe.chunk/apoe.txt"),stringsAsFactors=F)
  system.time(file.snps<-read.table(pipe(paste0("cut -d ' ' -f2 ../imputation/outputs/chr",chr,".imputed.",chunk,"_info")),header=T,stringsAsFactors=F))
  #print(head(file.snps[,1]))
  #
  snp.index <- which(file.snps[,1]%in%as.character(snp.list[,1]))
  if(length(snp.index)==0){
    stop(paste0("No SNP's in snp_list.txt were found in ../imputation/outputs/chr",chr,".imputed.",chunk,"_info"))
  } #else{
    #  print("Found it!!!")
  #}
    
}
print("APOE snp list")
print(head(snp.list[,1]))

# Load genotype info metric
geno_data_info<-read.table(pipe(paste0("cut -d ' ' -f5 ../imputation/outputs/chr",chr,".imputed.",chunk,"_info")),header=T)
dim(geno_data_info)
head(geno_data_info)
geno_data_info<-geno_data_info>info_threshold
dim(geno_data_info)
geno_data_i<-1
dim(geno_data_info)

# Get index values for the 3 genotypes
geno_2_idx<-(1:length(retained_samples))*3
geno_1_idx<-geno_2_idx-1
geno_0_idx<-geno_2_idx-2
#head(geno_0_idx)
#cat('hello3')

####################
# Initialize output files
####################
#file_list<-list()
#for (i in 1:length(phenos)) {
#  file_list[[i]]<-file(paste0(path,"/association_cv/imputed_chunks/chr",chr,".",chunk,".",phenos[i],".assoc"),open="w")
#  if(model.types[i]=="linear"){
#      cat("CHR\tSNP\tBP\tA1\tA2\tFRQ\tHWE_CHI2\tBETA\tSE\tP\tN\n",file=file_list[[i]])
#  }else{ #### THIS IS NOT AN OR
#     cat("CHR\tSNP\tBP\tA1\tA2\tFRQ\tHWE_CHI2\tOR\tSE\tP\tN\n",file=file_list[[i]])
# }
#
#}

print(snp.index)
if(!is.null(snp.index)){
  read.range <- snp.index-1
  read_size <- 1
}else{
  read.range <- 0:ceiling(1e9/read_size)
}

#read_i <- read.range[length(read.range)]
#read_i <- read.range[1]

#cat('hello4')

cat(f_gz)
#read.range <- 1
read_i = read.range[1]
####################
# Load genotype data
####################
geno_data<-try(fread(f_gz,header=F,verbose=F,showProgress=F,nrow=read_size,skip=read_i*read_size))
if (any(class(geno_data)=="try-error")){
  cat(paste0('try-error: ',f_gz,'\n'))
  break
}
geno_data = data.frame(geno_data)
tmp = geno_data
read_i = read.range[2]
geno_data[2,]<-try(fread(f_gz,header=F,verbose=F,showProgress=F,nrow=read_size,skip=read_i*read_size))
head(geno_data)[,1:10]

# Create full list for all samples with imputed genotypes
samples<-read.table(pipe(paste0("cut -d ' ' -f1,2 ../imputation/data/chr1.unphased.fam")),header=F,stringsAsFactors=F)
samples$order<-seq(len=nrow(samples))
names(samples)<-c("FID","IID","order")

head(samples)

#pheno_data<-merge(samples,data,sort=F,all.x=T)
#pheno_data<-pheno_data[sort.list(pheno_data$order), -3]


dim(samples)
dim(pheno_data)
head(pheno_data)
head(pheno_data[retained_samples,])
identical(pheno_data[,1],samples[,1])
head(pheno_data)
head(samples)
head(geno_data)[,1:10]

## offset the genocalls
geno_0_idx = geno_0_idx+5
geno_1_idx = geno_1_idx+5
geno_2_idx = geno_2_idx+5

#helloworld = c(.3,.7,0)
#helloworld = c(.6,.3,.2)
#max = 999
#  if (helloworld[1]>helloworld[2]){
#    if (helloworld[1] > helloworld[3]){
#      max = 1
#    }else{
#      max = 3
#    }
#  }else{
#    if (helloworld[2] > helloworld[3]){
#      max = 2
#    }else{
#      max = 3
#    }
#  }
#  print(max)

pheno_data$rs429358="Q"
pheno_data$rs7412="Q"
pheno_data$E2 = 99
pheno_data$E4 = 99

head(pheno_data)
head(samples)
head(geno_data)[,1:10]
geno_data$V4[1] = 'T'
geno_data$V5[2] = 'T'

i = 1
for (i in 1:nrow(pheno_data)){
  max = 999

  ## rs429358
  if (geno_data[1,geno_0_idx[i]]>geno_data[1,geno_1_idx[i]]){
    if (geno_data[1,geno_0_idx[i]] > geno_data[1,geno_2_idx[i]]){
      # Homo A1
      pheno_data$rs429358[i] = 'TT'
    }else{
      # Homo A2
      pheno_data$rs429358[i] = 'CC'
    }
  }else{
    if (geno_data[1,geno_1_idx[i]] > geno_data[1,geno_2_idx[i]]){
     # Hetero
    pheno_data$rs429358[i] = 'TC'
    }else{
      pheno_data$rs429358[i] = 'CC'
    }
  }

  ## rs7412
  if (geno_data[2,geno_0_idx[i]]>geno_data[2,geno_1_idx[i]]){
    if (geno_data[2,geno_0_idx[i]] > geno_data[2,geno_2_idx[i]]){
      # Homo A1
      pheno_data$rs7412[i] = 'CC'
    }else{
      # Homo A2
      pheno_data$rs7412[i] = 'TT'
    }
  }else{
    if (geno_data[2,geno_1_idx[i]] > geno_data[2,geno_2_idx[i]]){
     # Hetero
    pheno_data$rs7412[i] = 'CT'
    }else{
      pheno_data$rs7412[i] = 'TT'
    }
  }
  if( (length(grep("T",pheno_data$rs7412[i]))>0) & (length(grep("T",pheno_data$rs429358[i]))>0) ){
      pheno_data$E2[i] = 1
  }else{
      pheno_data$E2[i] = 0
  }
  if( (length(grep("C",pheno_data$rs7412[i]))>0) & (length(grep("C",pheno_data$rs429358[i]))>0) ){
      pheno_data$E4[i] = 1
  }else{
      pheno_data$E4[i] = 0
  }
}

## Drop columns you don't want to model
pheno_data$rs429358 = NULL
pheno_data$rs7412 = NULL

head(pheno_data)
pheno_data[8012,]

#data<-read.table(paste0(path,"/pheno_data/pheno_data_step2.txt"),header=T,stringsAsFactors=F)
#head(data)
#head(pheno_data)
dir.create(file.path(path, "apoe"), showWarnings = FALSE)
dir.create(file.path(path, "apoe","association_cv"), showWarnings = FALSE)
dir.create(file.path(path, "apoe","association_cv","imputed_chunks"), showWarnings = FALSE)
## Save new phenotype file
models<-unlist(read.table(paste0(path,"/pheno_data/models.txt"),stringsAsFactors=F))
models = paste0(models,"+E2+E4")
## Save new model file
write.table(models,paste0(path,"/apoe/models.txt"),col.names=F,row.names=F,quote=F)
write.table(pheno_data,paste0(path,"/apoe/pheno_data.apoe.txt"),col.names=T,row.names=F,quote=F,sep='\t')

## Make new pheno files
cov = read.table(file.path(path,"pheno_data",paste0("covar_",phename,".txt")),header=T)
#head(cov)
#table(cov$E2)
cov$E2 = 99
cov$E4 = 99
#head(pheno_data)
for (i in 1:nrow(cov)){
  cov$E2[i] = pheno_data$E2[which(pheno_data$FID %in% cov$FID[i])]
  cov$E4[i] = pheno_data$E4[which(pheno_data$FID %in% cov$FID[i])]
}

print(file.path(path,"apoe",paste0("covar_",phename,".txt")))
write.table(cov,file.path(path,"apoe",paste0("covar_",phename,".txt")),col.names=T,row.names=F,quote=F,sep='\t')


### ACR SETUP

dir.create(file.path(path, "acr"), showWarnings = FALSE)
dir.create(file.path(path, "acr","association_cv"), showWarnings = FALSE)
dir.create(file.path(path, "acr","association_cv","imputed_chunks"), showWarnings = FALSE)

## Make ACR file
models<-unlist(read.table(paste0(path,"/pheno_data/models.txt"),stringsAsFactors=F))
models = paste0(models,"+log10ACR")



#head(pheno_data)
pheno_data$E2 = NULL
pheno_data$E4 = NULL

d = read.table('~/data/accord/data/pheno_data/mindcog/accord_analysis_file04_16_2018.txt',header=T,fill=NA,sep='\t')
#head(d)

pheno_data$log10ACR = 99
i = 8012
tmp.i = which(pheno_data$FID %in% d$MaskID)
i = tmp.i[1]
#head(pheno_data[tmp.i,])
for (i in 1:nrow(pheno_data)){
  if (pheno_data$FID[i] %in% d$MaskID){
    pheno_data$log10ACR[i] = log10(d$acr[which(d$MaskID %in% as.integer(pheno_data$FID[i]))]+1)
  }
}

#typeof(d$MaskID)
#typeof(as.numeric(pheno_data$FID[i]))

write.table(models,paste0(path,"/acr/models.txt"),col.names=F,row.names=F,quote=F)
write.table(pheno_data,paste0(path,"/acr/pheno_data.acr.txt"),col.names=T,row.names=F,quote=F,sep='\t')

## Make new pheno files
cov = read.table(file.path(path,"pheno_data",paste0("covar_",phename,".txt")),header=T)
#head(cov)
cov$log10ACR = 99
#head(pheno_data)
for (i in 1:nrow(cov)){
  cov$log10ACR[i] = pheno_data$log10ACR[which(pheno_data$FID %in% cov$FID[i])]
}

print(file.path(path,"acr",paste0("covar_",phename,".txt")))
write.table(cov,file.path(path,"acr",paste0("covar_",phename,".txt")),col.names=T,row.names=F,quote=F,sep='\t')


