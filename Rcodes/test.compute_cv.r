rm(list=ls())

require(data.table)
debug <- FALSE
debug <- TRUE
if (debug){
    path <- '~/data/accord/data/analysis/icaps/icaps.thiazide.all_races/'
    setwd('~/data/accord/data/analysis/')
    chr<-8
    chunk<-1
}else{
    args <- commandArgs(T)
    path <- args[1]
    chr <- args[2]
    chunk <- args[3]   
}


# Check if imputed genotype file exists
f_gz<-paste0("../imputation/outputs/chr",chr,".imputed.",chunk,".gz")
if (!file.exists(f_gz)) {
  cat("File",f_gz,"does not exist\n")
  quit()
}

if (debug){
    f_gz2<-paste0("../imputation/outputs/test/chr",chr,".imputed.",chunk,"")
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
  file_list[[i]]<-file(paste0(path,"/association_cv/imputed_chunks/chr",chr,".",chunk,".",phenos[i],".assoc"),open="w")
  if(model.types[i]=="linear"){
      cat("CHR\tSNP\tBP\tA1\tA2\tFRQ\tHWE_CHI2\tBETA\tSE\tP\n",file=file_list[[i]])
  }else{
     cat("CHR\tSNP\tBP\tA1\tA2\tFRQ\tHWE_CHI2\tOR\tSE\tP\n",file=file_list[[i]])
 }

}

if(!is.null(snp.index)){
  read.range <- snp.index-1
  read_size <- 1
}else{
  read.range <- 0:ceiling(1e9/read_size)
}

#read_i <- read.range[-1]
#read_i <- read.range[length(read.range)]

total.fread.time <- 0
total.read.table.time <- 0
for (read_i in read.range) {
    ####################
    ## Load genotype data
    ####################
    ##obj<-try(geno_data<-read.table(gzfile(f_gz),header=F,stringsAsFactors=F,nrow=read_size,skip=read_i*read_size),silent=T)
    ##if (is(obj, "try-error")) break
    cat('*****************\n ')
    cat('fread: ')
    old <- Sys.time()
    obj <- try(geno_data<-fread(paste('zcat -c', f_gz),header=F,showProgress=T,nrow=read_size,skip=read_i*read_size),silent=T)
    ## (is(obj,"try-error")){
    #    print("error on fread...empty chunk or end of file. Most likely there is not a problem because all files have an end.")
    #                                    #continue#
    #}
    tmp.time <- Sys.time()-old
    total.fread.time <- total.fread.time + tmp.time
    cat(tmp.time,'\n')
    cat('Total fread time = ',total.fread.time,'\n')
    cat('readtable: ')
    old <- Sys.time()
    obj <- try(geno_data<-read.table(gzfile(f_gz),header=F,stringsAsFactors=F,nrow=read_size,skip=read_i*read_size),silent=T)
    if (is(obj, "try-error")){
        tmp.time <- Sys.time()-old
        total.read.table.time <- total.read.table.time + tmp.time
        print(Sys.time()-old)
        print("Finished.")
        cat('Total fread time = ',total.fread.time,'\n')
        cat('Total read.table time = ',total.read.table.time,'\n')
        break
    }
    tmp.time <- Sys.time()-old
    total.read.table.time <- total.read.table.time + tmp.time
    cat(tmp.time,'\n')
    cat('Total read.table time = ',total.read.table.time,'\n')
    cat('geno_data # of rows = ',nrow(geno_data),'\n')
}
