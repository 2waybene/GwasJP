rm(list=ls())
debug <- FALSE
#debug <- TRUE
if (debug){
    path <- '~/data/accord/data/analysis/mindcog/mindcog_cog_dsst/'
    setwd('~/data/accord/data/analysis/')
    adjusted <- FALSE
    
}else{
    args<-(commandArgs(TRUE))
    path<-args[[1]]
}

#path <- 'hgi_baseline/hgi_white/'
# Load initial phenotype data
data<-read.table(paste0(path,"/pheno_data/pheno_data_step1.txt"),header=T,stringsAsFactors=F)

# Remove selected variables
if (file.exists(paste0(path,"/remove_covars.txt"))) {
	discard<-unlist(read.table(paste0(path,"/remove_covars.txt"),stringsAsFactors=F))
	data<-data[,!names(data)%in%discard]
}


# Remove relatedness discard list
obj<-try(discard<-read.table(paste0(path,"/relatedness/discard.txt"),header=F,stringsAsFactors=F),silent=T)
if (!is(obj, "try-error")) {
	data<-data[!(data$MaskID%in%discard[,1] & data$LabID%in%discard[,2]),]
}

# Add first 10 principle components
pc_data<-read.table(paste0(path,"/pca/result.pca.evec"),header=F,stringsAsFactors=F)
pc_data<-cbind(data.frame(do.call('rbind', strsplit(pc_data[,1],'_'))),pc_data[,-c(1,dim(pc_data)[2])])
colnames(pc_data)<-c("MaskID","LabID","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10")
data<-merge(data,pc_data,sort=F)

# Print output
write.table(data,paste0(path,"/pheno_data/pheno_data_step2.txt"),col.names=T,row.names=F,quote=F,sep="\t")
