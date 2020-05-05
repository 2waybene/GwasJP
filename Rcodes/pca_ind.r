rm(list=ls())

print('#########')
print('pca_ind.r')
print('#########')

debug <- FALSE
#debug <- TRUE
if (debug){
    path <- '~/data/accord/data/analysis/icaps/icaps.aceia.arb.all_races/'
    #path <- '~/data/accord/data/analysis/mindcog/mindcog_cog_dsst'
    #fname = 'pheno_data_mindcog.txt'
    fname = 'pheno_data_icaps.txt'
    #ath <- '~/data/accord/data/analysis/Metformin/met_90d_allraces4/'
    adjusted <- FALSE
    getwd()
    setwd('~/data/accord/data/analysis/')
}else{
    args<-(commandArgs(TRUE))
   print("pca_ind.r args:")
    print(args)
    path<-args[[1]]
    fname <- args[[2]]
}


data<-read.table(paste0(path,"/pca/data_pruned.tfam"),header=F,stringsAsFactors=F)
#head(data)
data<-data[,-c(3,4,6)]

names(data)<-c("MaskID","LabID","gender")
data[data$gender==1,3]<-"M"
data[data$gender==2,3]<-"F"

samples<-read.table(fname,header=T,stringsAsFactors=F,sep='\t')
samples<-samples[,names(samples)%in%c("MaskID","LabID","ethnicity")]
#head(samples)

m<-merge(data,samples,sort=F)
#head(m)
out<-cbind(apply(m,1,function(x) paste(x[1],"_",sub("^\\s+","",x[2]),sep="")),m[,3:4])

print("Writing ind.txt")
write.table(out,paste0(path,"/pca/ind.txt"),col.names=F,row.names=F,quote=F,sep="\t")
