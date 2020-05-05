args<-(commandArgs(TRUE))
path<-args[[1]]

info<-read.table(paste0(path,"/relatedness/king.kin0"),header=T,stringsAsFactors=F)
info<-info[info$Kinship>0.1768,1:4]

set.seed(1485)
discard<-runif(dim(info)[1])>0.5

FID<-info[,1]
IID<-info[,2]
FID[discard]=info[discard,3]
IID[discard]=info[discard,4]
data<-data.frame(FID,IID)
write.table(data,paste0(path,"/relatedness/discard.txt"),col.names=F,row.names=F,quote=F,sep="\t")
