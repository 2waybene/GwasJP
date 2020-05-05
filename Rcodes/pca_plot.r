args<-(commandArgs(TRUE))
path<-args[[1]]

print('##########')
print('pca_plot.r')
print('##########')

info<-read.table(paste0(path,"/pca/result.pca.evec"),header=F,stringsAsFactors=F)
race<-as.factor(info[,dim(info)[2]])
n<-length(unique(race))

t<-rep(0,length(race))
t[race=="Asian"]    <-1
t[race=="Black"]    <-2
t[race=="CEU"]      <-3
t[race=="Hispanic"] <-4
t[race=="Non-White"]<-5
t[race=="Other"]    <-6
t[race=="White"]    <-7
t[race=="YRI"]      <-8

t2<-rep(0,length(levels(race)))
t2[levels(race)=="Asian"]    <-1
t2[levels(race)=="Black"]    <-2
t2[levels(race)=="CEU"]      <-3
t2[levels(race)=="Hispanic"] <-4
t2[levels(race)=="Non-White"]<-5
t2[levels(race)=="Other"]    <-6
t2[levels(race)=="White"]    <-7
t2[levels(race)=="YRI"]      <-8


pdf(paste0(path,"/outputs/PC1-2.pdf"),width=7.5, height=7.5)
par(mar = c(0,0,1.5,0));plot(info[,2],info[,3],col=t,pch=t,axes=F,xlab="",ylab="",main="PCA");
legend("bottomright",levels(race),col=t2,pch=t2)
trash<-dev.off()

pdf(paste0(path,"/outputs/PC1-5.pdf"),width=7.5, height=7.5)
M<-matrix(c(1,2,3,4,5,0,6,7,8,9,0,0,10,11,12,0,0,0,13,14,0,0,0,0,15),byrow=T,nrow=5)
layout(M)
par(mar = c(0,0,0,0));plot(info[,2],info[,2],col=t,pch=t,axes=F,xlab="",ylab="");box()
par(mar = c(0,0,0,0));plot(info[,2],info[,3],col=t,pch=t,axes=F,xlab="",ylab="");box()
par(mar = c(0,0,0,0));plot(info[,2],info[,4],col=t,pch=t,axes=F,xlab="",ylab="");box()
par(mar = c(0,0,0,0));plot(info[,2],info[,5],col=t,pch=t,axes=F,xlab="",ylab="");box()
par(mar = c(0,0,0,0));plot(info[,2],info[,6],col=t,pch=t,axes=F,xlab="",ylab="");box()
par(mar = c(0,0,0,0));plot(info[,3],info[,3],col=t,pch=t,axes=F,xlab="",ylab="");box()
par(mar = c(0,0,0,0));plot(info[,3],info[,4],col=t,pch=t,axes=F,xlab="",ylab="");box()
par(mar = c(0,0,0,0));plot(info[,3],info[,5],col=t,pch=t,axes=F,xlab="",ylab="");box()
par(mar = c(0,0,0,0));plot(info[,3],info[,6],col=t,pch=t,axes=F,xlab="",ylab="");box()
par(mar = c(0,0,0,0));plot(info[,4],info[,4],col=t,pch=t,axes=F,xlab="",ylab="");box()
par(mar = c(0,0,0,0));plot(info[,4],info[,5],col=t,pch=t,axes=F,xlab="",ylab="");box()
par(mar = c(0,0,0,0));plot(info[,4],info[,6],col=t,pch=t,axes=F,xlab="",ylab="");box()
par(mar = c(0,0,0,0));plot(info[,5],info[,5],col=t,pch=t,axes=F,xlab="",ylab="");box()
par(mar = c(0,0,0,0));plot(info[,5],info[,6],col=t,pch=t,axes=F,xlab="",ylab="");box()
par(mar = c(0,0,0,0));plot(info[,6],info[,6],col=t,pch=t,axes=F,xlab="",ylab="");box()
trash<-dev.off()
