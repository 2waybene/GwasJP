args<-commandArgs(T)
path<-args[1]
#pheno <- args[2]
#rm(list=ls())
#setwd('~/data/accord/data/analysis/')
#path <- '~/data/accord/data/analysis/hgi_baseline/hgi_black'
#path <- '~/data/accord/data/analysis/hgi_baseline/hgi_hispanic'
library(data.table)

pheno<-unlist(read.table(paste0(path,"/phenotypes.txt"),stringsAsFactors=F))
chr<-1:22
out<-file(paste0(path,"/outputs/gc/gc_info.txt"),open="w")

#########################
# Load data
#########################
# Get results from chip
# data<-read.table(pipe(paste0("awk '{if ($5==\"ADD\") print $1\"\t\"$2\"\t\"$9}' ",
# 	path,"/association_cv/chr0.",pheno[1],".assoc.linear")),header=F,stringsAsFactors=F)
# names(data)<-c("CHR","SNP",pheno[1])
# if (length(pheno)>1) {
# 	for (pheno_i in 2:length(pheno)) {
# 		data<-cbind(data,read.table(pipe(paste0("awk '{if ($5==\"ADD\") print $9}' ",
# 			path,"/association_cv/chr0.",pheno[pheno_i],".assoc.linear")),header=F,stringsAsFactors=F))
# 		names(data)[pheno_i+2]<-pheno[pheno_i]
# 	}
# }
#frq<-read.table(pipe(paste0("awk '{print $1,$2,$5}' ",path,"/association_cv/plink.frq")),header=T,stringsAsFactors=F)
#frq<-read.table(file.path(path,"association_cv","plink.frq"),header=T,stringsAsFactors=F)
# data<-merge(data,frq)
data <- fread(file.path(path,"association_cv",paste0("all_commonVariantAnalysis_",pheno[1],"_results.tsv")))
data <- data[,c("SNP","CHR","P_unadj","Type","MAF"),with=FALSE]
colnames(data)[which(colnames(data)=="P_unadj")] <- pheno[1] 
gc_probes<-unlist(read.table("/home/accord/data/geno_data/LD.prune.in",header=F,stringsAsFactors=F))
#gc_probes<-unlist(read.table("~/NCSU_Projects/Accord/Data/LD.prune.in",header=F,stringsAsFactors=F))
gc_idx<-data$SNP%in%gc_probes
data <- data[gc_idx,]
data<-data[data$Type%in%c("GENO","META") & data$MAF>0,]
#data<-data[complete.cases(data),]

if (length(pheno)>1) {
	for (pheno_i in 2:length(pheno)){
	  new.data <- fread(file.path(path,"association_cv",paste0("all_commonVariantAnalysis_",pheno[pheno_i],"_results.tsv")))
	  new.data <- new.data[,c("SNP","CHR","P_unadj","Type","MAF"),with=FALSE]
	  colnames(new.data)[which(colnames(new.data)=="P_unadj")] <- pheno[pheno_i] 
	  gc_idx<-new.data$SNP%in%gc_probes
	  new.data <- new.data[gc_idx,]
	  new.data<-new.data[new.data$Type%in%c("GENO","META") & new.data$MAF>0,]
	  #data<-data[complete.cases(data),]
	  data <-merge(data,new.data,by=c("SNP","CHR","Type","MAF"))
	}
}

gc()



#########################
# Plot genomic control vs MAF threshold
#########################
cutoffs<-seq(0,0.05,0.01)
gc<-list()
n<-c()
for (pheno_i in 1:length(pheno)){
	temp<-c()
	set(data,NULL,pheno[pheno_i],as.numeric(data[[pheno[pheno_i]]]))
	for (i in 1:length(cutoffs)) {
		chi2_stats<-qchisq(unlist(data[data$MAF>=cutoffs[i],pheno[pheno_i],with=FALSE]),df=1,lower.tail=F)
		n[i]<-length(chi2_stats)
		temp[i]<-median(chi2_stats)/qchisq(0.5,df=1)
	}
	gc[[pheno_i]]<-temp
}
cat("MAF threshold\n",
    "  MAF cutoffs: ",cutoffs,"\n",
    "  # probes:    ",n,"\n",file=out)
for (pheno_i in 1:length(pheno)) {
	cat ("  ",pheno[pheno_i],": ",gc[[pheno_i]],"\n",file=out)
}

pdf(paste0(path,"/outputs/gc/gc_all_maf_threshold.pdf"),width=9,height=4.5)
par(mar=c(3.5,3.5,2,1))
plot(cutoffs,gc[[1]],col=1,pch=1,type="b",
     ylim=c(min(1,unlist(gc)),max(1,unlist(gc))),xlab="",ylab="",main="")
if (length(pheno)>1) {
	for (p in 2:length(pheno)) {
		points(cutoffs,gc[[p]],col=p,pch=p,type="b")
	}
}
abline(a=1,b=0,lty=2)
mtext("MAF threshold",side=1,line=2)
mtext(bquote(widehat(lambda)),side=2,line=2)
legend("top",legend=pheno,col=c(1:length(pheno)),pch=c(1:length(pheno)),
       xpd=TRUE,horiz=TRUE,inset=-.1,bty="n")
graphics.off()

##Start

#########################
# Plot genomic control vs MAF quantile
#########################
n_bin<-7
cutoffs<-quantile(data$MAF,seq(0,1,length=n_bin+1))
#?quantile
#hist(data$MAF)
#graphics.off()
bin<-list()
bin[[1]]<-data[data$MAF<=cutoffs[2],]
#head(bin[[7]])
summary(data$MAF)
#plot(data$MAF)
for (i in 2:(n_bin-1)) {bin[[i]]<-data[data$MAF>cutoffs[i] & data$MAF<=cutoffs[i+1],]}
bin[[n_bin]]<-data[data$MAF>=cutoffs[n_bin],]
for (i in 1:n_bin) {
	bin[[i]]<-bin[[i]][complete.cases(bin[[i]]),]
	bin[[i]]<-cbind(bin[[i]],ppoints(dim(bin[[i]])[1]))
	names(bin[[i]])[dim(data)[2]+1]<-"x"
}
cat(paste0("MAF quantile bins:\n",
    "  MAF cutoffs: ",cutoffs,"\n",
    "  # probes:    ",sapply(lapply(bin,dim),"[",1),"\n",file=out))

for (pheno_i in 1:length(pheno)) {
  pdf(paste0(path,"/outputs/gc/gc_",pheno[pheno_i],"_quantile_qq.pdf"),width=7.5, height=7.5)
  layout(matrix(c(1:min(n_bin,9),rep(0,9-min(n_bin,9))),3,3,byrow=T))
  for (i in 1:min(n_bin,9)) {
    par(mar=c(0,0,1,0))
    plot(head(-log10(bin[[i]]$x),n=1500),head(-log10(sort(unlist(bin[[i]][,pheno[pheno_i],with=FALSE]))),n=1500),
         axes=F,xlab="",ylab="",main=paste0(cutoffs[i]," < MAF <= ",cutoffs[i+1]))
    abline(a=0,b=1,lty=2)
  }
  layout(1)
  graphics.off()
}

gc<-list()
for (i in 1:length(pheno)) {
  gc[[i]]<-unlist(lapply(bin,function(x) qchisq(median(unlist(x[,pheno[i],with=FALSE])),1,lower.tail=F)/qchisq(0.5,1)))
}

pdf(paste0(path,"/outputs/gc/gc_all_quantile.pdf"),width=9, height=4.5)
par(mar=c(3.5,3.5,2,1))
x<-log10(min(data$MAF[data$MAF>0])*cutoffs[2])/2
x<-c(x,log10(cutoffs[2:n_bin]*cutoffs[3:(n_bin+1)])/2)
plot(x,gc[[1]],ylim=c(min(1,unlist(gc)),max(1,unlist(gc))),
    type="b",xlab="",ylab="")
if (length(pheno)>1) {
	for (i in 2:length(pheno)) {
	  points(x,gc[[i]],type="b",col=i)
	}
}
abline(h=1,lty=2)
abline(v=log10(seq(0.01,0.1,0.01)),lty=2,col=2)
abline(v=log10(0.03),lwd=2,col=2)
abline(v=log10(cutoffs),lty=1)
mtext("MAF quantile bin",side=1,line=2)
mtext(bquote(widehat(lambda)),side=2,line=2)
legend("top",legend=pheno,col=c(1:length(pheno)),pch=c(1:length(pheno)),
       xpd=TRUE,horiz=TRUE,inset=-.1,bty="n")
graphics.off()

#########################
# Plot genomic control vs MAF bin
#########################
n_bin<-10
cutoffs<-seq(0,n_bin*0.01,0.01)
bin<-list()
bin[[1]]<-data[data$MAF<cutoffs[2],]
for (i in 2:(n_bin-1)) {bin[[i]]<-data[data$MAF>=cutoffs[i] & data$MAF<cutoffs[i+1],]}
bin[[n_bin]]<-data[data$MAF>=cutoffs[n_bin],]
for (i in 1:n_bin) {
  bin[[i]]<-bin[[i]][complete.cases(bin[[i]]),]
  bin[[i]]<-cbind(bin[[i]],ppoints(dim(bin[[i]])[1]))
  names(bin[[i]])[dim(data)[2]+1]<-"x"
}
cat(paste0("MAF bins:\n",
    "  MAF cutoffs: ",cutoffs,"\n",
    "  # probes:    ",sapply(lapply(bin,dim),"[",1),"\n",file=out))

pheno_i <- 1
for (pheno_i in 1:length(pheno)) {  
  pdf(paste0(path,"/outputs/gc/gc_",pheno[pheno_i],"_bin_qq.pdf"),width=7.5, height=7.5)
  layout(matrix(c(1:min(n_bin,9),rep(0,9-min(n_bin,9))),3,3,byrow=T))
  i <- 1
  for (i in 1:min(n_bin,9)) {
    par(mar=c(0,0,1,0))
    plot(head(-log10(bin[[i]]$x),n=1500),head(-log10(sort(unlist(bin[[i]][,pheno[pheno_i],with=FALSE]))),n=1500),
         axes=F,xlab="",ylab="",main=paste0(cutoffs[i]," <= MAF <= ",cutoffs[i+1]))
    abline(a=0,b=1,lty=2)
  }
  layout(1)
  graphics.off()
}

gc<-list()
for (i in 1:length(pheno)) {
  gc[[i]]<-unlist(lapply(bin,function(x) qchisq(median(unlist(x[,pheno[i],with=FALSE])),1,lower.tail=F)/qchisq(0.5,1)))
}

pdf(paste0(path,"/outputs/gc/gc_all_bin.pdf"),width=9, height=4.5)
par(mar=c(3.5,3.5,2,1))
plot(cutoffs[1:n_bin]+0.005,gc[[1]],ylim=c(min(1,unlist(gc)),max(1,unlist(gc))),
    type="b",xlab="",ylab="")
if (length(pheno)>1) {
	for (i in 2:length(pheno)) {
	  points(cutoffs[1:n_bin]+0.005,gc[[i]],type="b",col=i)
	}
}
abline(h=1,lty=2)
abline(v=cutoffs,lty=2)
abline(v=0.03,lwd=2,col=2)
mtext("MAF bin",side=1,line=2)
mtext(bquote(widehat(lambda)),side=2,line=2)
legend("top",legend=pheno,col=c(1:length(pheno)),pch=c(1:length(pheno)),
       xpd=TRUE,horiz=TRUE,inset=-.1,bty="n")
graphics.off()

close(out)
