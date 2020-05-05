args<-commandArgs(T)
path<-args[1]
pheno<-args[2]
# path<-"cholesterol/log_blr_all"
#setwd("Y://NCSU_Projects//Accord//Data//analysis")
#setwd("~/NCSU_Projects/Accord/Data/analysis/")
# path <- "Metformin/met_90d_allraces4"
# pheno<-"met_hba1c"

q_breaks<-c(0.2,0.15,0.1)
q_cols<-c("green","blue","red")

library(qvalue)
library(hexbin)
#library(qvalue,lib="~/R/library")
#library(hexbin,lib="~/R/library")

############################################################
# Combine p-value function
############################################################
# Correlated Lancaster procedure (see Dai et al. paper)
Lancaster<-function(p,w) {
  n<-dim(p)[2]
  # Rescale w for missing p-values, make all wi sum to 1 regardless of missing values
  w<-w/sum(w)
  ww<-t(apply(p,1,function(x) {temp<-w/sum(w[!is.na(x)]);temp[is.na(x)]<-NA;return(temp)}))
  # Compute Xsq statistics
  Xsq<-t(apply(cbind(p,ww),1,function(x) {qgamma(1-x[1:n],shape=x[(n+1):(2*n)]/2,scale=2)}))
  # Compute Lancaster test statistics
  t<-apply(Xsq,1,function(x) {sum(x,na.rm=T)})
  # Compute E[T] and Var[T] using complete observations
  rho<-cov(Xsq,use="complete.obs")
  et<-sum(w)
  vart<-2*sum(w)+(sum(rho)-sum(diag(rho)))
  # Compute scaled Xsq variables
  v<-2*et^2/vart
  c<-v/et
  return(pchisq(c*t,v,lower.tail=F))
}

############################################################
# Read data function
############################################################
GetData<-function(path,pheno) {
  h<-read.table(paste0(path,"/association_rv/chr1_",pheno,".txt"),nrows=1,stringsAsFactors=F)
  data<-read.table(pipe(paste0("tail -qn +2 ",path,"/association_rv/chr*_",pheno,".txt")),header=F,stringsAsFactors=F)
  names(data)<-h
  data<-data[,grepl("geneID|_p$",h)]
  return(data)
}

############################################################
# Run analysis
############################################################
# Load data
data<-GetData(path,pheno)
# Combine p-values
p<-as.matrix(data[,-1])
w<-rep(1,dim(p)[2])/dim(p)[2]
data<-cbind(data,Lancaster(p,w))
names(data)[dim(data)[2]]<-c("Lancaster")

p_val_source<-"Lancaster"
q_data<-data[!is.na(data[,p_val_source]),]
#print("Hello")
#qobj<-qvalue(as.numeric(q_data[,p_val_source])) ## JJ COMMENTED FOR p.adjust
qobj<-p.adjust(as.numeric(q_data[,p_val_source]))
#print("Hello2")
#q_data<-cbind(q_data,qobj$qvalues) ## JJ COMMENTED FOR p.adjust
q_data<-cbind(q_data,qobj)
names(q_data)[dim(q_data)[2]]<-"qvalue"
q_data<-q_data[order(q_data[p_val_source]),]

# qplot(qobj)
# title(paste0(pheno,"  ",p_val_source),outer=T,line=-1)
# hist(q_data[,p_val_source],freq=F,main=paste0(pheno,"  ",p_val_source))
# abline(h=1,lty=2)

pdf(paste0(path,"/outputs/",pheno,"_fdr.pdf"),width=10,height=7.5)
x<-ppoints(q_data[,p_val_source])
col<-rep("gray",length(x))
for (i in 1:length(q_breaks)) {
  col[q_data$qvalue<=q_breaks[i]]=q_cols[i]
}
plot(-log10(x),-log10(q_data[,p_val_source]),
     col=col,
     xlab="-log10(expected p-value)",
     ylab="-log10(observed p-value)",
     main=paste0(pheno,"  ",p_val_source))
abline(a=0,b=1,lty=2)
legend("topleft",legend=paste0("q<=",q_breaks),pch=1,col=q_cols)
trash<-dev.off()

#write.table(format(q_data[q_data$qvalue<=0.3,],digits=4),paste0(path,"/outputs/",pheno,"_fdr.txt"),row.names=F,quote=F,sep="\t")
write.table(format(q_data,digits=4),paste0(path,"/outputs/",pheno,"_fdr.txt"),row.names=F,quote=F,sep="\t")

############################################################
# Look at correlation between p-values from various tests
############################################################
z<-cor(data[,-1],use="complete.obs")
n<-dim(z)[1]
colPalette<-colorRampPalette(c("white","black"))
m<-c();
count<-1;for (i in 1:n) {for (j in 1:n) {if (i<=j) {m<-c(m,count);count<-count+1} else {m<-c(m,0)}}}
png(paste0(path,"/outputs/",pheno,"_fdr_cor.png"),width=7.5,height=7.5,units="in",res=150)
par(oma=c(0,2,2,0))
layout(matrix(m,n,n,byrow=T))
for (i in 1:n) {
  for (j in i:n) {
    if (i==j) {
      p<-data[!is.na(data[,i+1]),i+1]
      #qobj<-qvalue(p) ## JJ comment for p.adjust
      qobj <- p.adjust(p)
      col<-rep("gray",length(p))
      for (k in 1:length(q_breaks)) {
          #col[qobj$qvalues<=q_breaks[k]]=q_cols[k] ## JJ comment for p.adjust
          col[qobj<=q_breaks[k]]=q_cols[k] 
      }
      #q_data<-data.frame(p,qobj$qvalues,col,stringsAsFactors=F) ## JJ comment for p.adjust
      q_data<-data.frame(p,qobj,col,stringsAsFactors=F)
      names(q_data)<-c("p","q","col")
      q_data<-q_data[order(p),]
      y<- -log10(q_data$p)
      x<- -log10(ppoints(y))
      par(mar=c(0,0,0,0))
      plot(x,y,xlim=c(0,1.05*max(x)),ylim=c(0,1.05*max(y)),col=q_data$col,
           pch=16,cex=0.5,xlab="",ylab="",main="",xaxs="i",yaxs="i",axes=F);box()
      abline(a=0,b=1,lty=2)
#       text(0,0.98*max(y),as.character(length(x)),pos=4)
      par(xpd=NA);text(-max(x)*.1,max(y)/2,names(data)[i+1],srt=90);par(xpd=F)
    } else {
      y<- data[,i+1]
      x<- data[,j+1]
      fit<-lm(y~x)
      par(mar=c(0,0,0,0))
      bin<-hexbin(x,y,xbins=100)
      plot(bin@xcm,bin@ycm,pch=16,cex=0.5,
           col=colPalette(max(15))[sapply(as.numeric(bin@count),function(x) {min(x,15)})],
           xlab="",ylab="",main="",axes=F);box()
      abline(a=fit$coef[1],b=fit$coef[2],lty=2,col=2,lwd=2)
#       text(0,0.98*max(y),as.character(length(x)),pos=4)
      text(max(x,na.rm=T),0.05*max(y,na.rm=T),format(z[i,j],digits=4),pos=2)
    }
  }
}
title(paste0(path,"     ",pheno),outer=T,line=0.5,cex=2)
trash<-dev.off()
