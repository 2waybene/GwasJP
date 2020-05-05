rm(list=ls())
library(data.table) ##install.packages('data.table',type="source",repos = "http://Rdatatable.github.io/data.table")
library(scales)#install.packages('scales')
library(GenABEL)#install.packages('GenABEL')
debug <- FALSE
debug <- TRUE
if (debug){
    path <- '~/data/accord/data/analysis/Sulfonylurea/sulf_90d_allraces2'
    pheno <- 'sulf_hba1c'
    setwd('~/data/accord/data/analysis')
    maf_cutoff <- .03
    chr <- seq(1,22) 
}else{
    args<-commandArgs(T)
    path<-args[1] ##path <- '~/data/accord/data/analysis/Metformin/met_90d_allraces4'
    pheno<-args[2] ##pheno <- 'met_hba1c'
    maf_cutoff<-as.numeric(args[3]) ##maf_cutoff <- .03
    #chr<-sort(unique(as.numeric(args[-c(1,2,3)]))) ##chr <- seq(1,22) 
    chr<- seq(1,22) 
}



#########################
# Color palette
#########################
# Color for c(Genotypes,Imputed (All Markers/Subjects), Meta Analysis)
qq_colors<-c('black','steelblue1','violetred1')

# Colors to alternate in manhattan plot for imputed variants
manhattan_colors<-c('grey0','grey21','steelblue3','steelblue1','violetred3','violetred1')
## Decided not to change alpha blending
manhattan_alpha<-c(1,1)
#manhattan_background_color<-c('grey89','gray100')

##################################
## LOAD ALL ASSOCIATION RESULTS ##
##################################
data <- fread(paste0(path,"/association_cv/all_commonVariantAnalysis_",pheno,"_results.tsv"))

gc()

head(data)

#if (debug){
#    data[data$SNP=='rs112568157',]
#    data[data$SNP=='rs34094931',]
#    min(data$P_unadj)
#    length(which(data$P_unadj==0))
#    data[which(data$P_unadj==0),]
#    unique(data$Type)
#    dim(data[data$CHR<=22,])
#}

if (length(which(data$CHR>22))>0){
    data <- data[data$CHR<=22,]
    dim(data)
}
if (length(which(data$P_unadj==0))>0){
    data.zeros<- data[which(data$P_unadj==0),]
    dim(data.zeros)
    data <- data[-which(data$P_unadj==0),]
    dim(data)
}

## Filter by MAF (should be decided by looking at GC)
dim(data)
data <- data[MAF>.03,]
print(data[which(data$MAF>.03)&which(data$Type=='GENO'),])
print(dim(data[which(data$MAF>.03)&which(data$Type=='GENO'),]))
dim(data)
gc()

#############################
## Compute genomic control ##
#############################
#gc_probes<-unlist(read.table("/home2/accord/geno_data/LD.prune.in",header=F,stringsAsFactors=F))
#gc_probes<-unlist(read.table("/home/accord/data/geno_data/LD.prune.in",header=F,stringsAsFactors=F)) 
#gc_idx<-which(data[,SNP]%in%gc_probes)
#length(gc_idx)
#head(gc_idx)
#table(data$Type[gc_idx]) ## We do not want any 
#tmp <- intersect(gc_idx,which(data$Type=='IMPU'))
#length(tmp)
#if (length(tmp)>0){
#    gc_idx <- gc_idx[-which(gc_idx %in% intersect(gc_idx,which(data$Type=='IMPU')))]
#}
#length(gc_idx)
#gc_idx<-which(data[,SNP]%in%gc_probes)
#gc_p<-data[gc_idx,P_unadj]
#chi2_stats<-qchisq(na.omit(gc_p),df=1,lower.tail=F)
#gc <- median(chi2_stats)/qchisq(0.5,df=1)
#gc()

save.image('bin/jj.scripts/gc.correction.RData')

getwd()
load('gc.correction.RData')

head(data)
head(gc_probes)
## LD Pruned
gc_probes<-unlist(read.table("/home/accord/data/geno_data/LD.prune.in",header=F,stringsAsFactors=F)) 
gc_idx<-which(data[,SNP]%in%gc_probes)
gc_p <- data[gc_idx,P_unadj]
gc.ld <- estlambda(gc_p)$estimate
## 

head(gc_idx)
head(data)
gc <- estlambda(data[,"P_unadj"])$estimate
print(gc)
#gc <- gc$estimate
## Dataput info
output<-file(paste0(path,"/outputs/",pheno,"_info.txt"),open="w")
cat("Phenotype:         ",pheno,"\n",
    #"MAF cutoff:        ",maf_cutoff,"\n",
     #"Chromosomes:       ",paste(chr,collapse=" "),"\n",
     "Genomic inflation: ",format(gc,digits=4),"\n",
     "Genotyped:         ",nrow(data[Type=='GENO',]),"\n",
     "Imputed (all):           ",nrow(data[Type=='IMPU',]),"\n",
     "Meta:           ",nrow(data[Type=='META',]),"\n",
     sep="",file=output)
close(output)
gc()

#########################
# Adjust p-values
#########################
chi2_stats<-qchisq(data[,"P_unadj"],df=1,lower.tail=F)
head(chi2_stats)
chi2_stats<-chi2_stats/gc
adjP<-pchisq(chi2_stats,df=1,lower.tail=F)
plot(adjP)
data$P_adj <- adjP
gc()

length(which(data$P_unadj==1))
length(which(data$P_unadj==NA))
length(which(is.na(data$P_unadj)))
length(which(is.na(data$P_adj)))

#log(1)
min(data$P_adj)

###################################
## Write out files with p-values ##
###################################
write.csv(data,paste0(path,"/association_cv/allChrom_",pheno,"_gc_adj.csv"),row.names=FALSE)


##################################
##           Plot QQ            ##
##################################


system.time(data<-data[order(data$P_adj),])
#head(data)

#min(data$P_unadj)

#system.time(
#if(debug){
#    pt <- proc.time()
#}
x<-ppoints(dim(data)[1])
png(paste0(path,"/outputs/",pheno,"_qq.png"),width=4,height=4,units="in",res=300)
par(mar=c(4,4,1,1))
plot(-log10(x[data$Type=='IMPU']),-log10(data$P_adj[data$Type=='IMPU']),
     xlim=c(0,-1.05*log10(min(x))),ylim=c(0,-1.05*log10(min(data$P_adj))),
     col=qq_colors[2],pch=16,xlab="",ylab="",main="",xaxs="i",yaxs="i")
points(-log10(x[data$Type=='META']),-log10(data$P_adj[data$Type=='META']),col=qq_colors[3],pch=16)
points(-log10(x[data$Type=='GENO']),-log10(data$P_adj[data$Type=='GENO']),col=qq_colors[1],pch=16)
mtext("-log10(expected p-value)",side=1,line=2.5)
mtext("-log10(observed p-value)",side=2,line=2.5)
abline(a=0,b=1,lty=2)
legend("topleft",legend=c("Genotyped","Imputed","Meta"),
       col=c(qq_colors[1],qq_colors[2],qq_colors[3]),pch=16)
graphics.off()
#)
#if(debug){
#    print(proc.time()-pt)
#}

## Comparison plot
#system.time(
#if(debug){
#    pt <- proc.time()
#}
x<-ppoints(dim(data)[1])
png(paste0(path,"/outputs/",pheno,"_qq.compare.png"),width=4,height=4,units="in",res=300)
par(mar=c(4,4,1,1))
plot(-log10(x[data$Type=='IMPU']),-log10(data$P_adj[data$Type=='IMPU']),
     xlim=c(0,-1.05*log10(min(x))),ylim=c(0,-1.05*log10(min(data$P_adj))),
     col=qq_colors[2],pch=16,xlab="",ylab="",main="",xaxs="i",yaxs="i")
points(-log10(x[data$Type=='META']),-log10(data$P_adj[data$Type=='META']),col=qq_colors[3],pch=16)
points(-log10(x[data$Type=='GENO']),-log10(data$P_adj[data$Type=='GENO']),col=qq_colors[1],pch=16)
points(-log10(x[data$Type=='IMPU']),-log10(data$P_unadj[data$Type=='IMPU']),col=alpha(qq_colors[2],.2),pch=16)
points(-log10(x[data$Type=='META']),-log10(data$P_unadj[data$Type=='META']),col=alpha(qq_colors[3],.2),pch=16)
points(-log10(x[data$Type=='GENO']),-log10(data$P_unadj[data$Type=='GENO']),col=alpha(qq_colors[1],.2),pch=16)
mtext("-log10(expected p-value)",side=1,line=2.5)
mtext("-log10(observed p-value)",side=2,line=2.5)
abline(a=0,b=1,lty=2)
legend("topleft",legend=c("Genotyped","Imputed","Meta"),
       col=c(qq_colors[1],qq_colors[2],qq_colors[3]),pch=16)
graphics.off()
#)
#if(debug){
#    print(proc.time()-pt)
#}

#data <- data[data$CHR<=22,]


system.time(data<-data[order(data$P_unadj),])

## Plot only unadjusted p-val QQ plot
#if(debug){
#    pt <- proc.time()
#}
x<-ppoints(dim(data)[1])
png(paste0(path,"/outputs/",pheno,"_qq.unadjusted.png"),width=4,height=4,units="in",res=300)
par(mar=c(4,4,1,1))
plot(-log10(x[data$Type=='IMPU']),-log10(data$P_unadj[data$Type=='IMPU']),
     xlim=c(0,-1.05*log10(min(x))),ylim=c(0,-1.05*log10(min(data$P_unadj))),
     col=qq_colors[2],pch=16,xlab="",ylab="",main="",xaxs="i",yaxs="i")
points(-log10(x[data$Type=='META']),-log10(data$P_unadj[data$Type=='META']),col=qq_colors[3],pch=16)
points(-log10(x[data$Type=='GENO']),-log10(data$P_unadj[data$Type=='GENO']),col=qq_colors[1],pch=16)
mtext("-log10(expected p-value)",side=1,line=2.5)
mtext("-log10(observed p-value)",side=2,line=2.5)
abline(a=0,b=1,lty=2)
legend("topleft",legend=c("Genotyped","Imputed","Meta"),
       col=c(qq_colors[1],qq_colors[2],qq_colors[3]),pch=16)
graphics.off()
#)
#if(debug){
#    print(proc.time()-pt)
#}

data[which(data$P_unadj==min(data$P_unadj)),]

table(data$CHR)



##################################
##        Plot Manhattan        ##
##################################

#length(which(is.na(data$BP)))
#length(which(is.na(x)))
#i <- 1
#length(which(is.na(x[data$CHR==chr[i]])))
#head(data[which(is.na(x[data$CHR==chr[i]]))])
#data[99646]
#length(x)
#nrow(data)
#length(which(data$CHR==chr[i]))
#head(x[data$CHR==chr[i]]+x_end)
#data[data$CHR==chr[i]]
#tmp <- (which(data$CHR==chr[i]))
#length(tmp)
#tmp1 <- which(is.na(x[data$CHR==chr[i]]+x_end))
#length(tmp1)
#tmp <- data[data$CHR==chr[i]]
#gc()
#tmp[99646]
#tmp[which(is.na(tmp$BP))]
#head(data[is.na(data$BP)])
## Replace allImp "BP" values with genotyped BP values "BP_geno" if exists
#data[!is.na(BP_geno),"BP"] <- data[!is.na(BP_geno),"BP_geno"]
#dim(data[!is.na(BP_geno),"BP"])
#length(which(!(data[!is.na(BP_geno),"BP"] == data[!is.na(BP_geno),"BP_geno"])))

x<-data$BP
x_end<-0;x_mid<-c();x_edge<-0
head(x)
#length(which(data$BP %in% NA))
#summary(data[which(data[,CHR]==chr[i]),BP])
#length(data$BP[data[,CHR]==chr[i]])
                                        #i <- 1
table(data[,CHR])
# <- 4
for (i in 1:length(chr)) {
    print(i)
  x[which(data[,CHR]==chr[i])]<-x[which(data[,CHR]==chr[i])]+x_end
  print(paste0(x_end,' ',x_mid,' ',x_edge))
  temp<-x_end
  x_end<-x_end+max(data[which(data[,CHR]==chr[i]),BP])
  x_mid<-c(x_mid,(temp+x_end)/2)
  x_edge<-c(x_edge,x_end)
  #print(paste0(x_end,x_mid,x_edge))
}

#summary(x)
##head(data)
min(data$P_unadj)
y <- as.vector(-log10(data$P_adj))
y.unadj <- as.vector(-log10(data$P_unadj))
summary(y.unadj)
gc()

#length(y)
#head(y)
#head(x)

## Check top hits ## DEBUG ##
#data[intersect(which(y>8),which(data$CHR==10)),]
#data[intersect(which(y>6),which(data$CHR==13)),]
#data[intersect(which(y>6),which(data$CHR==18)),]
#data[intersect(which(y>6),which(data$CHR==4)),]

#head(data)
# Get cutoffs for Bonferroni corrections on alpha levels 0.1 and 0.05
#cutoff_significant<-0.05/dim(data)[1]
cutoff_significant<-1E-8
#cutoff_suggestive<-cutoff_significant*100
cutoff_suggestive<-5E-6
size<-rep(.5,dim(data)[1])
size[data$P_adj<=cutoff_suggestive]<-.5
size[data$P_adj<=cutoff_significant]<-.5
summary(size)

#summary(

## UNADJUSTED
if(debug){
    pt <- proc.time()
}
png(paste0(path,"/outputs/",pheno,"_manhattan.unadjusted.v2.png"),width=8,height=4,units="in",res=300)
par(mar=c(4,4,2,1))
# Paint background
plot(1:10,xlim=c(0,x_end),ylim=c(0,max(y.unadj)*1.05),
     xlab="",ylab="",main="",xaxt="n",xaxs="i",yaxs="i")
d<-c(1,x_edge[2:length(x_edge)]/x_end)
r<-par("usr")
for (i in 1:length(x_edge)) {
    #rect(r[(i>1)+1]*d[i],r[3],r[2]*d[i+1],r[4],col=manhattan_background_color[i%%2+1],border=NA)
    rect(r[(i>1)+1]*d[i],r[3],r[2]*d[i+1],r[4],col='white',border=NA)
}
# Plot imputed variants 
idx<-which(data$CHR==chr[1] & data$Type=='IMPU')
length(idx)
tmp.cols <- manhattan_colors[3:4]
#points(x[idx],y[idx],col=manhattan_color[1],pch=16,cex=size[idx])
points(x[idx],y.unadj[idx],col=alpha(tmp.cols[1],manhattan_alpha[1]),pch=16,cex=.4)
i <- 3
if (length(chr)>1) {
    for (i in 2:length(chr)) {
        print(i)
        idx<-which(data$CHR==chr[i] & data$Type=='IMPU')
        #length(idx)
        #head(y.unadj[idx])
        points(x[idx],y.unadj[idx],col=alpha(tmp.cols[((i+1)%%2)+1],manhattan_alpha[((i+1)%%2)+1]),pch=16,cex=.4)
  }
}
# Plot META analysisvariants
tmp.cols <- manhattan_colors[5:6]
idx<-which(data$CHR==chr[1] & data$Type=='META')
head(data)
table(data$Type)
                                        #points(x[idx],y[idx],col=qq_colors[3],pch=16,cex=size[idx])
head(y.unadj[idx])
points(x[idx],y.unadj[idx],col=alpha(tmp.cols[1],manhattan_alpha[1]),pch=16,cex=.4)
if (length(chr)>1) {
  for (i in 2:length(chr)) {
    idx<-which(data$CHR==chr[i] & data$Type=='META')
    points(x[idx],y.unadj[idx],col=alpha(tmp.cols[((i+1)%%2)+1],manhattan_alpha[((i+1)%%2)+1]),pch=16,cex=.4)
  }
}
# Plot genotyped variants
tmp.cols <- manhattan_colors[1:2]
idx<-which(data$CHR==chr[1] & data$Type=='GENO')
#points(x[idx],y[idx],col=qq_colors[1],pch=16,cex=size[idx])
points(x[idx],y.unadj[idx],col=alpha(tmp.cols[1],manhattan_alpha[1]),pch=16,cex=.4)
if (length(chr)>1) {
  for (i in 2:length(chr)) {
    idx<-which(data$CHR==chr[i] & data$Type=='GENO')
    points(x[idx],y.unadj[idx],col=alpha(tmp.cols[((i+1)%%2)+1],manhattan_alpha[((i+1)%%2)+1]),pch=16,cex=.4)
  }
}
abline(a=-log10(cutoff_significant),b=0,lty=2,lwd=.8,col="grey60")
abline(a=-log10(cutoff_suggestive),b=0,lty=2,lwd=.8,col="grey69")
mtext("Chromosome",side=1,line=2.5)
mtext("-log10(p-value)",side=2,line=2.5)
axis(1,at=x_mid,tick=F,labels=chr,las=2,cex.axis=0.8)
axis(1,at=x_edge,labels=F)
#legtext <- c("Genotyped","Imputed","Meta")
#xcoords <- c(0, 1000, 3000)
#secondvector <- (1:length(legtext)*30)-1
#textwidths <- xcoords/secondvector # this works for all but the first element
#textwidths[1] <- 0 # so replace element 1 with a finite number (any will do)
textwidths <- c(20,50,20)
legend("top",legend=c("Genotyped","Imputed","Meta"),
       col=c(manhattan_colors[1],manhattan_colors[3],manhattan_colors[5]),
       #text.width=c(2/3*strwidth("Genotyped"),0,-strwidth("Meta")/2),
       #text.width=textwidths,
       pch=16,bg=NULL,xpd=TRUE,horiz=TRUE,inset=-.15,bty="n")
trash<-graphics.off()
if(debug){
    print(proc.time()-pt)
}

#head(x[idx])
#summary(x)




## ADJUSTED
if(debug){
    pt <- proc.time()
}
png(paste0(path,"/outputs/",pheno,"_manhattan.png"),width=8,height=4,units="in",res=300)
par(mar=c(4,4,2,1))
# Paint background
plot(1:10,xlim=c(0,x_end),ylim=c(0,max(y)*1.05),
     xlab="",ylab="",main="",xaxt="n",xaxs="i",yaxs="i")
d<-c(1,x_edge[2:length(x_edge)]/x_end)
r<-par("usr")
for (i in 1:length(x_edge)) {
    #rect(r[(i>1)+1]*d[i],r[3],r[2]*d[i+1],r[4],col=manhattan_background_color[i%%2+1],border=NA)
    rect(r[(i>1)+1]*d[i],r[3],r[2]*d[i+1],r[4],col='white',border=NA)
}
# Plot imputed variants 
idx<-which(data$CHR==chr[1] & data$Type=='IMPU')
length(idx)
tmp.cols <- manhattan_colors[3:4]
#points(x[idx],y[idx],col=manhattan_color[1],pch=16,cex=size[idx])
points(x[idx],y[idx],col=alpha(tmp.cols[1],manhattan_alpha[1]),pch=16,cex=.4)
i <- 3
if (length(chr)>1) {
  for (i in 2:length(chr)) {
    idx<-which(data$CHR==chr[i] & data$Type=='IMPU')
    points(x[idx],y[idx],col=alpha(tmp.cols[((i+1)%%2)+1],manhattan_alpha[((i+1)%%2)+1]),pch=16,cex=.4)
  }
}
# Plot META analysisvariants
tmp.cols <- manhattan_colors[5:6]
idx<-which(data$CHR==chr[1] & data$Type=='META')
#points(x[idx],y[idx],col=qq_colors[3],pch=16,cex=size[idx])
points(x[idx],y[idx],col=alpha(tmp.cols[1],manhattan_alpha[1]),pch=16,cex=.4)
if (length(chr)>1) {
  for (i in 2:length(chr)) {
    idx<-which(data$CHR==chr[i] & data$Type=='META')
    points(x[idx],y[idx],col=alpha(tmp.cols[((i+1)%%2)+1],manhattan_alpha[((i+1)%%2)+1]),pch=16,cex=.4)
  }
}
# Plot genotyped variants
tmp.cols <- manhattan_colors[1:2]
idx<-which(data$CHR==chr[1] & data$Type=='GENO')
#points(x[idx],y[idx],col=qq_colors[1],pch=16,cex=size[idx])
points(x[idx],y[idx],col=alpha(tmp.cols[1],manhattan_alpha[1]),pch=16,cex=.4)
if (length(chr)>1) {
  for (i in 2:length(chr)) {
    idx<-which(data$CHR==chr[i] & data$Type=='GENO')
    points(x[idx],y[idx],col=alpha(tmp.cols[((i+1)%%2)+1],manhattan_alpha[((i+1)%%2)+1]),pch=16,cex=.4)
  }
}
abline(a=-log10(cutoff_significant),b=0,lty=2,lwd=.8,col="grey60")
abline(a=-log10(cutoff_suggestive),b=0,lty=2,lwd=.8,col="grey69")
mtext("Chromosome",side=1,line=2.5)
mtext("-log10(gc adjusted p-value)",side=2,line=2.5)
axis(1,at=x_mid,tick=F,labels=chr,las=2,cex.axis=0.8)
axis(1,at=x_edge,labels=F)
#legtext <- c("Genotyped","Imputed","Meta")
#xcoords <- c(0, 1000, 3000)
#secondvector <- (1:length(legtext)*30)-1
#textwidths <- xcoords/secondvector # this works for all but the first element
#textwidths[1] <- 0 # so replace element 1 with a finite number (any will do)
textwidths <- c(20,50,20)
legend("top",legend=c("Genotyped","Imputed","Meta"),
       col=c(manhattan_colors[1],manhattan_colors[3],manhattan_colors[5]),
       #text.width=c(2/3*strwidth("Genotyped"),0,-strwidth("Meta")/2),
       #text.width=textwidths,
       pch=16,bg=NULL,xpd=TRUE,horiz=TRUE,inset=-.15,bty="n")
trash<-graphics.off()
if(debug){
    print(proc.time()-pt)
}
