rm(list=ls())

###################
#### All Races ####

##args<-(commandArgs(TRUE))
#path<-args[[1]]
setwd('~/data/accord/data/analysis/')
path <- 'hgi_baseline/hgi_allraces'
getwd()
# Load full phenotype data file
data<-read.table("pheno_data.txt",header=T,stringsAsFactors=F)
#head(data)
#summary(data)
#data$hba1c
# Retain selected samples
if (file.exists(paste0(path,"/starting_samples.txt"))) {
    keep<-read.table(paste0(path,"/starting_samples.txt"),stringsAsFactors=F)
    data<-data[data$MaskID%in%keep[,1] & data$LabID%in%keep[,2],]
}
# Retain selected variables, force BLR values into model if looking at 4-month difference
keep<-unlist(read.table(paste0(path,"/starting_covars.txt"),stringsAsFactors=F))
pheno<-unlist(read.table(paste0(path,"/phenotypes.txt"),stringsAsFactors=F))
pheno_blr<-c()
for (i in 1:length(pheno)) {
    if (grepl("d4m",pheno[i])) {
        pheno_blr<-c(pheno_blr,strsplit(pheno[i],'_')[[1]][2])
    }
}
data<-data[,names(data)%in%c(keep,pheno,pheno_blr,"MaskID","LabID")]
#plot(data$hba1c ~ data$fpg)
#plot()
# Hemoglobin glycation index (HGI)
if ("hgi"%in%c(keep,pheno)) {
    # Compute HGI
    fit<-lm(data$hba1c ~ data$fpg, na.action=na.exclude)
    hgi<-round(residuals(fit),digits=4)
    data<-cbind(data,hgi)
    # Split into 3 groups
    grp<-rep(NA,dim(data)[1])
    cutoffs<-quantile(residuals(fit),c(1/3,2/3),na.rm=T)
    grp[residuals(fit)<cutoffs[1]]<-1
    grp[residuals(fit)>=cutoffs[1] & residuals(fit)<cutoffs[2]]<-2
    grp[residuals(fit)>=cutoffs[2]]<-3
    grp<-factor(grp)
    # Plot HGI
    pdf(file=paste0(path,"/outputs/HbA1c_vs_FPG.pdf"),width=10, height=7.5)
    plot(data$fpg,data$hba1c,col=grp,pch=16,cex=0.5,xlab="FPG",ylab="HbA1c",main="HbA1c vs FPG")
    abline(a=coef(fit)[1],b=coef(fit)[2],col=4,lwd=2)
    eq<-paste("y = ",format(coef(fit)[1], digits = 4)," + ",format(coef(fit)[2], digits = 4)," * x",sep="")
    text(max(data$fpg,na.rm=T),min(data$hba1c,na.rm=T),label=eq,col=4,adj=1)
    trash<-dev.off()
}
## Save all races fpg and hba1c data
d <- data
all.grp <- grp
all.fit <- fit
d$Race <- 'Other'



###############
#### White ####
###############

#args<-(commandArgs(TRUE))
#path<-args[[1]]
setwd('~/data/accord/data/analysis/')
path <- 'hgi_baseline/hgi_white'
getwd()
# Load full phenotype data file
data<-read.table("pheno_data.txt",header=T,stringsAsFactors=F)
#head(data)
#summary(data)
#data$hba1c
# Retain selected samples
if (file.exists(paste0(path,"/starting_samples.txt"))) {
    keep<-read.table(paste0(path,"/starting_samples.txt"),stringsAsFactors=F)
    data<-data[data$MaskID%in%keep[,1] & data$LabID%in%keep[,2],]
}
# Retain selected variables, force BLR values into model if looking at 4-month difference
keep<-unlist(read.table(paste0(path,"/starting_covars.txt"),stringsAsFactors=F))
pheno<-unlist(read.table(paste0(path,"/phenotypes.txt"),stringsAsFactors=F))
pheno_blr<-c()
for (i in 1:length(pheno)) {
    if (grepl("d4m",pheno[i])) {
        pheno_blr<-c(pheno_blr,strsplit(pheno[i],'_')[[1]][2])
    }
}
data<-data[,names(data)%in%c(keep,pheno,pheno_blr,"MaskID","LabID")]
#plot(data$hba1c ~ data$fpg)
#plot()
# Hemoglobin glycation index (HGI)
if ("hgi"%in%c(keep,pheno)) {
    # Compute HGI
    fit<-lm(data$hba1c ~ data$fpg, na.action=na.exclude)
    hgi<-round(residuals(fit),digits=4)
    data<-cbind(data,hgi)
    # Split into 3 groups
    grp<-rep(NA,dim(data)[1])
    cutoffs<-quantile(residuals(fit),c(1/3,2/3),na.rm=T)
    grp[residuals(fit)<cutoffs[1]]<-1
    grp[residuals(fit)>=cutoffs[1] & residuals(fit)<cutoffs[2]]<-2
    grp[residuals(fit)>=cutoffs[2]]<-3
    grp<-factor(grp)
    # Plot HGI
    pdf(file=paste0(path,"/outputs/HbA1c_vs_FPG.pdf"),width=10, height=7.5)
    plot(data$fpg,data$hba1c,col=grp,pch=16,cex=0.5,xlab="FPG",ylab="HbA1c",main="HbA1c vs FPG")
    abline(a=coef(fit)[1],b=coef(fit)[2],col=4,lwd=2)
    eq<-paste("y = ",format(coef(fit)[1], digits = 4)," + ",format(coef(fit)[2], digits = 4)," * x",sep="")
    text(max(data$fpg,na.rm=T),min(data$hba1c,na.rm=T),label=eq,col=4,adj=1)
    trash<-dev.off()
}
## Save all races fpg and hba1c data
white.d <- data
white.grp <- grp
white.fit <- fit
## Debug ## Make sure we're only replacing 'Other' Races
print(unique(d$Race[which(d$MaskID %in% data$MaskID)]))
## Save White race information
d$Race[which(d$MaskID %in% data$MaskID)] <- 'White'

###############
#### Black ####
###############

#args<-(commandArgs(TRUE))
#path<-args[[1]]
setwd('~/data/accord/data/analysis/')
path <- 'hgi_baseline/hgi_black'
getwd()
# Load full phenotype data file
data<-read.table("pheno_data.txt",header=T,stringsAsFactors=F)
#head(data)
#summary(data)
#data$hba1c
# Retain selected samples
if (file.exists(paste0(path,"/starting_samples.txt"))) {
    keep<-read.table(paste0(path,"/starting_samples.txt"),stringsAsFactors=F)
    data<-data[data$MaskID%in%keep[,1] & data$LabID%in%keep[,2],]
}
# Retain selected variables, force BLR values into model if looking at 4-month difference
keep<-unlist(read.table(paste0(path,"/starting_covars.txt"),stringsAsFactors=F))
pheno<-unlist(read.table(paste0(path,"/phenotypes.txt"),stringsAsFactors=F))
pheno_blr<-c()
for (i in 1:length(pheno)) {
    if (grepl("d4m",pheno[i])) {
        pheno_blr<-c(pheno_blr,strsplit(pheno[i],'_')[[1]][2])
    }
}
data<-data[,names(data)%in%c(keep,pheno,pheno_blr,"MaskID","LabID")]
#plot(data$hba1c ~ data$fpg)
#plot()
# Hemoglobin glycation index (HGI)
if ("hgi"%in%c(keep,pheno)) {
    # Compute HGI
    fit<-lm(data$hba1c ~ data$fpg, na.action=na.exclude)
    hgi<-round(residuals(fit),digits=4)
    data<-cbind(data,hgi)
    # Split into 3 groups
    grp<-rep(NA,dim(data)[1])
    cutoffs<-quantile(residuals(fit),c(1/3,2/3),na.rm=T)
    grp[residuals(fit)<cutoffs[1]]<-1
    grp[residuals(fit)>=cutoffs[1] & residuals(fit)<cutoffs[2]]<-2
    grp[residuals(fit)>=cutoffs[2]]<-3
    grp<-factor(grp)
    # Plot HGI
    pdf(file=paste0(path,"/outputs/HbA1c_vs_FPG.pdf"),width=10, height=7.5)
    plot(data$fpg,data$hba1c,col=grp,pch=16,cex=0.5,xlab="FPG",ylab="HbA1c",main="HbA1c vs FPG")
    abline(a=coef(fit)[1],b=coef(fit)[2],col=4,lwd=2)
    eq<-paste("y = ",format(coef(fit)[1], digits = 4)," + ",format(coef(fit)[2], digits = 4)," * x",sep="")
    text(max(data$fpg,na.rm=T),min(data$hba1c,na.rm=T),label=eq,col=4,adj=1)
    trash<-dev.off()
}
## Save all races fpg and hba1c data
black.d <- data
black.grp <- grp
black.fit <- fit
## Debug ## Make sure we're only replacing 'Other' Races
print(unique(d$Race[which(d$MaskID %in% data$MaskID)]))
## Save White race information
d$Race[which(d$MaskID %in% data$MaskID)] <- 'Black'


##################
#### Hispanic ####
##################

#args<-(commandArgs(TRUE))
#path<-args[[1]]
setwd('~/data/accord/data/analysis/')
path <- 'hgi_baseline/hgi_hispanic'
getwd()
# Load full phenotype data file
data<-read.table("pheno_data.txt",header=T,stringsAsFactors=F)
#head(data)
#summary(data)
#data$hba1c
# Retain selected samples
if (file.exists(paste0(path,"/starting_samples.txt"))) {
    keep<-read.table(paste0(path,"/starting_samples.txt"),stringsAsFactors=F)
    data<-data[data$MaskID%in%keep[,1] & data$LabID%in%keep[,2],]
}
# Retain selected variables, force BLR values into model if looking at 4-month difference
keep<-unlist(read.table(paste0(path,"/starting_covars.txt"),stringsAsFactors=F))
pheno<-unlist(read.table(paste0(path,"/phenotypes.txt"),stringsAsFactors=F))
pheno_blr<-c()
for (i in 1:length(pheno)) {
    if (grepl("d4m",pheno[i])) {
        pheno_blr<-c(pheno_blr,strsplit(pheno[i],'_')[[1]][2])
    }
}
data<-data[,names(data)%in%c(keep,pheno,pheno_blr,"MaskID","LabID")]
#plot(data$hba1c ~ data$fpg)
#plot()
# Hemoglobin glycation index (HGI)
if ("hgi"%in%c(keep,pheno)) {
    # Compute HGI
    fit<-lm(data$hba1c ~ data$fpg, na.action=na.exclude)
    hgi<-round(residuals(fit),digits=4)
    data<-cbind(data,hgi)
    # Split into 3 groups
    grp<-rep(NA,dim(data)[1])
    cutoffs<-quantile(residuals(fit),c(1/3,2/3),na.rm=T)
    grp[residuals(fit)<cutoffs[1]]<-1
    grp[residuals(fit)>=cutoffs[1] & residuals(fit)<cutoffs[2]]<-2
    grp[residuals(fit)>=cutoffs[2]]<-3
    grp<-factor(grp)
    # Plot HGI
    pdf(file=paste0(path,"/outputs/HbA1c_vs_FPG.pdf"),width=10, height=7.5)
    plot(data$fpg,data$hba1c,col=grp,pch=16,cex=0.5,xlab="FPG",ylab="HbA1c",main="HbA1c vs FPG")
    abline(a=coef(fit)[1],b=coef(fit)[2],col=4,lwd=2)
    eq<-paste("y = ",format(coef(fit)[1], digits = 4)," + ",format(coef(fit)[2], digits = 4)," * x",sep="")
    text(max(data$fpg,na.rm=T),min(data$hba1c,na.rm=T),label=eq,col=4,adj=1)
    trash<-dev.off()
}
## Save all races fpg and hba1c data
hisp.d <- data
hisp.grp <- grp
hisp.fit <- fit
## Debug ## Make sure we're only replacing 'Other' Races
print(unique(d$Race[which(d$MaskID %in% data$MaskID)]))
## Save White race information
d$Race[which(d$MaskID %in% data$MaskID)] <- 'Hispanic'


grp <- d$Race
## Set races
grp <- factor(grp)
table(grp)
table(as.numeric(grp))
## Black = 1, Hipanic = 2, Other = 3, White = 4

## white.fit
## Plot all
pdf(file=paste0('hgi_baseline/All.Regressions.HbA1c_vs_FPG.pdf'),width=10, height=7.5)
## All races data plot
plot(d$fpg,d$hba1c,col=grp,pch=16,cex=0.5,xlab="FPG",ylab="HbA1c",main="HbA1c vs FPG")
##All races line and eq
abline(a=coef(all.fit)[1],b=coef(all.fit)[2],col=5,lwd=2)
eq<-paste("y = ",format(coef(all.fit)[1], digits = 4)," + ",format(coef(all.fit)[2], digits = 4)," * x",sep="")
text(max(data$fpg,na.rm=T),min(data$hba1c,na.rm=T),label=eq,col=5,adj=c(1,0))
text(max(data$fpg,na.rm=T),min(data$hba1c,na.rm=T),label='All',col=5,adj=c(11.5,0))
## White
tmp.col <- unique(as.numeric(grp[which(grp %in% 'White')]))
abline(a=coef(white.fit)[1],b=coef(white.fit)[2],col=tmp.col,lwd=2)
eq<-paste("y = ",format(coef(white.fit)[1], digits = 4)," + ",format(coef(white.fit)[2], digits = 4)," * x",sep="")
text(max(data$fpg,na.rm=T),min(data$hba1c,na.rm=T),label=eq,col=tmp.col,adj=c(1,2))
text(max(data$fpg,na.rm=T),min(data$hba1c,na.rm=T),label='White',col=tmp.col,adj=c(5.5,2))
## Black
tmp.col <- unique(as.numeric(grp[which(grp %in% 'Black')]))
abline(a=coef(black.fit)[1],b=coef(black.fit)[2],col=tmp.col,lwd=2)
eq<-paste("y = ",format(coef(black.fit)[1], digits = 4)," + ",format(coef(black.fit)[2], digits = 4)," * x",sep="")
text(max(data$fpg,na.rm=T),min(data$hba1c,na.rm=T),label=eq,col=tmp.col,adj=c(1,4))
text(max(data$fpg,na.rm=T),min(data$hba1c,na.rm=T),label='Black',col=tmp.col,adj=c(5.7,4))
## Hispanic
tmp.col <- unique(as.numeric(grp[which(grp %in% 'Hispanic')]))
abline(a=coef(hisp.fit)[1],b=coef(hisp.fit)[2],col=tmp.col,lwd=2)
eq<-paste("y = ",format(coef(hisp.fit)[1], digits = 4)," + ",format(coef(hisp.fit)[2], digits = 4)," * x",sep="")
text(max(data$fpg,na.rm=T),min(data$hba1c,na.rm=T),label=eq,col=tmp.col,adj=c(1,6))
text(max(data$fpg,na.rm=T),min(data$hba1c,na.rm=T),label='Hispanic',col=tmp.col,adj=c(4,6))
trash<-dev.off()
