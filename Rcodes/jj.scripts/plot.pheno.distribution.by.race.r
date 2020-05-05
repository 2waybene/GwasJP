rm(list=ls())
setwd('~/data/accord/data/analysis/')   

tmp <- read.table(file='hgi_baseline/hgi_allraces/pheno_data/pheno_hgi.txt',header=T,stringsAsFactors=FALSE)
head(tmp)
tmp$Race='AllRaces'
length(unique(tmp$FID))
nrow(tmp)

#
d <- tmp
dim(d)

tmp <- read.table(file='hgi_baseline/hgi_white/pheno_data/pheno_hgi.txt',header=T,stringsAsFactors=FALSE)
head(tmp)
tmp$Race='White'
head(tmp)
tail(tmp)
d <- rbind(d,tmp)

tmp <- read.table(file='hgi_baseline/hgi_black/pheno_data/pheno_hgi.txt',header=T,stringsAsFactors=FALSE)
head(tmp)
tmp$Race='Black'
head(tmp)
tail(tmp)
d <- rbind(d,tmp)

tmp <- read.table(file='hgi_baseline/hgi_hispanic/pheno_data/pheno_hgi.txt',header=T,stringsAsFactors=FALSE)
head(tmp)
tmp$Race='Hispanic'
head(tmp)
tail(tmp)
d <- rbind(d,tmp)

png('hgi_baseline/Boxplot.By.Race.png',height=500,width=500)
par(mar=c(6,8,1,1))
boxplot(d$hgi~d$Race,xlab='HGI',cex.axis=1.5,cex.lab=2,horizontal=T,las=1)
graphics.off()

length()

png('hgi_baseline/Histogram.By.Race.png',height=1000,width=1000)
par(mfrow=c(2,2))
hist(d$hgi[d$Race=='AllRaces'],xlab='HGI',main=paste0('All Races\n(n=',length(d$hgi[d$Race=='AllRaces']),')'),xlim=c(-4,4),cex.lab=1.5,cex.axis=1.2,cex.main=2)
hist(d$hgi[d$Race=='White'],xlab='HGI',main=paste0('White\n(n=',length(d$hgi[d$Race=='White']),')'),xlim=c(-4,4),cex.lab=1.5,cex.axis=1.2,cex.main=2)
hist(d$hgi[d$Race=='Black'],xlab='HGI',main=paste0('Black\n(n=',length(d$hgi[d$Race=='Black']),')'),xlim=c(-4,4),cex.lab=1.5,cex.axis=1.2,cex.main=2)
hist(d$hgi[d$Race=='Hispanic'],xlab='HGI',main=paste0('Hispanic\n(n=',length(d$hgi[d$Race=='Hispanic']),')'),xlim=c(-4,4),cex.lab=1.5,cex.axis=1.2,cex.main=2)
graphics.off()


# PLOT HBA1C and HGI
d[1,]
tmp <- read.table(file='pheno_data.txt',header=T,stringsAsFactors=FALSE)
tmp[1,]
for (i in 1:nrow(d)){
    
}


rm(list=ls())
setwd('~/data/accord/data/analysis/')   
tmp <- read.table(file='pheno_data.txt',header=T,stringsAsFactors=FALSE)

head(tmp)
tmp$Race='AllRaces'
length(unique(tmp$FID))
nrow(tmp)
d <- tmp
dim(d)

tmp <- read.table(file='hgi_baseline/hgi_white/pheno_data/pheno_hgi.txt',header=T,stringsAsFactors=FALSE)
d$Race[which(d$MaskID %in% tmp$FID)] <- 'White'
tmp <- read.table(file='hgi_baseline/hgi_black/pheno_data/pheno_hgi.txt',header=T,stringsAsFactors=FALSE)
d$Race[which(d$MaskID %in% tmp$FID)] <- 'Black'
tmp <- read.table(file='hgi_baseline/hgi_hispanic/pheno_data/pheno_hgi.txt',header=T,stringsAsFactors=FALSE)
d$Race[which(d$MaskID %in% tmp$FID)] <- 'Hispanic'

dim(d)
d[1,]

png('hgi_baseline/Histogram.By.Race.HBA1C.png',height=1000,width=1000)
par(mfrow=c(2,2))
hist(d$hba1c[d$Race=='AllRaces'],xlab='hba1c',main=paste0('Other\n(n=',length(d$hba1c[d$Race=='AllRaces']),')'),cex.lab=1.5,cex.axis=1.2,cex.main=2)
hist(d$hba1c[d$Race=='White'],xlab='hba1c',main=paste0('White\n(n=',length(d$hba1c[d$Race=='White']),')'),cex.lab=1.5,cex.axis=1.2,cex.main=2)
hist(d$hba1c[d$Race=='Black'],xlab='hba1c',main=paste0('Black\n(n=',length(d$hba1c[d$Race=='Black']),')'),cex.lab=1.5,cex.axis=1.2,cex.main=2)
hist(d$hba1c[d$Race=='Hispanic'],xlab='hba1c',main=paste0('Hispanic\n(n=',length(d$hba1c[d$Race=='Hispanic']),')'),cex.lab=1.5,cex.axis=1.2,cex.main=2)
graphics.off()

png('hgi_baseline/Histogram.By.Race.FPG.png',height=1000,width=1000)
par(mfrow=c(2,2))
hist(d$fpg[d$Race=='AllRaces'],xlab='HGI',main=paste0('Other\n(n=',length(d$fpg[d$Race=='AllRaces']),')'),cex.lab=1.5,cex.axis=1.2,cex.main=2)
hist(d$fpg[d$Race=='White'],xlab='HGI',main=paste0('White\n(n=',length(d$fpg[d$Race=='White']),')'),cex.lab=1.5,cex.axis=1.2,cex.main=2)
hist(d$fpg[d$Race=='Black'],xlab='HGI',main=paste0('Black\n(n=',length(d$fpg[d$Race=='Black']),')'),cex.lab=1.5,cex.axis=1.2,cex.main=2)
hist(d$fpg[d$Race=='Hispanic'],xlab='HGI',main=paste0('Hispanic\n(n=',length(d$fpg[d$Race=='Hispanic']),')'),cex.lab=1.5,cex.axis=1.2,cex.main=2)
graphics.off()

png('hgi_baseline/Boxplot.By.Race.HBA1C.png',height=500,width=500)
par(mar=c(6,8,1,1))
boxplot(d$hba1c~d$Race,xlab='hba1c',cex.axis=1.5,cex.lab=2,horizontal=T,las=1)
graphics.off()

png('hgi_baseline/Boxplot.By.Race.FPG.png',height=500,width=500)
par(mar=c(6,8,1,1))
boxplot(d$fpg~d$Race,xlab='fpg',cex.axis=1.5,cex.lab=2,horizontal=T,las=1)
graphics.off()
