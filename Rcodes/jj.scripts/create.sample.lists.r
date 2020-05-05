## Get data
setwd('~/data/accord/data/analysis')

## Set up data matrix for accord analysis
dat <- read.table(file='pheno_data_hgi.txt',header=T)
dat <- dat[,c(1,2,4,5,21)]
#colnames(dat)
#head(dat)

## what are 'allraces'? ## DEBUG ##
#d <- read.table(file='GV/GV_allraces/pheno_data/sample_list.txt',header=T)
#dim(d)
#tmp <- setdiff(dat$MaskID,d[,1])
#length(tmp)
#table(dat$ethnicity[dat$MaskID %in% tmp])
#length(which(is.na(dat$ethnicity[dat$MaskID %in% tmp])))
## Based on this output, it looks like all races includes all races 

## Cycle arm and race and output files
unique(dat$ethnicity)
length(which(is.na(dat$ethnicity)))
dim(dat)
table(dat$ethnicity)

races <- ('allraces','Black','White','Hispanic')
arms <- unique(dat$arm)

write.table(dat,file='pheno_data_GV.txt',col.names=T,row.names=FALSE,quote=F)


dat <- read.table(file='../../pheno_data.txt',header=T)
dat <- dat[,c(1,2,4,5,21)]

