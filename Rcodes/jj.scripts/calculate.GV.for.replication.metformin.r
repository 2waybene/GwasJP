rm(list=ls())
## Calculate new phenotype, Glycaemic Variability (GV),
## based on text from "Request for replication betweeen rs8192675 and glycaemic variability:
## Definition of GV: Please use the standard error (SD) of ALL HbA1c values as your measure of GV. We suggest a minimum of 3 HbA1c values to be used as sample exclusion criteria. We acknowledge there are an array of different GV measures, of which SD and coefficient of variability (CV) are most commonly used, extended from metrics used in continuous glucose monitoring. Our preliminary heritability analysis suggested SD and CV are genetically identical (rg=0.99), therefore equivalent in variant discovery. The utility of other less commonly used measures of GV for genetic studies should be examined once larger data sets are available.

## I am assuming that standard error (SD) in the above is supposed to be "standard deviation (SD)" and
## "coefficient of variability (CV)" means "coefficient of variation (CV)"

## If length(HbA1c) > 3
## GV(sample) = SD(sample)

## GV = glycaemic varibility (SD)
## meanHba1c = mean of Hba1c
## cvHba1c = coefficient of variation of Hba1c

## Get data
setwd('~/data/accord/data/analysis')

#################################################
##### READ IN hba1c for GV calculations #########
#################################################
phe <- read.table(file='../pheno_data/hba1c.txt',header=T)
ind.list <- unique(phe$MaskID)
#length(ind.list)
#head(phe)

## Set up data matrix for GV, mean, and coefficient of variation
d <- data.frame(MaskID=ind.list,GV_sd=rep(999999999,length(ind.list)),GV_mean=rep(999999999,length(ind.list)),GV_cv=rep(999999999,length(ind.list)))
for (i in 1:nrow(d)){
    vals <- phe$hba1c[which(phe$MaskID %in% d$MaskID[i])]
    if (length(vals)>3){
        d$GV_sd[i] <- sd(vals)
        d$GV_mean[i] <- mean(vals)
        d$GV_cv[i] <- sd(vals)/mean(vals)
    }else{ ## If there are less than three values, do not use GV (set to NA)
        d$GV_sd[i] <- NA
        d$GV_mean[i] <- NA
        d$GV_cv[i] <- NA
    }
}

head(d)

save.image('../pheno_data/calculate.GV.RData')

########## START HERE ##########
rm(list=ls())
load('../pheno_data/calculate.GV.RData')
head(d)




###################################
##### READ IN pheno_data file #####
###################################

## Set up data matrix for accord analysis
dat <- read.table(file='pheno_data_metformin.txt',header=T)
#head(dat)

keep.cols <- c("MaskID","LabID","baseline_age","gender","arm","network","edu","smoking","alcohol","yrsdiab","hyptens","cvd_hx_baseline","x3malb","x3lvh","x3sten","eyedisea","neuropat","bmi","waist_cm","yrslipi","ethnicity","int_gly_arm","int_bp_arm","fib_arm")
dat <- dat[,which(colnames(dat)%in%keep.cols)]
colnames(dat)

#colnames(dat)
colnames(dat[c(1:11,18,20,21)]) ## only want to keep what is in 'starting_covars.txt'
dat <- dat[,c(1:11,18,20,21)]
dim(dat) 
dat <- merge(dat,d,by.x="MaskID",all.x=TRUE)
head(dat)
dim(dat)

#########################################
######## READ IN CVD outcomes ###########
#########################################
cvds <- read.table(file='../pheno_data/cvdoutcomes.txt',header=T,skip=33,sep='\t')
#head(cvds)
#table(cvds$fuyrs_po7p)
#length(which(is.na(cvds$fuyrs_po7p)))
colnames(cvds)[c(1,grep('censor',colnames(cvds)))]

cvds <- cvds[,c(1,grep('censor',colnames(cvds)))]
#head(cvds)

i <- 1
## Switch 0s and 1s for outcomes
for (i in 2:ncol(cvds)){
    ## Check NAs
    print(length(which(is.na(cvds[,i]))))
    print(table(cvds[,i]))
    cvds[which(cvds[,i]==1),i] <- -1
    cvds[which(cvds[,i]==0),i] <- 1
    cvds[which(cvds[,i]==-1),i] <- 0
    print(table(cvds[,i]))
}



dat <- merge(dat,cvds,by.x="MaskID",all.x=TRUE)
head(dat)
dim(dat)

## Make sure _GV corresponds to the path name when you run the pipeline
write.table(dat,file='pheno_data_GV.txt',col.names=T,row.names=FALSE,quote=F)


