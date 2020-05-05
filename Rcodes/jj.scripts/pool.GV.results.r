rm(list=ls())
setwd('~/data/accord/data/analysis')
start <- TRUE
path <- 'GV'
dset <- dir(path)[c(1,18)]
for (i in dset){
    print(file.path(path,i))
    if(start){
        d <- data.frame(Stratification=dset,stringsAsFactors=FALSE)
        thecols <- gsub('chr3.35.','',gsub('.assoc','',list.files(file.path(path,i,'association_cv/imputed_chunks/'),pattern='chr3.35.*')))
        thecols <- rep(thecols,each=4)
        thecols[seq(1,length(thecols),4)] <- paste0(thecols[seq(1,length(thecols),4)],rep('_P',length(thecols)/4))
        thecols[seq(2,length(thecols),4)] <- paste0(thecols[seq(2,length(thecols),4)],rep('_adjP',length(thecols)/4))
        thecols[seq(3,length(thecols),4)] <- paste0(thecols[seq(3,length(thecols),4)],rep('_BETA',length(thecols)/4))
        thecols[seq(4,length(thecols),4)] <- paste0(thecols[seq(4,length(thecols),4)],rep('_SE',length(thecols)/4))
        d[,thecols] <- NA
        start <- FALSE
    }
    list <- list.files(file.path(path,i,'association_cv/imputed_chunks/'),pattern='chr3.35.*')
    l <- list[1]
    for (l in list){
        acol <- gsub('chr3.35.','',gsub('.assoc','',l))
        tmp.d <- read.table(file.path(path,i,'association_cv/imputed_chunks/',l),header=T)
        d[d[,1]==i,colnames(d)==paste0(acol,'_P')] <- tmp.d$P
        d[d[,1]==i,colnames(d)==paste0(acol,'_adjP')] <- NA
        d[d[,1]==i,colnames(d)==paste0(acol,'_BETA')] <- tmp.d$BETA
        d[d[,1]==i,colnames(d)==paste0(acol,'_SE')] <- tmp.d$SE
    }    
}

## Get min p per start
d$lowestP <- apply(d[,grep('_P',colnames(d))],1,min)

dim(d[,grep('_P',colnames(d))])
bonf <- .05/(ncol(d[,grep('_P',colnames(d))])*nrow(d[,grep('_P',colnames(d))]))
d[which(d$lowestP<bonf),]


#library(pheatmap)

#pheatmap(-log10(d[,grep('_P',colnames(d))]))


#as.vector(unlist(d[,grep('_P',colnames(d))]))

#dim(d[,grep('_P',colnames(d))])

tmp.ps <- as.vector(unlist(d[,grep('_P',colnames(d))]))
length(tmp.ps)
## 24 tests

tmp.qs.bh <- p.adjust(tmp.ps,"BH")
min(tmp.qs.bh)
tmp.qs.fdr <- p.adjust(tmp.ps,"fdr")
min(tmp.qs.fdr)



d[,grep('_adjP',colnames(d))] <- matrix(tmp.qs.bh,nrow=dim(d[,grep('_adjP',colnames(d))])[1],ncol=dim(d[,grep('_P',colnames(d))])[2])
d[,grep('_P',colnames(d))]
d[,sort(c(1,grep('_P',colnames(d)),grep('_adjP',colnames(d))))]

write.table(file='GV.table.pvalues.only.tsv',d[,sort(c(1,grep('_P',colnames(d)),grep('_adjP',colnames(d))))],col.names=T,row.names=FALSE,quote=F,sep='\t')

write.table(file='GV.table.tsv',d,col.names=T,row.names=FALSE,quote=F,sep='\t')

pheatmap(-log10(d[,grep('_Q',colnames(d))]))



d <- read.table(file='pheno_data.txt',header=T)
dim(d[d$ethnicity=='White',])
dim(d)
table(d$censor_po)
table(d$censor_tm)
table(d$censor_cm)
table(d$censor_nmi)
table(d$censor_nst)
table(d$censor_tst)
table(d$censor_chf)
table(d$censor_ex)
table(d$censor_maj)

length(na.omit(d$GV_mean))
median(na.omit(d$GV_mean))
length(na.omit(d$GV_sd))

d <- read.table(file='pheno_data.txt',header=T)
dim(d)
head(d)
table(d$censor_po[d$ethnicity=='White'])
table(d$censor_tm[d$ethnicity=='White'])
table(d$censor_cm[d$ethnicity=='White'])
table(d$censor_nmi[d$ethnicity=='White'])
table(d$censor_nst[d$ethnicity=='White'])
table(d$censor_tst[d$ethnicity=='White'])
table(d$censor_chf[d$ethnicity=='White'])
table(d$censor_ex[d$ethnicity=='White'])
table(d$censor_maj[d$ethnicity=='White'])

length(na.omit(d$GV_mean[d$ethnicity=='White']))
median(na.omit(d$GV_mean))
median(na.omit(d$GV_mean[d$ethnicity=='White']))

length(na.omit(d$GV_sd[d$ethnicity=='White']))
