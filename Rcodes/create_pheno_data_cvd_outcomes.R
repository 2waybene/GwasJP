setwd("~/Research/kirk/NCSU_Projects/Accord/Data/")

library(data.table)
cvd <- fread("pheno_data/cvdoutcomes.txt",data.table=FALSE)

out.cols <- grep("censor",colnames(cvd))
cvd.out <- cvd[,c(1,out.cols)]
swap <- function(x) ifelse(x==1,0,1)
cvd.temp <- apply(cvd.out[,2:ncol(cvd.out)],2,swap)
cvd.out <- cbind(cvd.out[,1],cvd.temp)
colnames(cvd.out)[1] <- "MaskID"

pheno.data <- fread("analysis/pheno_data_dtw_tmExtreme_k3.txt",data.table=FALSE)

pheno <- merge(pheno.data,cvd.out,by="MaskID")
pheno <- pheno[,-grep("Extremes",colnames(pheno))]

write.table(pheno,"analysis/pheno_data_cvd_outcomes.txt",row.names=FALSE,quote=FALSE)
