# Script for caitrin.
## Input: ACCORD data
## Output: EasyQC file

rm(list=ls())
## Load data
ACCORD <- read.table ("~/downloads/allChrom_ICAPS_po_gc_adj.csv", header=TRUE, sep =",")
summary(ACCORD$OR_allImp)
ACCORD_1 <- subset(ACCORD, select=c(CHR, BP, A1_allImp, A2_allImp, FRQ_allImp, OR_allImp, SE_allImp, P_allImp))
ACCORD_1$SNPID2 <- paste(ACCORD_1$CHR,ACCORD_1$BP,sep=":")
ACCORD_1$n <- rep(1645,nrow(ACCORD_1))
ACCORD_1$Rsq <- rep(0.5,nrow(ACCORD_1))
head(ACCORD_1, n=10)

###############################################
## Addition by John Jack
## NEW CODE: Convert OR and SEs to Beta
ACCORD_1$OR_allImp = log(ACCORD_1$OR_allImp)
ACCORD_1$SE_allImp = log(ACCORD_1$SE_allImp)
###############################################

## The next line was throwing an error. Column "snp" did not exist.
#colnames(ACCORD_1) <-c("Chr", "Position", "Coded_all", "Noncoded_all", "AF_coded_all", "Beta", "SE", "Pvalue", "snp", "n_total", "Rsq", "SNPID")
colnames(ACCORD_1) <-c("Chr", "Position", "Coded_all", "Noncoded_all", "AF_coded_all", "Beta", "SE", "Pvalue", "SNPID", "n_total", "Rsq" )
ACCORD_2 <- ACCORD_1[c("SNPID", "Chr", "Position", "Coded_all", "Noncoded_all", "n_total", "AF_coded_all", "Beta", "SE", "Pvalue", "Rsq")]
head(ACCORD_2)
write.table(ACCORD_2, "ACCPRD_1000G_CCB_1_22_imputed_easyQC.txt", quote=FALSE, row.names = F, sep= "\t")


###############################################
## Additon by John Jack.
## Calculate Z score to recalculate P values from betas (as does EasyQC for the PZ plot)
## Z score is Betas/SEs
ACCORD_2$Z=ACCORD_2$Beta/ACCORD_2$SE
## Calculate P value from Z score
ACCORD_2$PZ=2*pnorm(-abs(ACCORD_2$Z))
head(ACCORD_2)
png('~/Desktop/caitrin.plot.png',height=1000,width=1000)
plot(ACCORD_2$Pvalue,ACCORD_2$PZ,xlab='P value',ylab='PZ value')
graphics.off()
##############################################