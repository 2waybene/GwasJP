rm(list=ls())
#for eash qc for ACCORD
ACCORD <- read.table ("~/downloads/allChrom_ICAPS_po_gc_adj.csv", header=TRUE, sep =",")
dim(ACCORD)
ACCORD_1 <- subset(ACCORD, select=c(CHR, BP, A1_allImp, A2_allImp, FRQ_allImp, OR_allImp, SE_allImp, P_allImp))
dim(ACCORD_1)
head(ACCORD_1)
ACCORD_1$SNPID2 <- paste(ACCORD_1$CHR,ACCORD_1$BP,sep=":")
dim(ACCORD_1)
head(ACCORD_1)
ACCORD_1$n <- rep(1645,nrow(ACCORD_1))
ACCORD_1$Rsq <- rep(0.5,nrow(ACCORD_1))

head(ACCORD_1, n=10)

#colnames(ACCORD_1) <-c("Chr", "Position", "Coded_all", "Noncoded_all", "AF_coded_all", "Beta", "SE", "Pvalue", "snp", "n_total", "Rsq", "SNPID")
colnames(ACCORD_1) <-c("Chr", "Position", "Coded_all", "Noncoded_all", "AF_coded_all", "Beta", "SE", "Pvalue", "SNPID", "n_total", "Rsq" )
ACCORD_2 <- ACCORD_1[c("SNPID", "Chr", "Position", "Coded_all", "Noncoded_all", "n_total", "AF_coded_all", "Beta", "SE", "Pvalue", "Rsq")]
head(ACCORD_2)
write.table(ACCORD_2, "ACCPRD_1000G_CCB_1_22_imputed_easyQC.txt", quote=FALSE, row.names = F, sep= "\t")


d = ACCORD_2

head(d)

d = read.table("~/data/accord/data/analysis/icaps.sent.DATA/icaps.ccb.all_races/association_cv/allChrImputed.ICAPS_po.assoc", header=TRUE, sep =" ")
colnames(d) = c("CHR","SNP","BP","A1","A2","FRQ","HWE_CHI2","OR","SE","P")

head(d)
dim(ACCORD_2)
head(ACCORD)
dim(d)
head(d)

tmp.d = d
rm(list=c('tmp.d'))
## Get only RSIDs from ACCORD_2

d = d[which(d$SNP %in% 'rs115033199')]
d = d[which(d$SNP %in% ACCORD$SNP),]
dim(d)
dim(ACCORD)

d$PZ = convert.z.score(d$Zstat,)

head(ACCORD)
tmp.d = ACCORD[which(ACCORD$SNP %in% d$SNP),] 


tmp.d = subset(tmp.d, select=c(SNP,CHR, BP, A1_allImp, A2_allImp, FRQ_allImp, OR_allImp, SE_allImp, P_allImp))


head(tmp.d)
head(d)
d = d[order(d$CHR,d$BP),]
tmp.d = tmp.d[order(tmp.d$CHR,tmp.d$BP),]
dim(tmp.d)
tmp.d = tmp.d[which(tmp.d$SNP %in% d$SNP),]
dim(tmp.d)
dim(d)
d = d[which(d$SNP %in% tmp.d$SNP),]
dim(d)
which(d$SNP %in% 'rs1650199')
d[5407466,]

tmp.d <- tmp.d[-is.na(tmp.d$A1_allImp),]
dim(tmp.d)


d = d[order(d$CHR,d$SNP),]
tmp.d = tmp.d[order(tmp.d$CHR,tmp.d$SNP),]


head(which(duplicated(tmp.d$SNP[-which(tmp.d$SNP=='.')])))

tmp.d[8488599,]
tmp.d[tmp.d$SNP=='AX-40870011',]
tmp.d[tmp.d$SNP==tmp.d$SNP[8488599],]
tmp.d[tmp.d$SNP==tmp.d$SNP[8525449],]

d$Zstat = d$Beta/d$SE

head(d)
plot(d$PZ,d$Pvalue)



convert.z.score<-function(z, one.sided=NULL) {
    if(is.null(one.sided)) {
        pval = pnorm(-abs(z));
        pval = 2 * pval
    } else if(one.sided=="-") {
        pval = pnorm(z);
    } else {
        pval = pnorm(-z);                                                                                 
    }
    return(pval);
}   

x = rnorm(100)
# normalize
min.x = min(x)
max.x = max(x)
t2 = (x - min.x)/(max.x - min.x)
print(t2)

sample(c(0,1),size=10,replace=T)
rnorm(n=10,0,1)
########### Random data
## Generate "response" 0s and 1s
t1 = sample(c(0,1),size=100,replace=T)
## Generate "dosage" between 0...1
x = rnorm(100)
# normalize
min.x = min(x)
max.x = max(x)
t2 = (x - min.x)/(max.x - min.x)
print(t2)
## Model
f = "t1=t2"
fit <- glm(t1~t2,family="binomial")

summary(fit)
## Lifted from compute_cv.r
coef<-summary(fit)$coef
beta<-as.numeric(format(coef["t2","Estimate"],digits=6))
se<-as.numeric(format(coef["t2","Std. Error"],digits=6))
p<-format(coef["t2","Pr(>|z|)"],digits=6)
n<- length(fit$fitted.values)

## Calc Z score
z = beta/se

## Get p value
convert.z.score<-function(z, one.sided=NULL) {
    if(is.null(one.sided)) {
        pval = pnorm(-abs(z));
        pval = 2 * pval
    } else if(one.sided=="-") {
        pval = pnorm(z);
    } else {
        pval = pnorm(-z);                                                                                 
    }
    return(pval);
}   
convert.z.score(z)





ACCORD <- read.table ("~/downloads/allChrom_ICAPS_po_gc_adj.csv", header=TRUE, sep =",")
d = read.table("~/data/accord/data/analysis/icaps.sent.DATA/icaps.ccb.all_races/association_cv/allChrImputed.ICAPS_po.assoc", header=TRUE, sep =" ")
colnames(d) = c("CHR","SNP","BP","A1","A2","FRQ","HWE_CHI2","OR","SE","P")
head(d)
head(ACCORD)
## Drop p value


d = d[which(d$SNP %in% ACCORD$SNP),]
tmp.d = ACCORD[which(ACCORD$SNP %in% d$SNP),] 

## grab one chromosome
d1 = d[d$CHR==1,]
dim(d1)

a1 = ACCORD[ACCORD$CHR==1,]

a1 = a1[-which(a1$SNP=='.'),]
a1 = a1[-which(is.na(a1$P_allImp)),]
dim(a1)
d1 = d1[which(d1$SNP %in% a1$SNP),]
dim(d1)

d1 = d1[order(d1$BP),]
head(d1)
a1 = a1[order(a1$BP),]
head(a1)

plot(d1$P,a1$P_allImp)

d1 = d1[,c("SNP",,"OR","SE","P")]
head(d1)
a1 = a1[,c("SNP","OR_allImp","SE_allImp","P_allImp")]
head(a1)
dd = merge(d1,a1,by="SNP")
head(dd)

plot(dd$P,dd$P_allImp)

dd$BetaA=log(dd$OR_allImp)
dd$SEA=log(dd$SE_allImp)
head(dd)

dd$Z=dd$BetaA/dd$SEA

dd$PZ=convert.z.score(dd$Z)

plot(dd$PZ,dd$P_allImp)


dim(d)
dim(ACCORD)

ACCORD = ACCORD[,c("SNP","OR_allImp","SE_allImp","P_allImp")]
d = d[,c("SNP",,"OR","SE","P")]

dd = ACCORD
dd$BetaA=log(dd$OR_allImp)
dd$SEA=log(dd$SE_allImp)
head(dd)

dd$Z=dd$BetaA/dd$SEA

dd$PZ=convert.z.score(dd$Z)

png('~/Desktop/caitrin.plot.png',height=1000,width=1000)
plot(dd$PZ,dd$P_allImp)
graphics.off()
