## compare list of snps from T2D clustering - Li and Dudley et al. 
library(data.table)
##read dudley snps by cluster
dclustA<- read.csv("~/NCSU_Projects/Accord/Data/analysis/Clustering/dudley_et_al/sig_snps_cluster_A.csv")
dclustB<- read.csv("~/NCSU_Projects/Accord/Data/analysis/Clustering/dudley_et_al/sig_snps_cluster_B.csv")
dclustC<- read.csv("~/NCSU_Projects/Accord/Data/analysis/Clustering/dudley_et_al/sig_snps_cluster_C.csv")

d.all.snps <- unique(c(as.character(dclustA$dbSNP),as.character(dclustB$dbSNP),as.character(dclustC$dbSNP)))


##read in SNPs from ACCORD
aclust1 <- fread("~/NCSU_Projects/Accord/accord_data/analysis/clinical_cluster_gwas/association_cv/all_commonVariantAnalysis_clust_1_results.tsv")
a.snps.sig1 <- aclust1[which(aclust1$P_unadj<1E-5),]
aclust1 <- aclust1[which(as.character(aclust1$SNP) %in% c(d.all.snps)),]

aclust2 <- fread("~/NCSU_Projects/Accord/accord_data/analysis/clinical_cluster_gwas/association_cv/all_commonVariantAnalysis_clust_2_results.tsv")
a.snps.sig2 <- aclust2[which(aclust2$P_unadj<1E-5),]
aclust2 <- aclust2[which(as.character(aclust2$SNP) %in% c(d.all.snps)),]

aclust3 <- fread("~/NCSU_Projects/Accord/accord_data/analysis/clinical_cluster_gwas/association_cv/all_commonVariantAnalysis_clust_3_results.tsv")
a.snps.sig3 <- aclust3[which(aclust3$P_unadj<1E-5),]
aclust3 <- aclust3[which(as.character(aclust3$SNP) %in% c(d.all.snps)),]


hist(aclust1$P_unadj)
qqplot(aclust1$P_unadj,runif(length(aclust1$P_unadj)))
abline(a=0,b=1)

hist(aclust2$P_unadj)
qqplot(aclust2$P_unadj,runif(length(aclust2$P_unadj)))
abline(a=0,b=1)

hist(aclust3$P_unadj)
qqplot(aclust3$P_unadj,runif(length(aclust3$P_unadj)))
abline(a=0,b=1)

## test dudley cluster A against each Accord cluster
pdf("~/NCSU_Projects/Accord/accord_data/analysis/clinical_cluster_gwas/compare_Dudley/dudley_vs_accord_snps_ks_test.pdf",
    width=12,height =12)
layout(matrix(c(1:9),ncol=3,byrow=TRUE))

DaA1 <- aclust1[aclust1$SNP %in% as.character(dclustA$dbSNP),]
DaA1.ks <- ks.test(DaA1$P_unadj,"punif",alternative="two.sided")
hist(DaA1$P_unadj,main="Li et al_A vs Accord_A")
legend("topright",paste0("ks p: ",signif(DaA1.ks$p.value,3)))

DaA2 <- aclust2[aclust2$SNP %in% as.character(dclustA$dbSNP),]
DaA2.ks <- ks.test(DaA2$P_unadj,"punif",alternative="two.sided")
hist(DaA2$P_unadj,main="Li et al_A vs Accord_B")
legend("topright",paste0("ks p: ",signif(DaA2.ks$p.value,3)))


DaA3 <- aclust3[aclust3$SNP %in% as.character(dclustA$dbSNP),]
DaA3.ks <- ks.test(DaA3$P_unadj,"punif",alternative="two.sided")
hist(DaA3$P_unadj,main="Li et al_A vs Accord_C")
legend("topright",paste0("ks p: ",signif(DaA3.ks$p.value,3)))


## test dudley cluster B against each Accord cluster

DbA1 <- aclust1[aclust1$SNP %in% as.character(dclustB$dbSNP),]
DbA1.ks <- ks.test(DbA1$P_unadj,"punif",alternative="two.sided")
hist(DbA1$P_unadj,main="Li et al_B vs Accord_A")
legend("topright",paste0("ks p: ",signif(DbA1.ks$p.value,3)))


DbA2 <- aclust2[aclust2$SNP %in% as.character(dclustB$dbSNP),]
DbA2.ks <- ks.test(DbA2$P_unadj,"punif",alternative="two.sided")
hist(DbA2$P_unadj,main="Li et al_B vs Accord_B")
legend("topright",paste0("ks p: ",signif(DbA2.ks$p.value,3)))

DbA3 <- aclust3[aclust3$SNP %in% as.character(dclustB$dbSNP),]
DbA3.ks <- ks.test(DbA3$P_unadj,"punif",alternative="two.sided")
hist(DbA3$P_unadj,main="Li et al_B vs Accord_C")
legend("topright",paste0("ks p: ",signif(DbA3.ks$p.value,3)))


## test dudley cluster C against each Accord cluster

DcA1 <- aclust1[aclust1$SNP %in% as.character(dclustC$dbSNP),]
DcA1.ks <- ks.test(DcA1$P_unadj,"punif",alternative="two.sided")
hist(DcA1$P_unadj,main="Li et al_C vs Accord_A")
legend("topright",paste0("ks p: ",signif(DcA1.ks$p.value,3)))


DcA2 <- aclust2[aclust2$SNP %in% as.character(dclustC$dbSNP),]
DcA2.ks <- ks.test(DcA2$P_unadj,"punif",alternative="two.sided")
hist(DcA2$P_unadj,main="Li et al_C vs Accord_B")
legend("topright",paste0("ks p: ",signif(DcA2.ks$p.value,3)))

DcA3 <- aclust3[aclust3$SNP %in% as.character(dclustC$dbSNP),]
DcA3.ks <- ks.test(DcA3$P_unadj,"punif",alternative="two.sided")
hist(DcA3$P_unadj,main="Li et al_C vs Accord_C")
legend("topright",paste0("ks p: ",signif(DcA3.ks$p.value,3)))
graphics.off()

ks.results <- rbind(c("Dudley_A","Accord_A",DaA1.ks$p.value),
                    c("Dudley_A","Accord_B",DaA2.ks$p.value),
                    c("Dudley_A","Accord_C",DaA3.ks$p.value),
                    c("Dudley_B","Accord_A",DbA1.ks$p.value),
                    c("Dudley_B","Accord_B",DbA2.ks$p.value),
                    c("Dudley_B","Accord_C",DbA3.ks$p.value),
                    c("Dudley_C","Accord_A",DcA1.ks$p.value),
                    c("Dudley_C","Accord_B",DcA2.ks$p.value),
                    c("Dudley_C","Accord_C",DcA3.ks$p.value))
                    

#write csv
write.csv(ks.results,"~/NCSU_Projects/Accord/accord_data/analysis/clinical_cluster_gwas/compare_Dudley/dudley_vs_accord_snps_ks_test.csv",row.names=FALSE)


### look for snp enrichment with accord snps <1e-5

all_acc <- rbind(a.snps.sig1,a.snps.sig2,a.snps.sig3)
sig.snps <- unique(test$SNP)
all_dud <- rbind(dclustA,dclustB,dclustC)
sig.dud <- as.character(unique(all_dud$dbSNP))
intersect(sig.snps,sig.dud)

# No overlap so code below not needed:

# a1_da <- intersect(a.snps.sig1$SNP,as.character(dclustA$dbSNP))
# a1_db <- intersect(a.snps.sig1$SNP,as.character(dclustB$dbSNP))
# a1_dc <- intersect(a.snps.sig1$SNP,as.character(dclustC$dbSNP))
# 
# a2_da <- intersect(a.snps.sig2$SNP,as.character(dclustA$dbSNP))
# a2_db <- intersect(a.snps.sig2$SNP,as.character(dclustB$dbSNP))
# a2_dc <- intersect(a.snps.sig2$SNP,as.character(dclustC$dbSNP))
# 
# a3_da <- intersect(a.snps.sig3$SNP,as.character(dclustA$dbSNP))
# a3_db <- intersect(a.snps.sig3$SNP,as.character(dclustB$dbSNP))
# a3_dc <- intersect(a.snps.sig3$SNP,as.character(dclustC$dbSNP))






