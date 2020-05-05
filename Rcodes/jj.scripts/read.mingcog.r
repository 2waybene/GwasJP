d = read.table('~/data/accord/data/pheno_data/mindcog/accord_analysis_file_13March2018.txt',header=T,fill=NA,sep='\t')

snps = read.table('~/data/accord/data/pheno_data/mindcog/AADHS_replication_snps_11302017.csv',header=T,fill=NA)

dat = read.table('~/data/accord/data/analysis/pheno_data',header=T,fill=NA,sep='\t')


head(snps)
head(d)

table(d$techschool)



