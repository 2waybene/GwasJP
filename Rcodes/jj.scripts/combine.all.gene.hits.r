setwd('~/data/accord/data/analysis/')   

tmp <- read.csv(file='hgi_baseline/hgi_allraces/peak_data/hgi_snp_gene_mapping.csv',stringsAsFactors=FALSE)
d <- data.frame(GRP=rep('AllRaces',nrow(tmp)),tmp,stringsAsFactors=FALSE)

tmp <- read.csv(file='hgi_baseline/hgi_hispanic/peak_data/hgi_snp_gene_mapping.csv',stringsAsFactors=FALSE)
d <- data.frame(GRP=rep('AllRaces',nrow(tmp)),tmp,stringsAsFactors=FALSE)

d[1,]
