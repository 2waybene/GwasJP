setwd('~/data/accord/data/analysis')

paths <-c('icaps/icaps.aceia.arb.all_races/association_cv/imputed_chunks','icaps/icaps.bb.all_races/association_cv/imputed_chunks','icaps/icaps.ccb.all_races/association_cv/imputed_chunks','icaps/icaps.thiazide.all_races/association_cv/imputed_chunks')
#p <- paths[4]


for (p in paths){
    print(p)
    filez <- c(list.files(p,pattern='chr22'),list.files(p,pattern='chr12'))
    #f <- filez[1]
    for (f in filez){
        print(f)
        tmp <- read.table(file.path(p,f),sep='\t',header=T,fill=NA)
        tmp$N <- NULL
        head(tmp)
        colnames(tmp)[8] <- 'BETA'
        write.table(file=file.path(p,f),tmp,sep='\t',row.names=FALSE,quote=F)
    }
}


