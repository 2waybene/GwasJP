rm(list=ls())
debug <- FALSE
#debug <- TRUE
if (debug){
    #path <- '~/data/accord/data/analysis/icaps/icaps.aceia.arb.all_races/'
    #getwd()
    #setwd('~/data/accord/data/analysis/')
    path <- '~/data/accord/data/analysis/eye_disease/eye_dis_allraces/'
    #pheno <- ''
    setwd('~/data/accord/data/analysis/')
    ##setwd("~/NCSU_Projects/Accord/Data/analysis")
}else{
    args<-(commandArgs(TRUE))
    path<-args[1]
    pheno <- args[2]
}

library(data.table)

file.list <- list.files(file.path(path,"association_cv"))
file.list <- file.list[grep(".assoc",file.list)]
file.list <- file.list[grep("chr",file.list)]
file.list <- file.list[grep("chr0",file.list,invert=TRUE)]
chr.files <- t(data.frame(sapply(file.list,FUN = function(x) strsplit(x,split=".",fixed=TRUE))))
#chr.files <- file.table[which(file.table[,1]==paste0("chr",chr)),]

#pheno.list <- unique(chr.files[,2])
#pheno <- pheno.list[1]
#for(pheno in pheno.list){

## Get file list (per chromosome)
chr.files.pheno <- chr.files[chr.files[,2]==pheno,]                                        
merge.files <- rownames(chr.files.pheno)
out.name <- paste0("allChrImputed.",pheno,".assoc")
#if(file.exists(file.path(path,"association_cv",out.name))){
#    file.remove(file.path(path,"association_cv",out.name))
#}

cat('** Merge results for phenotype: ',pheno,' **\n')

## Loop through chr files and output to single file
for(f in merge.files){
    cat("Merging:",f,"-- File #",which(merge.files==f)," of ",length(merge.files),"\n")
    temp.file <- try(fread(file.path(path,"association_cv",f),header=TRUE),silent = TRUE)
    if(any(class(temp.file)!="try-error")){
        if(which(merge.files==f)==1){ ## This overwrites the current file
            write.table(temp.file,file.path(path,"association_cv",out.name),row.names=FALSE)
        }else{
            write.table(temp.file,file.path(path,"association_cv",out.name),row.names=FALSE,col.names=FALSE,append=TRUE)
        }              
    }
}
cat("Done\n")
