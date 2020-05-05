rm(list=ls())
debug <- FALSE
#debug <- TRUE
if (debug){
    #path <- '~/data/accord/data/analysis/icaps/icaps.thiazide.all_races/'
    path <- '~/data/accord/data/analysis/mindcog/mindcog_cog_dsst'
    setwd('~/data/accord/data/analysis/')
    pheno = "dsst"
    chr <- 21
}else{
    args<-(commandArgs(TRUE))
    path<-args[1]
    pheno <- args[2]
    chr <- args[3]
}


dostuff = function(subpath){
    cat('Combining all chunks for chromosome ',chr,' and phenotype ',pheno,'\n')
    ## Prepare filenames (for chunks) and prepare output file (delete if exists)
    file.list <- list.files(file.path(path,subpath,"association_cv","imputed_chunks"))
    file.list <- file.list[grep(".assoc",file.list)]
    file.table <- t(data.frame(sapply(file.list,FUN = function(x) strsplit(x,split=".",fixed=TRUE))))
    chr.files <- file.table[which(file.table[,1]==paste0("chr",chr)),]
    chr.files.pheno <- chr.files[chr.files[,3]==pheno,]
    chr.files.pheno <- chr.files.pheno[order(as.numeric(chr.files.pheno[,2])),]
    merge.files <- rownames(chr.files.pheno)
    out.name <- paste(chr.files.pheno[1,c(1,3,4)],collapse=".")
    if(file.exists(file.path(path,subpath,"association_cv",out.name))){
        file.remove(file.path(path,subpath,"association_cv",out.name))
    }

    ## Build data frame for whole chromosome (combine chunks)
    out.file <- NULL
    f <- merge.files[2]
    for(f in merge.files){
        cat("Merging:",f," File #",which(merge.files==f)," of ",length(merge.files),"\n")
        temp.file <- try(read.table(file.path(path,subpath,"association_cv","imputed_chunks",f),header=TRUE,sep='\t'),silent = TRUE)

        if(class(temp.file)!="try-error"){
            out.file <-rbind(out.file,temp.file)
        }else{
    	cat('**ERROR** Potential Error, try-error, read error for',file.path(path,subpath,"association_cv","imputed_chunks",f),'\n')	
        }
    }
    if (debug){
        print('Printing...')
        print(out.name)
    }
    print(paste0("Writing results to ",out.name))
    write.table(out.file,file.path(path,subpath,"association_cv",out.name),row.names=FALSE)
    cat("Done\n")
}


dostuff("apoe")
dostuff("acr")