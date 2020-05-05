rm(list=ls())

debug <- FALSE
#debug <- TRUE
if (debug){
    #path <- '~/data/accord/data/analysis/GV/GV_allraces'
    path <- '~/data/accord/data/analysis/GV/GV_white'
}else{
    args<-(commandArgs(TRUE))
    path<-args[[1]]
}

## If modeltypes exists, erase it so we can create+append in the for loop below.
if (file.exists(paste0(path,'/modeltypes.txt'))){
    unlink(paste0(path,'/modeltypes.txt'))
}

## Load initial phenotype data
data <- read.table(paste0(path,"/pheno_data/pheno_data_step1.txt"),header=T)
pheno <- read.table(paste0(path,"/phenotypes.txt"),stringsAsFactors=F)

#head(data)
#head(pheno)

if (debug){
    ## Test logistic fake pheno
    #data$hellothere <- c(rep(0,1/2*nrow(data)),rep(1,round(1/2*nrow(data))))
    #pheno <- rbind(pheno  ,'hellothere')
    pheno_i <- 1
}

for (pheno_i in 1:nrow(pheno)) {
    ## Check number of unique phenotype values in the pheno_data file
    ## If 2 values, then logistic (glm) model is used, else linear model is used.
    phen.vals <- length(unique(na.omit(data[,pheno[pheno_i,]])))

    if (!file.exists(paste0(path,'/modeltypes.txt'))){
        if (phen.vals>2){
            write(file=paste0(path,'/modeltypes.txt'),'linear')
        }else{
            write(file=paste0(path,'/modeltypes.txt'),'logistic')
        }    
    }else{
        if (phen.vals>2){
            write(file=paste0(path,'/modeltypes.txt'),'linear',append=TRUE)
        }else{
            write(file=paste0(path,'/modeltypes.txt'),'logistic',append=TRUE)
        }    
    }
}

