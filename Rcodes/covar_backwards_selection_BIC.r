rm(list=ls())

debug <- FALSE
#debug <- TRUE
if (debug){
    #path <- '~/data/accord/data/analysis/GV/GV_allraces'
    path <- '~/data/accord/data/analysis/GV/GV_arm1_allraces'
}else{
    args<-(commandArgs(TRUE))
    path<-args[[1]]
}



#cat("path:",path,"\n")

# Load phenotypes and forced covariate list
pheno<-unlist(read.table(paste0(path,"/phenotypes.txt"),stringsAsFactors=F))
model.types <- unlist(read.table(paste0(path,"/modeltypes.txt"),stringsAsFactors=F))
forced_covars<-unlist(read.table(paste0(path,"/forced_covars.txt"),stringsAsFactors=F))

#path <- 'Metformin/met_90d_allraces4'
#setwd('../')
#getwd()
# Load phenotype data
source("bin/load_pheno_data.r")
data<-load_pheno_data(path)
data<-data[complete.cases(data),]

num.unique <- apply(data,2,function(x)length(unique(x)))
## remove covars with <2 values
if(any(num.unique<2)){
  rm.vals <- names(num.unique[num.unique<2])
  forced_covars <- forced_covars[!forced_covars%in%rm.vals]
  data <- data[,!colnames(data)%in%rm.vals]
}

# Open model output file
model_out<-file(paste0(path,"/pheno_data/models.txt"),open="w")

#pheno_i <- 1
# Backwards select
#pheno_i <- 9
for (pheno_i in 1:length(pheno)) {
  
  model.type <- model.types[pheno_i]
  cat("Phenotype,",pheno[pheno_i]," will be analyzed with model type = ",model.type,"\n")
  out<-file(paste0(path,"/outputs/covar_backwards_selection_BIC_",pheno[pheno_i],".txt"),open="w")
  pheno_blr<-c()
  for (i in 1:length(pheno)) {
    if (grepl("_lr",pheno[i])) {
      pheno_blr<-c(pheno_blr,paste0("pre_",strsplit(pheno[i],'_')[[1]][1]))
    }
    if (grepl("prt_adj_hba1c",pheno[i])) {
      pheno_blr<-c(pheno_blr,"Pretreatment_hba1c")
    }
  }
  covars<-names(data)[!names(data)%in%c("FID","IID",pheno,forced_covars,pheno_blr)]
  pheno_blr<-c()
  if (grepl("_lr",pheno[pheno_i])) {
    pheno_blr<-c(pheno_blr,paste0("pre_",strsplit(pheno[pheno_i],'_')[[1]][1]))
  }
  if (grepl("prt_adj_hba1c",pheno[i])) {
    pheno_blr<-c(pheno_blr,"Pretreatment_hba1c")
  }
  bic_old<-100001
  bic_new<-100000
  while(bic_new<bic_old && length(covars)>0) {
    f<-as.formula(paste0(pheno[pheno_i],"~",paste(c(forced_covars,pheno_blr,covars),collapse="+")))
    if(model.type=="linear"){
      fit<-lm(f,data)
    }else if(model.type=="logistic"){
      fit <- glm(f,data=data,family="binomial")
    }
    bic_old<-BIC(fit)
    f_str<-paste(as.character(f)[2],as.character(f)[1],as.character(f)[3])
    cat("--------------------------------------------------\n",file=out)
    cat(paste0("BIC ",signif(bic_old,6),": ",f_str,"\n"),file=out)
    cat("--------------------------------------------------\n",file=out)
    for (i in 1:length(covars)) {
      test_covars<-covars[-i]
      f<-as.formula(paste0(pheno[pheno_i],"~",paste(c(forced_covars,pheno_blr,test_covars),collapse="+")))
      if(model.type=="linear"){
        fit<-lm(f,data)
      }else if(model.type=="logistic"){
        fit <- glm(f,data=data,family="binomial")
      }      
      bic_test<-BIC(fit)
      cat(paste0("  BIC=",signif(bic_test,6)," if remove ",covars[i],"\n"),file=out)
      if (bic_test<bic_new) {
        bic_new<-bic_test
        remove_covar<-i
      }
    }
    if (bic_new<bic_old) {
      covars<-covars[-remove_covar]
    }    
  }
  close(out)
  # Output final fit
  f<-as.formula(paste0(pheno[pheno_i],"~",paste(c(forced_covars,pheno_blr,covars),collapse="+")))
  if(model.type=="linear"){
    fit<-lm(f,data)
  }else if(model.type=="logistic"){
    fit <- glm(f,data=data,family="binomial")
  }
  sink(paste0(path,"/outputs/covar_backwards_selection_BIC_",pheno[pheno_i],".txt"),append=T)
  print(summary(fit))
  sink()
  # Write files for PLINK
  pheno_data<-data[,c("FID","IID",pheno[pheno_i])]
  covar_data<-data[,c("FID","IID")]
  covar_names<-c("FID","IID")
  temp_covar_data<-data[,c(forced_covars,pheno_blr,covars)]
  for (i in 1:dim(temp_covar_data)[2]) {
    if (is.factor(temp_covar_data[,i])) {
      lev<-levels(temp_covar_data[,i])
      temp<-data.frame(as.numeric(temp_covar_data[,i]==lev[1]))
      for (j in 2:length(lev)) {
        temp<-cbind(temp,as.numeric(temp_covar_data[,i]==lev[j]))
        covar_names<-c(covar_names,paste0(names(temp_covar_data)[i],lev[j]))
      }
      covar_data<-cbind(covar_data,temp[,-1])
    } else {
      covar_names<-c(covar_names,names(temp_covar_data)[i])
      covar_data<-cbind(covar_data,temp_covar_data[,i])
    }
  }
  names(covar_data)<-covar_names
  write.table(pheno_data,paste0(path,"/pheno_data/pheno_",pheno[pheno_i],".txt"),col.names=T,row.names=F,quote=F,sep="\t")
  write.table(covar_data,paste0(path,"/pheno_data/covar_",pheno[pheno_i],".txt"),col.names=T,row.names=F,quote=F,sep="\t")
  # Write files for GCTA
  covar_data<-data[,c("FID","IID",forced_covars,pheno_blr,covars)]
  data_class<-sapply(covar_data,class)
  scale_idx<-data_class%in%c("integer","factor") | c(1,1,rep(0,dim(covar_data)[2]-2))
  d_data<-covar_data[,scale_idx]
  scale_idx<-!data_class%in%c("integer","factor") | c(1,1,rep(0,dim(covar_data)[2]-2))
  q_data<-covar_data[,scale_idx]
  write.table(d_data,paste0(path,"/gcta/dcovar_",pheno[pheno_i],".txt"),col.names=F,row.names=F,quote=F,sep="\t")
  write.table(q_data,paste0(path,"/gcta/qcovar_",pheno[pheno_i],".txt"),col.names=F,row.names=F,quote=F,sep="\t")
  write.table(pheno_data,paste0(path,"/gcta/pheno_",pheno[pheno_i],".txt"),col.names=F,row.names=F,quote=F,sep="\t")
  # Write model
  cat(paste0(pheno[pheno_i],"~",paste(c(forced_covars,pheno_blr,covars),collapse="+"),"\n"),file=model_out)
}
close(model_out)

