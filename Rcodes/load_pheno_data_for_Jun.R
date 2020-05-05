load_pheno_data<-function(path) {
  # run this from the accord analysis directory
  # this is the file with your risk alleles
  data<-read.table("RA file path"),header=T,stringsAsFactors=F)
  names(data)[1]="FID"
  names(data)[2]="IID"
  
  #pull phenotypes so they don't get centered/scaled #this is what is new from old version
  #phenos <- unlist(read.table(file.path(path,"phenotypes.txt"),stringsAsFactors=FALSE))

  # Scale continuous variables
#   data_class<-sapply(data,class)
#   data_levels<-sapply(data,function(x) length(levels(factor(x))))
#   scale_idx<-data_class=="numeric" | (data_class=="integer" & data_levels>8) & c(0,0,rep(1,dim(data)[2]-2))
#   scale_idx[which(names(data_class)%in%phenos)] <- FALSE
#   std<-function(x){(x-mean(x,na.rm=T))/sd(x,na.rm=T)}
#   data[,scale_idx]<-data.frame(sapply(data[,scale_idx],std))

#   # Define categorical variables
#   vars<-names(data)[data_class=="integer" & data_levels<=8]
#   if(length(grep("prt_hba1c_timegap",vars))>0) vars <- vars[which(vars!="prt_hba1c_timegap")]
#   for (i in 1:length(vars)) {data[,vars[i]]<-factor(data[,vars[i]])}

  # Create full list for all samples with imputed genotypes
  samples<-read.table(pipe(paste0("cut -d ' ' -f1,2 ../imputation/data/chr1.unphased.fam")),header=F,stringsAsFactors=F)
  samples$order<-seq(len=nrow(samples))
  names(samples)<-c("FID","IID","order")
  pheno_data<-merge(samples,data,sort=F,all.x=T)
  pheno_data<-pheno_data[sort.list(pheno_data$order), -3]

  if(any(grepl("AX",colnames(pheno_data)))){
      class(pheno_data[,grep("AX",colnames(pheno_data))]) <- "numeric"
  }
  return(pheno_data)
}
