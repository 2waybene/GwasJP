args<-commandArgs(T)
path<-args[1]
chr<-args[2]
maf_cutoff<-as.numeric(args[3])

#cat("path:",path,"\n")
#cat("chr:",chr,"\n")
#cat("maf:",maf_cutoff,"\n")

#library(SKAT,lib="~/R/library")
library(SKAT)
########################################
# Which tests to perform
########################################
test_ind   <-1 # burden      - indicator
test_prop  <-1 # burden      - proportion
test_weight<-1 # burden      - weighted
test_skat  <-1 # non-burden  - SKAT
test_skato <-1 # combination - SKAT-O = SKAT + weighted burden

########################################
# Options
########################################
missingness_cutoff<-0.15 # Variants with missingness above cutoff are not included

########################################
# Generic burden test function
#   Inputs:
#     f - model function
#     X - phenotype and covaraite data
#     z - burden score
#   Outputs:
#     beta, se and p-value of z
########################################
burden<-function(f,X,z) {
  fit<-lm(f,data.frame(X,z))
  coef<-summary(fit)$coef
  if ("z"%in%row.names(coef)) {
    beta<-format(coef["z","Estimate"],digits=6)
    se<-format(coef["z","Std. Error"],digits=6)
    p<-format(coef["z","Pr(>|t|)"],digits=6)
  } else {beta<-NA;se<-NA;p<-NA}
  return(c(beta,se,p))
}

########################################
# Load genotype data (in SKAT SSD format)
########################################
File.SSD<-paste0("../imputation/merge.p13.functional/chr",chr,".SSD")
File.Info<-paste0("../imputation/merge.p13.functional/chr",chr,".SSD.info")
SSD.INFO<-Open_SSD(File.SSD, File.Info)

########################################
# Load models and phenotype data
########################################
models<-unlist(read.table(paste0(path,"/pheno_data/models.txt"),stringsAsFactors=F))
phenos<-sapply(strsplit(models,"~"),"[",1)
source("bin/load_pheno_data.r")
pheno_data<-load_pheno_data(path)
#make sure phenotypes are not factors
for(i in phenos){
  pheno_data[,i]<-as.numeric(pheno_data[,i])
}
attach(pheno_data)

########################################
# Only look at samples with complete data
########################################
sample_idx<-complete.cases(pheno_data)

########################################
# Initialize phenotype-specific models/files
########################################
null_model<-list()
f<-list()
file_out<-list()
for (i in 1:length(phenos)) {
  # Fit SKAT null model
  suppressWarnings(null_model[[i]]<-SKAT_Null_Model(as.formula(models[i]),out_type="C"))
  
  # Formula for burden tests
  f[[i]]<-as.formula(paste0(models[i],"+z"))
  
  # Start output file
  file_out[[i]]<-file(paste0(path,"/association_rv/chr",chr,"_",phenos[i],".txt"),open="w")
  cat("geneID",file=file_out[[i]])
  if (test_ind)    {cat(paste0("\tind_",   c("beta","se","p")),sep="",file=file_out[[i]])}
  if (test_prop)   {cat(paste0("\tprop_",  c("beta","se","p")),sep="",file=file_out[[i]])}
  if (test_weight) {cat(paste0("\tweight_",c("beta","se","p")),sep="",file=file_out[[i]])}
  if (test_skat)   {cat(paste0("\tskat_",  c("stat","p")),sep="",file=file_out[[i]])}
  if (test_skato)  {cat(paste0("\tskato_", c("stat","p")),sep="",file=file_out[[i]])}
  cat("\n",file=file_out[[i]])
}

########################################
# Loop through genes
########################################
for (gene_i in 1:SSD.INFO$nSets) {
  #cat("gene_i",gene_i,"\n")
#  gene_name<-sapply(strsplit(SSD.INFO$SetInfo$SetID[gene_i],":"),"[",2)
  gene_name<-SSD.INFO$SetInfo$SetID[gene_i]
  
  # Dont include large loci
#  if (!grepl("^TRA$|^IGH$|@$",gene_name)) {
  if (SSD.INFO$SetInfo$SetSize[gene_i]<1500) {
    geno_data<-Get_Genotypes_SSD(SSD.INFO,gene_i)
    geno_data[geno_data==9]<-NA
    
    # Need at least 2 variants to proceed
    if (dim(geno_data)[2]>1) {
      # Compute MAF
      maf<-apply(geno_data[sample_idx,],2,function(x) sum(x,na.rm=T)/2/sum(!is.na(x)))
      
      # Compute variant missingness
      missingness<-apply(geno_data[sample_idx,],2,function(x) sum(is.na(x))/dim(geno_data)[1])
      
      # Only use variants with missingness <= cutoff and maf <= maf_cutoff
      variant_idx<- missingness<=missingness_cutoff & maf<=maf_cutoff
      variant_n<-sum(variant_idx)
      
      # Need at least 2 variants passing cutoffs to proceed
      if (variant_n>1) {
        for (i in 1:length(phenos)) {
          cat(gene_name,file=file_out[[i]])
        }
        
        ########################################
        # Burden tests
        ########################################
        # Indicator of any rare allele
        if (test_ind) {
          z<-apply(geno_data[,variant_idx],1,function(x) 1*any(x>0,na.rm=T))
          for (i in 1:length(phenos)) {
            test<-burden(f[[i]],pheno_data,z)
            cat(paste0("\t",test),sep="",file=file_out[[i]])
          }
        }
        
        # Proportion of rare alleles
        if (test_prop) {
          z<-apply(geno_data[,variant_idx],1,function(x) sum(x,na.rm=T)/2/sum(!is.na(x)))
          for (i in 1:length(phenos)) {
            test<-burden(f[[i]],pheno_data,z)
            cat(paste0("\t",test),sep="",file=file_out[[i]])
          }
        }
        
        # Weighted sum
        if (test_weight) {
          w<-(maf*(1-maf))^-.5
          w[maf==0 | maf>maf_cutoff]<-0
          z<-apply(geno_data[,variant_idx],1,function(x) sum(x*w[variant_idx],na.rm=T))
          for (i in 1:length(phenos)) {
            test<-burden(f[[i]],pheno_data,z)
            cat(paste0("\t",test),sep="",file=file_out[[i]])
          }
        }
        
        ########################################
        # Non-burden tests
        ########################################
        w<-dbeta(maf,1,25)
        w[maf==0 | maf>maf_cutoff]<-0
        
        # SKAT - Q value is statistic
        if (test_skat) {
          for (i in 1:length(phenos)) {
            suppressWarnings(try_obj<-try(test<-SKAT(geno_data,null_model[[i]],weights=w,estimate_MAF=2),silent=T))
            if (is(try_obj, "try-error")) {
              cat("\tNA\tNA",file=file_out[[i]])
            } else {
              cat(paste0("\t",format(c(test$Q,test$p.value),digits=6,trim=T)),sep="",file=file_out[[i]])
            }
          }
        }
        
        # SKAT-O - min p value is statistic
        if (test_skato) {
          #print("testing skato")
          #print(dim(geno_data))
          #print(head(geno_data))
          for (i in 1:length(phenos)) {
            #print(i)
            suppressWarnings(try_obj<-try(test<-SKAT(geno_data,null_model[[i]],weights=w,estimate_MAF=2,method="optimal.adj"),silent=T))
            #try_obj<-try(test<-SKAT(geno_data,null_model[[i]],weights=w,estimate_MAF=2,method="optimal"),silent=F)
            #print(test)
            if (is(try_obj, "try-error")) {
              cat("\tNA\tNA",file=file_out[[i]])
              #print(test$Test.Type)
            } else if (test$param$n.marker.test <=1){#(is.na(test$Test.Type) | test$Test.Type!="optimal.adj") { #from SKAT <v1.3.2.1, no longer works in recent versions because test$Test.Type no longer exists
              # There may be a chance the missingness from SKAT is different than the one I computed
              # This can result in a single variant remaining and optimal.adj will revert back to davies method
              # If that occurs just print NA
              cat("\tNA\tNA",file=file_out[[i]])
            } else {
              cat(paste0("\t",format(c(NA,test$p.value),digits=6,trim=T)),sep="",file=file_out[[i]])
            }
          }
        }
        
        for (i in 1:length(phenos)) {
          cat("\n",file=file_out[[i]])
        }
      }
    }
  }
}

for (i in 1:length(phenos)) {
  close(file_out[[i]])
}
Close_SSD()

