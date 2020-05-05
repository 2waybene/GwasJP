rm(list=ls())

#library(NCBI2R) #install.packages('https://cran.r-project.org/src/contrib/Archive/NCBI2R/NCBI2R_1.4.7.tar.gz',repos=NULL,type='source')
library(data.table)
library(biomaRt)
#library(rsnps) #install.packages("rsnps")


debug <- FALSE
#debug <- TRUE
if (debug){
    path <- '~/data/accord/data/analysis/Metformin/met_90d_allraces4'
    ##path <- 'Metformin/met_90d_allraces4'
    ##setwd('~/NCSU_Projects/Accord/Data/analysis')
    adjusted <- FALSE
}else{
    args<-commandArgs(T)
    path<-args[1]
    adjusted <- args[2]
}

p_value_cutoff<-1e-5

phenotypes<-unlist(read.table(paste0(path,"/phenotypes.txt"),stringsAsFactors=F))
model.types <- unlist(read.table(paste0(path,"/modeltypes.txt"),stringsAsFactors=F))
#pheno_i <- 1
for (pheno_i in 1:length(phenotypes)) {
  file.name <- paste0(phenotypes[pheno_i],"_gc_adj.assoc.linear")
  peak.vals <- fread(file.path(path,"peak_data",file.name))
  
  if(adjusted){
    peak.vals <- peak.vals[which(peak.vals[,"P_adj"]<p_value_cutoff),]
    #loc<-data$BP[data$P_adj<=p_value_cutoff]
  }else{
    peak.vals <- peak.vals[which(peak.vals[,"P_unadj"]<p_value_cutoff),]
    #loc<-data$BP[data$P_unadj<=p_value_cutoff]        
  }
  
  
  
  
  #dim(peak.vals)
  #head(peak.vals)
    
  snplist<-peak.vals[,"SNP"]
  #head(snplist)
  #snplist[1:20]
  #dim(snplist)
  #snplist <- snplist[which(snplist!=".")]
  #snplist <- snplist[grep("AX",snplist[,1],invert=TRUE)]
  #snplist <- snplist[grep("es",snplist[,1],invert=TRUE)]
  #snplist <- snplist[grep("ex",snplist[,1],invert=TRUE)]
  #snplist <- snplist[grep("kg",snplist[,1],invert=TRUE)]
  snplist <- snplist[like(SNP,"rs"),]
  snplist <- snplist[,SNP]
 
  
  
  
  snpmart = useMart(biomart = "ENSEMBL_MART_SNP", dataset = "hsapiens_snp")
  
  # Check which filters and attributes we wan't to use:
  # listAttributes(snpmart)
  # listFilters(snpmart)
  
  # result
  output <- getBM(attributes = c("refsnp_id", "chr_name", "chrom_start", "chrom_end", "ensembl_gene_stable_id"), 
        filters = c("snp_filter"), 
        values = snplist, 
        mart = snpmart)
  
  #get gene names
  genemart = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
  
  
  gene.output <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "external_gene_name"), 
                  filters = c("ensembl_gene_id"), 
                  values = output[,"ensembl_gene_stable_id"], 
                  mart = genemart)
  
  
  merge.out <- merge(output,gene.output, by.x="ensembl_gene_stable_id",by.y="ensembl_gene_id",all.x=TRUE)
  
   #output <- GetSNPInfo(snplist[,SNP])
  #query cannot handle large lists
  # output <- NULL
  # if(length(snplist)>300){
  #   start <- 1
  #   while(start<length(snplist)){
  #     
  #     if((start+299)>length(snplist)){
  #       end <- length(snplist)
  #       #print("here1")
  #     }else{
  #       end <- start+299
  #       #print("here2")
  #     }
  #     tmp.snplist <- snplist[start:end]
  #     #output <- rbind(output,NCBI_snp_query(tmp.snplist))
  #     output <- rbind(output,ncbi_snp_query(tmp.snplist))
  #     start <- end+1
  #     #print(start)
  #   }
  # }else{
  #   #output <- NCBI_snp_query(snplist)
  #   output <- ncbi_snp_query(snplist)
  # }
  # 
  # 
  colnames(merge.out) <- paste0("ENSEMBL_",colnames(merge.out))
  colnames(merge.out)[2] <- "SNP"
  dim(merge.out)
  dim(peak.vals)
  
  peak.vals[,"BETA"] <- 999
  if(model.types[pheno_i]=="logistic"){
    peak.vals[Type=='IMPU',"BETA"] <- peak.vals[Type=='IMPU',"OR_allImp"]
    peak.vals[Type=='GENO',"BETA"] <- peak.vals[Type=='GENO',"OR_geno"]
  }else{
    ## Collapse on beta
    peak.vals[Type=='IMPU',"BETA"] <- peak.vals[Type=='IMPU',"BETA_allImp"]
    peak.vals[Type=='GENO',"BETA"] <- peak.vals[Type=='GENO',"BETA_geno"]
  }
  peak.vals[Type=='META',"BETA"] <- peak.vals[Type=='META',"OR_meta"]
  peak.vals[Type=='IMPU',"A1"] <- peak.vals[Type=='IMPU',"A1_allImp"]
  peak.vals[Type=='IMPU',"A2"] <- peak.vals[Type=='IMPU',"A2_allImp"]
  
  
  matches <- match(merge.out[,1],unlist(peak.vals[,SNP]))
  length(matches)
  dim(peak.vals)

  #tail(output)
  #tail(peak.vals)
  #output[30:35,]
  #peak.vals[30:35,]
  join.data <- peak.vals[,c("SNP","BETA","Type","A1","A2","MAF","P_unadj","P_adj")]
  # beta <- peak.vals[,"BETA"]
  # type <- peak.vals[,"Type"]#peak.vals[matches,as.numeric(Type=='GENO'|Type=='META')]
  # A1 <- peak.vals[,"A1"]
  # A2 <- peak.vals[,"A2"]
  # maf <- peak.vals[,"MAF"]
  # p.vals <- peak.vals[,"P_unadj"]
  # p.vals.adj <- peak.vals[,"P_adj"]
  # output <- cbind(output,A1,A2,maf,beta,type,p.vals,p.vals.adj)
  final.out <- merge(join.data,merge.out,by="SNP",all.x=TRUE)
  dim(final.out)
  output.name <- paste0(phenotypes[pheno_i],"_snp_gene_mapping.csv")
  write.csv(final.out,file.path(path,"peak_data",output.name),row.names=FALSE)
  
}
