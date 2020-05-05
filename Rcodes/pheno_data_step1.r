debug <- FALSE
#debug <- TRUE
if (debug){
    path <- '~/data/accord/data/analysis/icaps/icaps.thiazide.all_races/'
    setwd('~/data/accord/data/analysis/')
    adjusted <- FALSE
    
}else{
    args<-(commandArgs(TRUE))
    path<-args[1]
    fname <- args[2]
}


# Load full phenotype data file
data<-read.table(fname,header=T,stringsAsFactors=F,sep='\t')

#head(data)
#dim(data)

# Retain selected samples
if (file.exists(paste0(path,"/starting_samples.txt"))) {
    keep<-read.table(paste0(path,"/starting_samples.txt"),stringsAsFactors=F)
    data<-data[data$MaskID%in%keep[,1] & data$LabID%in%keep[,2],]
}else{
	print("WARNING: NO starting_samples.txt file was provided!")
}

# Retain selected variables, force BLR values into model if looking at 4-month difference
keep<-unlist(read.table(paste0(path,"/starting_covars.txt"),stringsAsFactors=F))
pheno<-unlist(read.table(paste0(path,"/phenotypes.txt"),stringsAsFactors=F))
pheno_blr<-c()
for (i in 1:length(pheno)) {
    if (grepl("_lr",pheno[i])) {
        pheno_blr<-c(pheno_blr,paste0("pre_",strsplit(pheno[i],'_')[[1]][1]))
    }
    if (grepl("prt_adj_hba1c",pheno[i])) {
      pheno_blr<-c(pheno_blr,"Pretreatment_hba1c")
    }
}

if (length(setdiff(intersect(pheno,names(data)),pheno))>0){
	cat("WARNING: Some of your phenos in phenotypes.txt are missing in your data (",fname,")\n")       
}

data<-data[,names(data)%in%c(keep,pheno,pheno_blr,"MaskID","LabID")]


# Hemoglobin glycation index (HGI)
if ("hgi"%in%c(keep,pheno)) {
    # Compute HGI
    fit<-lm(data$hba1c ~ data$fpg, na.action=na.exclude)
    hgi<-round(residuals(fit),digits=4)
    data<-cbind(data,hgi)
    # Split into 3 groups
    grp<-rep(NA,dim(data)[1])
    cutoffs<-quantile(residuals(fit),c(1/3,2/3),na.rm=T)
    grp[residuals(fit)<cutoffs[1]]<-1
    grp[residuals(fit)>=cutoffs[1] & residuals(fit)<cutoffs[2]]<-2
    grp[residuals(fit)>=cutoffs[2]]<-3
    grp<-factor(grp)
    # Plot HGI
    pdf(file=paste0(path,"/outputs/HbA1c_vs_FPG.pdf"),width=10, height=7.5)
    plot(data$fpg,data$hba1c,col=grp,pch=16,cex=0.5,xlab="FPG",ylab="HbA1c",main="HbA1c vs FPG")
    abline(a=coef(fit)[1],b=coef(fit)[2],col=4,lwd=2)
    eq<-paste("y = ",format(coef(fit)[1], digits = 4)," + ",format(coef(fit)[2], digits = 4)," * x",sep="")
    text(max(data$fpg,na.rm=T),min(data$hba1c,na.rm=T),label=eq,col=4,adj=1)
    trash<-dev.off()
}


## incorporate scores
score.covars <- c(keep,pheno)[grep("score",c(keep,pheno))]
if(length(score.covars)>0){
  for(cov in unique(score.covars)){
    if(grepl("g2avetid_time_adjusted_score",cov))next
    score.dir <- unlist(strsplit(cov,split="/"))[1]
    score.name <- unlist(strsplit(cov,split="/"))[2]
    score.file <- paste0(cov,".txt") 
    drug.class <- unlist(strsplit(score.name,"_"))[1]
    score.data <- read.table(file.path("..","pheno_data",paste0(cov,".txt")),header=TRUE,sep="\t")
    score0.data <- read.table(file.path("..","pheno_data",score.dir,paste0(drug.class,"_0.txt")),header=TRUE,sep="\t")
    score.keep <- c(score0.data[which(score0.data[,2]==1),"MaskID"],score.data[which(score.data[,2]==1),"MaskID"])
    score.sub <- score.data[which(score.data[,"MaskID"]%in%score.keep),] 
    data <- cbind(data,score.sub[match(data[,"MaskID"],score.sub[,"MaskID"]),2])
    colnames(data)[ncol(data)] <- score.name
  }
}


# Average Serum Creatinine [time frame matched to metformin score]
if ("avgCreatClr"%in%c(keep,pheno)) {
  # read in data file
  creat.data <- read.table(file.path("../pheno_data/avgCreatClr.txt"),header=TRUE,sep="\t")
  sub.creat <- creat.data[na.omit(match(data[,"MaskID"],creat.data[,"MaskID"])),]
  data <- cbind(data,sub.creat[,"avgCreatClr"])
  colnames(data)[ncol(data)] <- "avgCreatClr"
}

if (length(grep("g2avetid_time_adjusted_score",c(keep,pheno)))>0) {
  cov <- c(keep,pheno)[grep("g2avetid_time_adjusted_score",c(keep,pheno))]
  ins.data <- read.table(file.path("..","pheno_data",paste0(cov,".txt")),header=TRUE,sep="\t")
  sub.ins <- ins.data[match(data[,"MaskID"],ins.data[,"MaskID"]),]
  data <- cbind(data,sub.ins[,"Score"])
  colnames(data)[ncol(data)] <- "g2avetid_time_adjusted_score"
}

# Add in time gap from pretreatment hba1c  [time frame matched to metformin score]
if ("prt_hba1c_timegap"%in%c(keep,pheno)) {
  # read in data file
  a1c.data <- read.table(file.path("../pheno_data/hba1c_time_adjusted.txt"),header=TRUE,sep="\t")
  sub.a1c <- a1c.data[na.omit(match(data[,"MaskID"],a1c.data[,"MaskID"])),]
  data <- cbind(data,sub.a1c[,"prt_hba1c_timegap"])
  colnames(data)[ncol(data)] <- "prt_hba1c_timegap"
}



# Additional information in case we want to keep samples with some missing data
cat("\n# of missing covariates and corresponding # of samples")
table(apply(data,1,function(x) sum(is.na(x))))
cat("\n# of missing samples for each covariate\n")
apply(data,2,function(x) sum(is.na(x)))

# Mean-impute yrsdiab and yrslipi
if(!is.null(data$yrsdiab)) data[is.na(data$yrsdiab),names(data)=="yrsdiab"]<-mean(data$yrsdiab,na.rm=T)
if(!is.null(data$yrslipi))data[is.na(data$yrslipi),names(data)=="yrslipi"]<-mean(data$yrslipi,na.rm=T)

# Remove samples with missing data
data<-data[complete.cases(data),]

#head(data)

#i <- 5

# Needs to be fixed.  Causes errors.
for(i in 1:ncol(data)){
  if(length(unique(data[,i]))==1){
    # cat(names(data)[i], "\n")
    if(!file.exists(file.path(path,"remove_covars.txt"))) file.create(file.path(path,"remove_covars.txt"))
    tmp.remove <- unlist(read.table(file.path(path,"remove_covars.txt")))
    if (!any(tmp.remove==colnames(data)[i])){
        cat(paste0("#remove due to only one value across subjects\n",colnames(data)[i],"\n"),file=file.path(path,"remove_covars.txt"),append=TRUE)
    }
  }
}

# Plot correlation matrix
n<-dim(data)[2]-2
ColorRamp<-rgb(c(seq(0,1,length=21),seq(1,1,length=20)),       # Red
               c(seq(0,1,length=21),seq(1-1/20,0,length=20)),  # Green
               c(seq(1,1,length=21),seq(1-1/20,0,length=20)))  # Blue
n_color<-length(ColorRamp)
ColorLevels<-seq(-1,1,length=length(ColorRamp))
pdf(paste0(path,"/outputs/covar_correlation.pdf"),width=7, height=6)
layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(6,1), heights=c(1,1))
# Data Map
par(mar = c(4,4,.5,.5))
image(1:n,1:n,suppressWarnings(cor(data.matrix(data[,-c(1,2)]))),
    col=ColorRamp,xlab="",ylab="",axes=FALSE,zlim=c(-1,1))
for (i in 1:(n+1)) {abline(v=i-0.5)}
for (i in 1:(n+1)) {abline(h=i-0.5)}
axis(1,at=1:n,tick=F,labels=names(data)[-c(1,2)],cex.axis=0.5,las=2,pos=1.5)
axis(2,at=1:n,tick=F,labels=names(data)[-c(1,2)],cex.axis=0.5,las=2,pos=1.8)
# Color Scale
par(mar = c(4,2,.5,.5))
image(1,ColorLevels,matrix(data=ColorLevels,ncol=n_color,nrow=1),
      col=ColorRamp,xlab="",ylab="",yaxt="n",xaxt="n")
axis(2,at=ColorLevels,tick=F,labels=format(ColorLevels,digits=2),las=1,cex.axis=0.7,pos=0.8)
layout(1)
trash<-dev.off()

# Print highly correlated pairs
cor_matrix<-suppressWarnings(cor(data.matrix(data[,-c(1,2)])))
cor_matrix<-cor_matrix*lower.tri(cor_matrix)
idx<-which(abs(cor_matrix)>0.5,arr.ind=T)
cat("\ncorrelated variables\n")
if (dim(idx)[1]>0) {
    for (i in 1:dim(idx)[1]) {
        ii<-idx[i,1]
        jj<-idx[i,2]
        cat(names(data)[ii+2],"\t",names(data)[jj+2],"\t",cor_matrix[ii,jj],"\n")
    }
}

# Write output
write.table(data,paste0(path,"/pheno_data/pheno_data_step1.txt"),col.names=T,row.names=F,quote=F,sep="\t")
