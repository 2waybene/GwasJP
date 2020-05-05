 ### update original pheno_data.txt file
rm(list=ls())
#setwd("Y://NCSU_Projects//Accord//Data")
#setwd("/Volumes/Kirk 1/NCSU_Projects/Accord/Data")
#setwd("~/NCSU_Projects/Accord/Data")
#setwd("~/NCSU_Projects/Accord/accord_data/")
# setwd("~/data/accord/data")

p.data <- read.table(file.path("pheno_data","original_pheno_data_backup_v2.txt"),stringsAsFactors=FALSE,sep="\t",header=TRUE)
o.data <- p.data
#columns to keep from original
keep.cols <- c("MaskID","LabID","baseline_age","gender","arm","network","edu","smoking","alcohol","yrsdiab","hyptens","cvd_hx_baseline","x3malb","x3lvh","x3sten","eyedisea","neuropat","bmi","waist_cm","yrslipi","ethnicity")
o.data <- o.data[,which(colnames(o.data)%in%keep.cols)]


#add var for intensive arms
o.data <- cbind(o.data,as.numeric(p.data$arm%in%c(3,4,7,8)))
colnames(o.data)[ncol(o.data)] <- "int_gly_arm"
o.data <- cbind(o.data,as.numeric(p.data$arm%in%c(3,1)))
colnames(o.data)[ncol(o.data)] <- "int_bp_arm"
o.data <- cbind(o.data,as.numeric(p.data$arm%in%c(7,5)))
colnames(o.data)[ncol(o.data)] <- "fib_arm"

#########################################################################
### METFORMIN ###
met.scores <- read.csv(file.path("pheno_data","metformin_phenos","90d","metforminscores.csv"),header=TRUE,stringsAsFactors=FALSE)
o.data <- cbind(o.data,met.scores[match(o.data[,"MaskID"],met.scores[,"MaskID"]),grep("Score",colnames(met.scores))])
colnames(o.data)[ncol(o.data)] <- "metformin_score"

#incorporate other drug data for fibrate
d.files <- list.files(file.path("pheno_data","metformin_phenos","90d"))
#d.files <- d.files[grep("90d",d.files)]
d.files <- d.files[grep("metformin",d.files,invert=TRUE)]
#d.files <- d.files[grep("hba1cscores",d.files,invert=TRUE)]
d.files <- d.files[grep("scores",d.files)]
#d.files <- d.files[grep("phenos",d.files,invert=TRUE)]
for(d in d.files){
  if(grepl(".csv",d)){
    d.scores <- read.csv(file.path("pheno_data","metformin_phenos","90d",d),header=TRUE,stringsAsFactors=FALSE)
  }else{
    d.scores <- read.table(file.path("pheno_data","metformin_phenos","90d",d),header=TRUE,stringsAsFactors=FALSE)
  }
  
  score.cols <- grep("Score",colnames(d.scores))
  if(length(score.cols)==0)score.cols <- grep("hba1c",colnames(d.scores))

  if(colnames(d.scores[score.cols])=="Score"){
    colnames(d.scores)[score.cols] <- paste0("met_",unlist(strsplit(d,"scores"))[1],"Score")
  }else{
    colnames(d.scores)[score.cols] <- paste0("met_",colnames(d.scores)[score.cols])
  }
  
  
  o.data <- cbind(o.data,d.scores[match(o.data[,"MaskID"],d.scores[,"MaskID"]),score.cols])
  if(length(score.cols)==1) colnames(o.data)[ncol(o.data)] <- colnames(d.scores)[score.cols]
}



#incorporate other lab data for Metformin

o.lab.met <- read.table(file.path("pheno_data","metformin_phenos","90d","other_labs.txt"),header=TRUE,stringsAsFactors=FALSE)
o.data.cols <- which(colnames(o.lab.met)%in%c("fpg","screat","gfr"))
#colnames(o.lab.met)[which(colnames(o.lab.met)=="hba1c")] <- "pre_hba1c"
colnames(o.lab.met)[o.data.cols] <- paste0("met_",colnames(o.lab.met)[o.data.cols])
o.data <- cbind(o.data,o.lab.met[match(o.data[,"MaskID"],o.lab.met[,"MaskID"]),o.data.cols])

#incorporate BP data for Metformin
bp.met <- read.table(file.path("pheno_data","metformin_phenos","90d","bp.txt"),header=TRUE,stringsAsFactors=FALSE)
bp.cols <- which(colnames(bp.met)%in%c("sbp","dbp"))
colnames(bp.met)[bp.cols] <- paste0("met_",colnames(bp.met)[bp.cols])
o.data <- cbind(o.data,bp.met[match(o.data[,"MaskID"],bp.met[,"MaskID"]),bp.cols])


#incorporate insulin for Metformin
insulin <- read.table(file.path("pheno_data","metformin_phenos","90d","insulin.txt"),header=T,stringsAsFactors=FALSE)
ins.data.cols <- which(colnames(insulin)=="insulinScore")
o.data <- cbind(o.data,insulin[match(o.data[,"MaskID"],insulin[,"MaskID"]),ins.data.cols])
colnames(o.data)[ncol(o.data)] <- "met_insulin"

write.table(o.data,file.path("analysis","pheno_data_metformin.txt"),row.names=FALSE,sep="\t",quote=FALSE)
#########################################################################

### FIBRATE ###

fib.scores <- read.table(file.path("pheno_data","fibrate_phenos","90d","fibrate_score90d.txt"),header=TRUE,stringsAsFactors=FALSE)
o.data <- cbind(o.data,fib.scores[match(o.data[,"MaskID"],fib.scores[,"MaskID"]),grep("Score",colnames(fib.scores))])
colnames(o.data)[ncol(o.data)] <- "fibrate_score"

#incorporate other drug data for fibrate
d.files <- list.files(file.path("pheno_data","fibrate_phenos","90d"))
d.files <- d.files[grep("90d",d.files)]
d.files <- d.files[grep("fibrate",d.files,invert=TRUE)]
d.files <- d.files[grep("phenos",d.files,invert=TRUE)]
for(d in d.files){
  if(grepl(".csv",d)){
    d.scores <- read.csv(file.path("pheno_data","fibrate_phenos","90d",d),header=TRUE,stringsAsFactors=FALSE)
  }else{
    d.scores <- read.table(file.path("pheno_data","fibrate_phenos","90d",d),header=TRUE,stringsAsFactors=FALSE)
  }
  
  score.cols <- grep("Score",colnames(d.scores))
  
  if(colnames(d.scores[score.cols])=="Score"){
    colnames(d.scores)[score.cols] <- paste0("fib_",unlist(strsplit(d,"scores"))[1],"Score")
  }else{
    colnames(d.scores)[score.cols] <- paste0("fib_",colnames(d.scores)[score.cols])
  }
  
  
  o.data <- cbind(o.data,d.scores[match(o.data[,"MaskID"],d.scores[,"MaskID"]),score.cols])
  if(length(score.cols)==1) colnames(o.data)[ncol(o.data)] <- colnames(d.scores)[score.cols]
}



#incorporate other lab data for fibrate
o.lab.fib <- read.table(file.path("pheno_data","fibrate_phenos","90d","other_labs.txt"),header=TRUE,stringsAsFactors=FALSE)
o.data.cols <- which(colnames(o.lab.fib)%in%c("hba1c","fpg","screat","gfr"))
colnames(o.lab.fib)[o.data.cols] <- paste0("fib_",colnames(o.lab.fib)[o.data.cols])
o.data <- cbind(o.data,o.lab.fib[match(o.data[,"MaskID"],o.lab.fib[,"MaskID"]),o.data.cols])

#incorporate BP data for fibrate
bp.fib <- read.table(file.path("pheno_data","fibrate_phenos","90d","bp.txt"),header=TRUE,stringsAsFactors=FALSE)
bp.cols <- which(colnames(bp.fib)%in%c("sbp","dbp"))
colnames(bp.fib)[bp.cols] <- paste0("fib_",colnames(bp.fib)[bp.cols])
o.data <- cbind(o.data,bp.fib[match(o.data[,"MaskID"],bp.fib[,"MaskID"]),bp.cols])

#incorporate phenotype data for fibrate
pheno.fib <- read.table(file.path("pheno_data","fibrate_phenos","90d","lipid_phenos90d.txt"),header=TRUE,stringsAsFactors=FALSE)
pheno.data.cols <- which(colnames(pheno.fib)%in%c("pre_ldl","pre_hdl","pre_trig","pre_chol","ldl_lr","hdl_lr","trig_lr","chol_lr"))
o.data <- cbind(o.data,pheno.fib[match(o.data[,"MaskID"],pheno.fib[,"MaskID"]),pheno.data.cols])

write.table(o.data,file.path("analysis","pheno_data_fibrate.txt"),row.names=FALSE,sep="\t",quote=FALSE)
#########################################################################


### TZD ###
fib.scores <- read.csv(file.path("pheno_data","tzd_phenos","90d","tzdscores.csv"),header=TRUE,stringsAsFactors=FALSE)
o.data <- cbind(o.data,fib.scores[match(o.data[,"MaskID"],fib.scores[,"MaskID"]),grep("Score",colnames(fib.scores))])
colnames(o.data)[ncol(o.data)] <- "tzd_score"

#incorporate other drug data for fibrate
d.files <- list.files(file.path("pheno_data","tzd_phenos","90d"))
#d.files <- d.files[grep("90d",d.files)]
d.files <- d.files[grep("tzd",d.files,invert=TRUE)]
d.files <- d.files[grep("scores",d.files)]
#d.files <- d.files[grep("phenos",d.files,invert=TRUE)]
for(d in d.files){
  if(grepl(".csv",d)){
    d.scores <- read.csv(file.path("pheno_data","tzd_phenos","90d",d),header=TRUE,stringsAsFactors=FALSE)
  }else{
    d.scores <- read.table(file.path("pheno_data","tzd_phenos","90d",d),header=TRUE,stringsAsFactors=FALSE)
  }
  
  score.cols <- grep("Score",colnames(d.scores))
  if(length(score.cols)==0)score.cols <- grep("hba1c",colnames(d.scores))
  
  
  if(colnames(d.scores[score.cols])=="Score"){
    colnames(d.scores)[score.cols] <- paste0("tzd_",unlist(strsplit(d,"scores"))[1],"Score")
  }else{
    colnames(d.scores)[score.cols] <- paste0("tzd_",colnames(d.scores)[score.cols])
  }
  
  
  o.data <- cbind(o.data,d.scores[match(o.data[,"MaskID"],d.scores[,"MaskID"]),score.cols])
  if(length(score.cols)==1) colnames(o.data)[ncol(o.data)] <- colnames(d.scores)[score.cols]
}



#incorporate other lab data for TZD

o.lab.fib <- read.table(file.path("pheno_data","tzd_phenos","90d","other_labs.txt"),header=TRUE,stringsAsFactors=FALSE)
o.data.cols <- which(colnames(o.lab.fib)%in%c("fpg","screat","gfr"))
#colnames(o.lab.fib)[which(colnames(o.lab.fib)=="hba1c")] <- "pre_hba1c"
colnames(o.lab.fib)[o.data.cols] <- paste0("tzd_",colnames(o.lab.fib)[o.data.cols])
o.data <- cbind(o.data,o.lab.fib[match(o.data[,"MaskID"],o.lab.fib[,"MaskID"]),o.data.cols])

#incorporate BP data for TZD
bp.fib <- read.table(file.path("pheno_data","tzd_phenos","90d","bp.txt"),header=TRUE,stringsAsFactors=FALSE)
bp.cols <- which(colnames(bp.fib)%in%c("sbp","dbp"))
colnames(bp.fib)[bp.cols] <- paste0("tzd_",colnames(bp.fib)[bp.cols])
o.data <- cbind(o.data,bp.fib[match(o.data[,"MaskID"],bp.fib[,"MaskID"]),bp.cols])

#incorporate phenotype data for TZD### ALREADY INCORPORATED IN D.FILES LOOP ABOVE
#pheno.fib <- read.csv(file.path("pheno_data","tzd_phenos","90d","hba1cscores.csv"),stringsAsFactors=FALSE)
#pheno.data.cols <- which(colnames(pheno.fib)=="hba1c")
#o.data <- cbind(o.data,pheno.fib[match(o.data[,"MaskID"],pheno.fib[,"MaskID"]),pheno.data.cols])
#colnames(o.data)[ncol(o.data)] <- "tzd_hba1c"

#incorporate insulin for TZD
insulin <- read.table(file.path("pheno_data","tzd_phenos","90d","insulin.txt"),header=T,stringsAsFactors=FALSE)
ins.data.cols <- which(colnames(insulin)=="insulinScore")
o.data <- cbind(o.data,insulin[match(o.data[,"MaskID"],insulin[,"MaskID"]),ins.data.cols])
colnames(o.data)[ncol(o.data)] <- "tzd_insulin"

write.table(o.data,file.path("analysis","pheno_data_tzd.txt"),row.names=FALSE,sep="\t",quote=FALSE)

#########################################################################
### TZD Weight Gain  ###

# Cluster filepath
phenopath <- file.path("pheno_data","tzd_phenos","12mo")
analysispath <- file.path("/home", "accord", "data", "analysis")

# Import 12month weight gain scores for TZD and add to o.data
tzd.scores <- read.csv(file.path(phenopath,"tzdscores_weight_gain.csv"),header = TRUE, stringsAsFactors = FALSE)
o.data <- cbind(o.data, tzd.scores[match(o.data[,"MaskID"], tzd.scores[,"MaskID"]), grep("Score", colnames(tzd.scores))])
colnames(o.data)[ncol(o.data)] <- "tzd_score_weight_gain"

#incorporate other drug data for TZD
d.files <- list.files(phenopath)
#d.files <- d.files[grep("90d",d.files)]
d.files <- d.files[grep("tzd",d.files,invert=TRUE)]
d.files <- d.files[grep("age",d.files,invert=TRUE)]
d.files <- d.files[grep("scores",d.files)]
#d.files <- d.files[grep("phenos",d.files,invert=TRUE)]
for(d in d.files){
  if(grepl(".csv",d)){
    d.scores <- read.csv(file.path(phenopath, d),header = TRUE, stringsAsFactors = FALSE)
  }else{
    d.scores <- read.table(file.path(phenopath, d),header = TRUE, stringsAsFactors = FALSE)
  }
  
  score.cols <- grep("Score",colnames(d.scores))
  # weight is the only file where the column of interst doesn't say "Score"
  if(length(score.cols)==0)score.cols <- grep("wt_kg",colnames(d.scores))
  
  
  if(colnames(d.scores[score.cols])=="Score"){
    colnames(d.scores)[score.cols] <- paste0("tzd_",unlist(strsplit(d,"scores"))[1],"Score")
  }else{
    colnames(d.scores)[score.cols] <- paste0("tzd_",colnames(d.scores)[score.cols])
  }
  
  
  o.data <- cbind(o.data,d.scores[match(o.data[,"MaskID"],d.scores[,"MaskID"]),score.cols])
  if(length(score.cols)==1) colnames(o.data)[ncol(o.data)] <- colnames(d.scores)[score.cols]
}

# Add in Age Starting TZD  
age <- read.csv(file.path(phenopath,"agescores.csv"),header=T,stringsAsFactors=FALSE)
ins.data.cols <- which(colnames(age)=="age_start_tzd")
o.data <- cbind(o.data,age[match(o.data[,"MaskID"],age[,"MaskID"]),ins.data.cols])
colnames(o.data)[ncol(o.data)] <- "tzd_age_start_tzd"

#incorporate insulin for TZD
insulin <- read.table(file.path(phenopath,"insulin.txt"),header=T,stringsAsFactors=FALSE)
ins.data.cols <- which(colnames(insulin)=="insulinScore")
o.data <- cbind(o.data,insulin[match(o.data[,"MaskID"],insulin[,"MaskID"]),ins.data.cols])
colnames(o.data)[ncol(o.data)] <- "tzd_insulin"



write.table(o.data, file.path(analysispath, "pheno_data_tzd_weight_gain.txt"),row.names = FALSE, sep = "\t", quote = FALSE)

#########################################################################




### STATIN ###

fib.scores <- read.table(file.path("pheno_data","statin_phenos","120d","statinscores.txt"),header=TRUE,stringsAsFactors=FALSE,sep="\t")
o.data <- cbind(o.data,fib.scores[match(o.data[,"MaskID"],fib.scores[,"MaskID"]),grep("Score",colnames(fib.scores))])
colnames(o.data)[ncol(o.data)] <- "statin_score"

#incorporate other drug data for fibrate
d.files <- list.files(file.path("pheno_data","statin_phenos","120d"))
#d.files <- d.files[grep("90d",d.files)]
d.files <- d.files[grep("statin",d.files,invert=TRUE)]
d.files <- d.files[grep("scores",d.files)]
#d.files <- d.files[grep("phenos",d.files,invert=TRUE)]
for(d in d.files){
  if(grepl(".csv",d)){
    d.scores <- read.csv(file.path("pheno_data","statin_phenos","120d",d),header=TRUE,stringsAsFactors=FALSE)
  }else{
    d.scores <- read.table(file.path("pheno_data","statin_phenos","120d",d),header=TRUE,stringsAsFactors=FALSE)
  }
  
  score.cols <- grep("Score",colnames(d.scores))
  if(colnames(d.scores[score.cols])=="Score"){
    colnames(d.scores)[score.cols] <- paste0("stat_",unlist(strsplit(d,"scores"))[1],"Score")
  }else{
    colnames(d.scores)[score.cols] <- paste0("stat_",colnames(d.scores)[score.cols])
  }
  
  
 
  
  o.data <- cbind(o.data,d.scores[match(o.data[,"MaskID"],d.scores[,"MaskID"]),score.cols])
  if(length(score.cols)==1) colnames(o.data)[ncol(o.data)] <- colnames(d.scores)[score.cols]
}



#incorporate other lab data for STATIN

o.lab.fib <- read.table(file.path("pheno_data","statin_phenos","120d","other_labs.txt"),header=TRUE,stringsAsFactors=FALSE)
o.data.cols <- which(colnames(o.lab.fib)%in%c("hba1c","fpg","screat","gfr"))
colnames(o.lab.fib)[o.data.cols] <- paste0("stat_",colnames(o.lab.fib)[o.data.cols])
o.data <- cbind(o.data,o.lab.fib[match(o.data[,"MaskID"],o.lab.fib[,"MaskID"]),o.data.cols])

#incorporate BP data for fibrate
bp.fib <- read.table(file.path("pheno_data","statin_phenos","120d","bp.txt"),header=TRUE,stringsAsFactors=FALSE)
bp.cols <- which(colnames(bp.fib)%in%c("sbp","dbp"))
colnames(bp.fib)[bp.cols] <- paste0("stat_",colnames(bp.fib)[bp.cols])
o.data <- cbind(o.data,bp.fib[match(o.data[,"MaskID"],bp.fib[,"MaskID"]),bp.cols])

#incorporate phenotype data for fibrate
pheno.fib <- read.table(file.path("pheno_data","statin_phenos","120d","lipid_phenos120d.txt"),header=TRUE,stringsAsFactors=FALSE)
pheno.data.cols <- which(colnames(pheno.fib)%in%c("pre_ldl","pre_hdl","pre_trig","pre_chol","ldl_lr","hdl_lr","trig_lr","chol_lr"))
o.data <- cbind(o.data,pheno.fib[match(o.data[,"MaskID"],pheno.fib[,"MaskID"]),pheno.data.cols])

write.table(o.data,file.path("analysis","pheno_data_statin.txt"),row.names=FALSE,sep="\t",quote=FALSE)
#########################################################################

### Sulfonylurea ###
fib.scores <- read.csv(file.path("pheno_data","sulfonylurea_phenos","90d","sulfscores.csv"),header=TRUE,stringsAsFactors=FALSE)
o.data <- cbind(o.data,fib.scores[match(o.data[,"MaskID"],fib.scores[,"MaskID"]),grep("Score",colnames(fib.scores))])
colnames(o.data)[ncol(o.data)] <- "sulf_score"

#incorporate other drug data for fibrate
d.files <- list.files(file.path("pheno_data","sulfonylurea_phenos","90d"))
#d.files <- d.files[grep("90d",d.files)]
d.files <- d.files[grep("sulf",d.files,invert=TRUE)]
d.files <- d.files[grep("scores",d.files)]
#d.files <- d.files[grep("phenos",d.files,invert=TRUE)]
for(d in d.files){
  if(grepl(".csv",d)){
    d.scores <- read.csv(file.path("pheno_data","sulfonylurea_phenos","90d",d),header=TRUE,stringsAsFactors=FALSE)
  }else{
    d.scores <- read.table(file.path("pheno_data","sulfonylurea_phenos","90d",d),header=TRUE,stringsAsFactors=FALSE)
  }
  
  score.cols <- grep("Score",colnames(d.scores))
  #print(d)
  #print(score.cols)
  if(length(score.cols)==0)score.cols <- grep("hba1c",colnames(d.scores))
  
  
  if(colnames(d.scores[score.cols])=="Score"){
    colnames(d.scores)[score.cols] <- paste0("sulf_",unlist(strsplit(d,"scores"))[1],"Score")
  }else{
    colnames(d.scores)[score.cols] <- paste0("sulf_",colnames(d.scores)[score.cols])
  }
  
  
  o.data <- cbind(o.data,d.scores[match(o.data[,"MaskID"],d.scores[,"MaskID"]),score.cols])
  if(length(score.cols)==1) colnames(o.data)[ncol(o.data)] <- colnames(d.scores)[score.cols]
}



#incorporate other lab data for sulfonylurea

o.lab.fib <- read.table(file.path("pheno_data","sulfonylurea_phenos","90d","other_labs.txt"),header=TRUE,stringsAsFactors=FALSE)
o.data.cols <- which(colnames(o.lab.fib)%in%c("fpg","screat","gfr"))
#colnames(o.lab.fib)[which(colnames(o.lab.fib)=="hba1c")] <- "pre_hba1c"
colnames(o.lab.fib)[o.data.cols] <- paste0("sulf_",colnames(o.lab.fib)[o.data.cols])
o.data <- cbind(o.data,o.lab.fib[match(o.data[,"MaskID"],o.lab.fib[,"MaskID"]),o.data.cols])

#incorporate BP data for sulfonylurea
bp.fib <- read.table(file.path("pheno_data","sulfonylurea_phenos","90d","bp.txt"),header=TRUE,stringsAsFactors=FALSE)
bp.cols <- which(colnames(bp.fib)%in%c("sbp","dbp"))
colnames(bp.fib)[bp.cols] <- paste0("sulf_",colnames(bp.fib)[bp.cols])
o.data <- cbind(o.data,bp.fib[match(o.data[,"MaskID"],bp.fib[,"MaskID"]),bp.cols])

#incorporate phenotype data for sulfonylurea### ALREADY INCORPORATED IN D.FILES LOOP ABOVE
#pheno.fib <- read.csv(file.path("pheno_data","tzd_phenos","90d","hba1cscores.csv"),stringsAsFactors=FALSE)
#pheno.data.cols <- which(colnames(pheno.fib)=="hba1c")
#o.data <- cbind(o.data,pheno.fib[match(o.data[,"MaskID"],pheno.fib[,"MaskID"]),pheno.data.cols])
#colnames(o.data)[ncol(o.data)] <- "tzd_hba1c"

#incorporate insulin for sulfonylurea
insulin <- read.table(file.path("pheno_data","sulfonylurea_phenos","90d","insulin.txt"),header=T,stringsAsFactors=FALSE)
ins.data.cols <- which(colnames(insulin)=="insulinScore")
o.data <- cbind(o.data,insulin[match(o.data[,"MaskID"],insulin[,"MaskID"]),ins.data.cols])
colnames(o.data)[ncol(o.data)] <- "sulf_insulin"

write.table(o.data,file.path("analysis","pheno_data_sulf.txt"),row.names=FALSE,sep="\t",quote=FALSE)




#########################################################################
### STATIN for drug scoring paper ###

fib.scores <- read.csv(file.path("pheno_data","statin_phenos","120d_drugScorePaper","statinscores.csv"))
o.data <- cbind(o.data,fib.scores[match(o.data[,"MaskID"],fib.scores[,"MaskID"]),grep("Score",colnames(fib.scores))])
colnames(o.data)[ncol(o.data)] <- "statin_score"

#incorporate other drug data for fibrate
d.files <- list.files(file.path("pheno_data","statin_phenos","120d_drugScorePaper"))
#d.files <- d.files[grep("90d",d.files)]
d.files <- d.files[grep("statin",d.files,invert=TRUE)]
d.files <- d.files[grep("scores",d.files)]
#d.files <- d.files[grep("phenos",d.files,invert=TRUE)]
for(d in d.files){
  if(grepl(".csv",d)){
    d.scores <- read.csv(file.path("pheno_data","statin_phenos","120d_drugScorePaper",d),header=TRUE,stringsAsFactors=FALSE)
  }else{
    d.scores <- read.table(file.path("pheno_data","statin_phenos","120d_drugScorePaper",d),header=TRUE,stringsAsFactors=FALSE)
  }
  
  score.cols <- grep("Score",colnames(d.scores))
  if(colnames(d.scores[score.cols])=="Score"){
    colnames(d.scores)[score.cols] <- paste0("stat_",unlist(strsplit(d,"scores"))[1],"Score")
  }else{
    colnames(d.scores)[score.cols] <- paste0("stat_",colnames(d.scores)[score.cols])
  }
  
  
  o.data <- cbind(o.data,d.scores[match(o.data[,"MaskID"],d.scores[,"MaskID"]),score.cols])
  if(length(score.cols)==1) colnames(o.data)[ncol(o.data)] <- colnames(d.scores)[score.cols]
}



#incorporate other lab data for STATIN

o.lab.fib <- read.table(file.path("pheno_data","statin_phenos","120d_drugScorePaper","other_labs.txt"),header=TRUE,stringsAsFactors=FALSE)
o.data.cols <- which(colnames(o.lab.fib)%in%c("hba1c","fpg","screat","gfr"))
colnames(o.lab.fib)[o.data.cols] <- paste0("stat_",colnames(o.lab.fib)[o.data.cols])
o.data <- cbind(o.data,o.lab.fib[match(o.data[,"MaskID"],o.lab.fib[,"MaskID"]),o.data.cols])

#incorporate BP data for fibrate
bp.fib <- read.table(file.path("pheno_data","statin_phenos","120d_drugScorePaper","bp.txt"),header=TRUE,stringsAsFactors=FALSE)
bp.cols <- which(colnames(bp.fib)%in%c("sbp","dbp"))
colnames(bp.fib)[bp.cols] <- paste0("stat_",colnames(bp.fib)[bp.cols])
o.data <- cbind(o.data,bp.fib[match(o.data[,"MaskID"],bp.fib[,"MaskID"]),bp.cols])

#incorporate phenotype data for fibrate
pheno.fib <- read.table(file.path("pheno_data","statin_phenos","120d_drugScorePaper","lipid_phenos120d.txt"),header=TRUE,stringsAsFactors=FALSE)
pheno.data.cols <- which(colnames(pheno.fib)%in%c("pre_ldl","pre_hdl","pre_trig","pre_chol","ldl_lr","hdl_lr","trig_lr","chol_lr"))
o.data <- cbind(o.data,pheno.fib[match(o.data[,"MaskID"],pheno.fib[,"MaskID"]),pheno.data.cols])

write.table(o.data,file.path("analysis","pheno_data_statin_drugScorePaper.txt"),row.names=FALSE,sep="\t",quote=FALSE)


#########################################################################
### Dieatary phenos ###
#incorporate other drug data for fibrate
d.files <- read.csv(file.path("pheno_data","dietary_scores.csv"))
o.data <- merge(o.data,d.files,by="MaskID",all.x = TRUE)
write.table(o.data,file.path("analysis","pheno_data_diet.txt"),row.names=FALSE,sep="\t",quote=FALSE)


#########################################################################
### GV phenos ###
getwd()
gv.scores <- read.csv(file.path("pheno_data","GV_phenos","GVscores.csv"))
o.data <- merge(o.data,gv.scores, by="MaskID",all.x = TRUE)


#incorporate other lab data for GV
o.lab<- read.table(file.path("pheno_data","otherlabs.txt"),header=TRUE,stringsAsFactors=FALSE)
o.lab <- o.lab[o.lab$Visit=="BLR",]
o.lab <- o.lab[,c("MaskID","fpg","screat","gfr")]
hba1c <- read.table(file.path("pheno_data","hba1c.txt"),header=TRUE,stringsAsFactors=FALSE)
hba1c <- hba1c[hba1c$Visit=="BLR",]
o.data <- merge(o.data,o.lab,by="MaskID",all.x=TRUE)
o.data <- merge(o.data,hba1c[,c("MaskID","hba1c")],by="MaskID",all.x=TRUE)

#incorporate BP data for GV
bp <- read.table(file.path("pheno_data","bloodpressure.txt"),header=TRUE,stringsAsFactors=FALSE)
bp <- bp[bp$Visit=="BLR",]
bp <- bp[,c("MaskID","sbp","dbp")]
o.data <- merge(o.data,bp,by="MaskID",all.x=TRUE)



cvds <- read.table(file='pheno_data/cvdoutcomes.txt',header=T,skip=33,sep='\t')
#head(cvds)
#table(cvds$fuyrs_po7p)
#length(which(is.na(cvds$fuyrs_po7p)))
colnames(cvds)[c(1,grep('censor',colnames(cvds)))]

cvds <- cvds[,c(1,grep('censor',colnames(cvds)))]
#head(cvds)

i <- 1
## Switch 0s and 1s for outcomes
for (i in 2:ncol(cvds)){
  ## Check NAs
  print(length(which(is.na(cvds[,i]))))
  print(table(cvds[,i]))
  cvds[which(cvds[,i]==1),i] <- -1
  cvds[which(cvds[,i]==0),i] <- 1
  cvds[which(cvds[,i]==-1),i] <- 0
  print(table(cvds[,i]))
}



o.data <- merge(o.data,cvds,by.x="MaskID",all.x=TRUE)
#head(dat)
#dim(dat)

write.table(o.data,file.path("analysis","pheno_data_GV.txt"),row.names=FALSE,sep="\t",quote=FALSE)


#########################################################################
### HGI phenos: hgi gets calculated in model_step1 ###

#incorporate other lab data for HGI
o.lab<- read.table(file.path("pheno_data","otherlabs.txt"),header=TRUE,stringsAsFactors=FALSE)
o.lab <- o.lab[o.lab$Visit=="BLR",]
o.lab <- o.lab[,c("MaskID","fpg","screat","gfr")]
hba1c <- read.table(file.path("pheno_data","hba1c.txt"),header=TRUE,stringsAsFactors=FALSE)
hba1c <- hba1c[hba1c$Visit=="BLR",]
o.data <- merge(o.data,o.lab,by="MaskID",all.x=TRUE)
o.data <- merge(o.data,hba1c[,c("MaskID","hba1c")],by="MaskID",all.x=TRUE)

#incorporate BP data for HGI
bp <- read.table(file.path("pheno_data","bloodpressure.txt"),header=TRUE,stringsAsFactors=FALSE)
bp <- bp[bp$Visit=="BLR",]
bp <- bp[,c("MaskID","sbp","dbp")]
o.data <- merge(o.data,bp,by="MaskID",all.x=TRUE)

write.table(o.data,file.path("analysis","pheno_data_hgi.txt"),row.names=FALSE,sep="\t",quote=FALSE)


#########################################################################
### ICAPS: data from michael

#incorporate other lab data for GV
o.lab<- read.table(file.path("pheno_data","otherlabs.txt"),header=TRUE,stringsAsFactors=FALSE)
o.lab <- o.lab[o.lab$Visit=="BLR",]
o.lab <- o.lab[,c("MaskID","fpg","screat","gfr")]
hba1c <- read.table(file.path("pheno_data","hba1c.txt"),header=TRUE,stringsAsFactors=FALSE)
hba1c <- hba1c[hba1c$Visit=="BLR",]
o.data <- merge(o.data,o.lab,by="MaskID",all.x=TRUE)
o.data <- merge(o.data,hba1c[,c("MaskID","hba1c")],by="MaskID",all.x=TRUE)

#incorporate BP data for GV
bp <- read.table(file.path("pheno_data","bloodpressure.txt"),header=TRUE,stringsAsFactors=FALSE)
bp <- bp[bp$Visit=="BLR",]
bp <- bp[,c("MaskID","sbp","dbp")]
o.data <- merge(o.data,bp,by="MaskID",all.x=TRUE)

getwd()
icaps.scores <- read.table(file.path("pheno_data/icaps/","summary.bp.file.csv"),sep='\t',header=T)

## need to read this into create.sample.list.r
colnames(icaps.scores)[13:17] <- paste0(colnames(icaps.scores)[8:12],'.exclusion')

head(icaps.scores)
o.data <- merge(o.data,icaps.scores, by="MaskID",all.x = TRUE)

head(o.data[1:10])
dim(o.data)
dim(o.data[which(!is.na(o.data$Other.exclusion)),])
## 3667 individuals left
dim(icaps.scores)
tmp <- setdiff(icaps.scores$MaskID,o.data$MaskID)
length(tmp)
## There were 1065 ICAPS individuals not in pheno_data
head(o.data)
dim(o.data)

head(o.data[which(o.data$on_trial..6mo.==0),])
length(which(o.data$on_trial..6mo.==0))

## For anyone not on trial long enough, change PO to 0
o.data$ICAPS_po[which(o.data$on_trial..6mo.==0)] <- 0

## Do this in create.sample.list.r
#o.data <- o.data[which(!is.na(o.data$Other)),]
#dim(o.data)

write.table(o.data,file.path("analysis","pheno_data_icaps.txt"),row.names=FALSE,sep="\t",quote=FALSE)
write.table(o.data,file.path("analysis","pheno_data.txt"),row.names=FALSE,sep="\t",quote=FALSE)

colnames(o.data)

##############################################################
## pheno data for clinical clusters
pheno.data <- read.table("~/NCSU_Projects/Accord/accord_data/analysis/pheno_data_clustering.txt",header=TRUE,sep="\t")
clust.data <- read.csv("~/NCSU_Projects/Accord/Data/analysis/Clustering/martin_final_summary/ACCORD_Baseline_Phenotype_Analysis_Cluster.csv")

o.data <- merge(pheno.data,clust.data[,c("mid","A","B","C")],by.x="MaskID",by.y="mid",all.x=TRUE)
colnames(o.data)[which(colnames(o.data)%in%c("A","B","C"))] <- c("clust_1","clust_2","clust_3")
write.table(o.data,"~/NCSU_Projects/Accord/accord_data/analysis/pheno_data_clustering.txt",row.names=FALSE,sep="\t",quote=FALSE)


##############################################################
## pheno data for hba1c baseline
blr.hba1c <- read.table("pheno_data/hba1c.txt",header=TRUE)
blr.hba1c <- blr.hba1c[blr.hba1c$Visit=="BLR",-2]

#incorporate other labs
o.lab<- read.table(file.path("pheno_data","otherlabs.txt"),header=TRUE,stringsAsFactors=FALSE)
o.lab <- o.lab[o.lab$Visit=="BLR",]
o.lab <- o.lab[,c("MaskID","fpg","screat","gfr")]
o.data <- merge(o.data,o.lab,by="MaskID",all.x=TRUE)

#incorporate BP data
bp <- read.table(file.path("pheno_data","bloodpressure.txt"),header=TRUE,stringsAsFactors=FALSE)
bp <- bp[bp$Visit=="BLR",]
bp <- bp[,c("MaskID","sbp","dbp")]
o.data <- merge(o.data,bp,by="MaskID",all.x=TRUE)

#incorporate concomitant med data
med.data <- read.table(file.path("pheno_data","concomitantmeds.txt"),header=TRUE,stringsAsFactors=FALSE)
med.data <- med.data[med.data$Visit=="BLR",-2]
o.data <- merge(o.data,med.data,by="MaskID",all.x=TRUE)

o.data <- merge(o.data,blr.hba1c,by="MaskID",all.x=TRUE)

o.data[which(o.data==".",arr.ind=TRUE)] <- NA

write.table(o.data,"~/NCSU_Projects/Accord/accord_data/analysis/pheno_data_hba1c.txt",row.names=FALSE,sep="\t",quote=FALSE)

############################################################
## BEGIN MIND COG Jasmin

## Read in MIND Data
d.old = read.table('~/data/accord/data/pheno_data/mindcog/accord_analysis_file_13March2018.txt',header=T,fill=NA,sep='\t')
d = read.table('~/data/accord/data/pheno_data/mindcog/accord_analysis_file04_16_2018.txt',header=T,fill=NA,sep='\t')

colnames(d.old)
colnames(d)
head(d)
#table(d$lipid)
#table(d$arm)
#length(d$MaskID) # 492
#length(unique(d$MaskID)) # 492
#length(which(o.data$MaskID %in% d$MaskID)) #389
#length(which(d$MaskID %!in% o.data$MaskID )) #103
#length(which(d$MaskID %in% o.data$MaskID )) #389
#head(tmp.d[,c('MaskID','bp')])
#head(tmp.o[,c('MaskID','int_bp_arm')])
#head(tmp.d)

## Find Genotyped IDs in MIND
tmp.o = o.data[which(o.data$MaskID %in% d$MaskID),]
tmp.d = d[which(d$MaskID %in% o.data$MaskID),]
## Sort
tmp.o = tmp.o[order(tmp.o$MaskID),]
tmp.d = tmp.d[order(tmp.d$MaskID),]
dim(tmp.d)
colnames(tmp.d)
dim(tmp.o)
head(tmp.o)
head(tmp.d)

tmp.d = merge(tmp.o[,1:2],tmp.d,by=c("MaskID"))
head(tmp.d)
table(tmp.d$smoking)
table(tmp.d$edu)
table(tmp.d$cvd)
head(o.data)
table(o.data$cvd_hx_baseline)

## make column ckd (ACR > 30 or eGFR < 60)
tmp.d$ckd = ifelse((tmp.d$acr30==1) | (tmp.d$gfr60==1),1,0)

head(tmp.d)
colnames(tmp.d)

head(tmp.o)
tmp.d$ethnicity = "Unknown"
for (i in 1:nrow(tmp.d)){
	tmp.d$ethnicity[i] = tmp.o$ethnicity[which(tmp.o$MaskID==tmp.d$MaskID[i])]
}
table(tmp.d$ethnicity)

head(tmp.d)


## Set up interaction covariates
# The snps rs429358 and rs7412 are in the imputed dataset, chr19 chunk 10.
#rs429358  rs7412  Name
# C           T      ε1
# T           T      ε2
# T           C      ε3
# C           C      ε4
system('mkdir /tmp/accord.party.19.10')
#f_gz="chr$chr.imputed.$chunk"
system('zcat ~/data/accord/data/imputation/outputs/chr19.imputed.10.gz > /tmp/accord.party.19.10/chr19.imputed.10')
library(data.table)
apoe = fread('/tmp/accord.party.19.10/chr19.imputed.10')
dim(apoe)
head(apoe)[,1:10]
keep.i = which(apoe[,2] %in% c('rs7412','rs429358'))
keep.i = grep('rs7412',unlist(apoe[,2]),fixed=T)[1]
tmp = grep('rs429358',unlist(apoe[,2]),fixed=T)
keep.i = c(keep.i,grep('rs429358',unlist(apoe[,2]),fixed=T))
apoe[keep.i,1:10]
apoe = apoe[keep.i,]
dim(apoe)
gc()

## From compute_cv.r
chr<-19
chunk<-10
# Parameters
info_threshold<-0.5
read_size<-1000
source("~/data/accord/data/analysis/bin/load_pheno_data.r")
path = "~/data/accord/data/analysis/mindcog/mindcog_cog_dsst"
pheno_data<-load_pheno_data(path)


#time R --slave --vanilla --file=bin/compute_cv.r --args $p $chr $chunk "/tmp/accord.party.$me.$chr.$chunk.$tmprand/$f_gz"
#echo "Remove tmp directory on node: /tmp/accord.party.$me.$chr.$chunk.$tmprand"
#rm -r "/tmp/accord.party.$me.$chr.$chunk.$tmprand"

## Create sample list (bypassing bin/create.sample.list.r)
setwd('~/data/accord/data/analysis/mindcog/')
paths = c('mindcog_mri_wm','mindcog_mri_gm','mindcog_mri_wmlv',
          'mindcog_cog_mse','mindcog_cog_ravlt',
          'mindcog_cog_stroop','mindcog_cog_dsst')

## Write out the phenotype file with all data
write.table(tmp.d,file="pheno_data_mindcog.txt",row.names=FALSE,sep="\t",quote=FALSE)

getwd()
################################
### Write out snp_list
d = read.csv('~/data/accord/data/pheno_data/mindcog/AADHS_replication_snps_11302017.csv')
head(d)
for (f in paths){
	write.table(d[,1],file=paste0(f,"/snp_list.txt"),row.names=FALSE,sep="\t",quote=FALSE,col.names=F)
}


####### Outcomes:
#### MRI
## white matter (wm)
## gray matter (gm)
## white matter lesion volumes (wmlv)
## DO NOT HAVE THESE:
## hippocampal WM, hippocampal GM, and total hippocampal volume
head(tmp.d)

## white matter (wm)
dim(tmp.d[-which(is.na(tmp.d$wm)),])
tmp.save = tmp.d[-which(is.na(tmp.d$wm)),1:2]
#head(tmp.save)
print(paths[1])
write.table(file=paste0("mindcog/",paths[1],"/starting_samples.txt"),tmp.save,col.names=F,row.names=F,quote=F)

## gray matter (gm)
dim(tmp.d[-which(is.na(tmp.d$gm)),])
tmp.save = tmp.d[-which(is.na(tmp.d$gm)),1:2]
#head(tmp.save)
print(paths[2])
write.table(file=paste0("mindcog/",paths[2],"/starting_samples.txt"),tmp.save,col.names=F,row.names=F,quote=F)

## white matter lesion volumes (wmlv)
dim(tmp.d[-which(is.na(tmp.d$wmlv)),])
tmp.save = tmp.d[-which(is.na(tmp.d$wmlv)),1:2]
#head(tmp.save)
print(paths[3])
write.table(file=paste0("mindcog/",paths[3],"/starting_samples.txt"),tmp.save,col.names=F,row.names=F,quote=F)

#### Cognitive function measures:
## MSE (base_mmse)
## RAVLT (base_ravlt)
## Stroop (base_stroop)
## inference (base_stroop?)
## total digit substitution
head(tmp.d)

## MSE (base_mmse) 
## NONE MISSING
#dim(tmp.d[-which(is.na(tmp.d$base_mmse)),])
tmp.save = tmp.d[,1:2]
#head(tmp.save)
print(paths[4])
write.table(file=paste0("mindcog/",paths[4],"/starting_samples.txt"),tmp.save,col.names=F,row.names=F,quote=F)


## RAVLT (base_ravlt)
## NONE MISSING 
#dim(tmp.d[-which(is.na(tmp.d$base_ravlt)),])
tmp.save = tmp.d[,1:2]
#head(tmp.save)
print(paths[5])
write.table(file=paste0("mindcog/",paths[5],"/starting_samples.txt"),tmp.save,col.names=F,row.names=F,quote=F)

## Stroop (base_stroop)
dim(tmp.d[-which(is.na(tmp.d$base_stroop)),])
tmp.save = tmp.d[,1:2]
#head(tmp.save)
print(paths[6])
write.table(file=paste0("mindcog/",paths[6],"/starting_samples.txt"),tmp.save,col.names=F,row.names=F,quote=F)

## starting_covars.txt created manually
## forced_covars.txt created manually
## MRI: 
#Generalized linear models adjusting for ancestry proportion (or PCs) scanner strength (1.5T vs. 3.0T), 
#total intra-cranial volume (TICV), age, BMI, hba1c, statin use, smoking (past/current, never), 
#CVD, CKD (ACR>30 or eGFR<60), hypertension status. 

## Cognitive function: 
#GLM model with negative binomial link adjusting for ancestry proportion (or PCs), 
#education, age, BMI, hba1c, statin use, smoking (past/current, never), CVD, 
#CKD (defined as ACR>30 or eGFR<60), hypertension status. 

## remove_covars.txt created manually


######### Compare variables: UNFINISHED

########## They don't even use all these variables. EDU is redundant

'%!in%' <- function(x,y)!('%in%'(x,y))

head(o.data)
d = read.table('~/data/accord/data/pheno_data/mindcog/accord_analysis_file_13March2018.txt',header=T,fill=NA,sep='\t')
head(d)

table(d$lipid)
table(d$arm)

length(d$MaskID) # 492
length(unique(d$MaskID)) # 492
length(which(o.data$MaskID %in% d$MaskID)) #389
length(which(d$MaskID %!in% o.data$MaskID )) #103
length(which(d$MaskID %in% o.data$MaskID )) #389

## Subset data by MaskID from MIND-COG
tmp.o = o.data[which(o.data$MaskID %in% d$MaskID),]
tmp.d = d[which(d$MaskID %in% o.data$MaskID),]
# Should be 389 Matches
length(which(tmp.d$MaskID %!in% tmp.o$MaskID )) #0
length(which(tmp.d$MaskID %in% tmp.o$MaskID )) #389
dim(tmp.d) ## 389 x 48
dim(tmp.o) ## 389 x 24

same = colnames(d)[which(colnames(d) %in% colnames(o.data))]
same ## MaskID, arm, edu. smoking
## Rename columns from MIND-COG to match our col names in o.data
colnames(tmp.d)[which(colnames(tmp.d) %in% 'female')] = 'gender'
colnames(tmp.d)[which(colnames(tmp.d) %in% 'gly')] = 'int_gly_arm' 
colnames(tmp.d)[which(colnames(tmp.d) %in% 'bp')] = 'int_bp_arm' 
colnames(tmp.d)[which(colnames(tmp.d) %in% 'diab_dur')] = 'yrsdiab' 
colnames(tmp.d)[which(colnames(tmp.d) %in% 'race')] = 'ethnicity' 
colnames(tmp.d)[which(colnames(tmp.d) %in% 'BMI')] = 'bmi' 
colnames(tmp.d)[which(colnames(tmp.d) %in% 'antihyp')] = 'hyptens' 
colnames(tmp.d)[which(colnames(tmp.d) %in% 'age')] = 'baseline_age' ## They rounded to integer
tmp.d[,which(colnames(tmp.d) %in% 'NonWhite')] = NULL ## They are all NonWhite

#tmp.find = colnames(d)[which(colnames(d) %!in% colnames(o.data))]
#tmp.find
#tmp.look = colnames(o.data)[which(colnames(o.data) %!in% colnames(d) )]
#tmp.look


tmp.o = tmp.o[order(tmp.o$MaskID),]
tmp.d = tmp.d[order(tmp.d$MaskID),]
head(tmp.d[,c(1,which(colnames(tmp.d) %!in% colnames(tmp.o)))])
#head(tmp.o[,c(1,which(colnames(tmp.o) %!in% colnames(tmp.d) ))])

#drop = which(colnames(tmp.o) %!in% colnames(tmp.d) )

## Need to find these columns
colnames(tmp.d[,c(which(colnames(tmp.d) %!in% colnames(tmp.o)))])
#colnames(tmp.d[,drop])
 [1] "lipid"        "cvd_events"   "sbp"          "dbp"          "hba1c"       
 [6] "glucose"      "ldl"          "hdl"          "trig"         "chol"        
[11] "screat"       "ualb"         "uacr"         "ckd_gfr"      "kidfail"     
[16] "statins"      "hypertension" "LipidLower"   "prev_mi"      "base_dsst"   
[21] "base_mmse"    "base_stroop"  "base_ravlt"   "gm"           "wm"          
[26] "ticv"         "wmlv"         "ace"          "arb"          "ace_arb"     
[31] "insu"         "lths"         "hs"           "techschool"   "colgrad"     

drop = c()

## !!!!!!!!!!!!!!!!!  "lipid"


## !!!!!!!!!!!!!!!!!  "cvd_events"
tmp = read.table('~/data/accord/data/pheno_data/cvdoutcomes.txt',skip=33,header=T,sep='\t')
tmp = read.table('~/data/accord/data/pheno_data/f01_inclusionexclusionsummary.txt',header=T,sep='\t')
tmp = read.table('~/data/accord/data/pheno_data/original_pheno_data_backup_v2.txt',header=T,sep='\t')
tmp = tmp[which(tmp$MaskID %in% tmp.d$MaskID),]
tmp = tmp[order(tmp$MaskID),]
## Get BLR
tmp = tmp[which(tmp$Visit %in% 'BLR'),]
#table(tmp$Visit)
dim(tmp)
head(tmp.d[,c(1,which(colnames(tmp.d) %!in% colnames(tmp.o)))])
head(tmp)

head(tmp.d[,c(1,which(colnames(tmp.d) %in% c('cvd_events')))])
tmp[which(tmp$MaskID %in% c(100235)),]
head(tmp[,c(1,which(colnames(tmp) %in% c('cvd_hx_baseline')))])

tail(tmp[,c(1,which(colnames(tmp) %in% c('statin','acei','a2rb','reg_insulin')))])
tail(tmp.d[,c(1,which(colnames(tmp.d) %in% c('statins','ace','arb','insu')))])


##############
# sdp dbp
tmp = read.table('~/data/accord/data/pheno_data/bloodpressure.txt',skip=4,header=T,sep='\t')
## Get IDs
tmp = tmp[which(tmp$MaskID %in% tmp.d$MaskID),]
tmp = tmp[order(tmp$MaskID),]
## Get BLR
tmp = tmp[which(tmp$Visit %in% 'BLR'),]
dim(tmp)
head(tmp.d[,c(1,which(colnames(tmp.d) %!in% colnames(tmp.o)))])
head(tmp)

drop = c(drop,'sbp','dbp')
## From now on drop these
head(tmp.d[,c(1,which(colnames(tmp.d) %!in% c(colnames(tmp.o),drop)))])


#############
## hba1c: glycosylated hemoglobin (%)
tmp = read.table('~/data/accord/data/pheno_data/hba1c.txt',skip=1,header=T,sep='\t')
## Get IDs
tmp = tmp[which(tmp$MaskID %in% tmp.d$MaskID),]
tmp = tmp[order(tmp$MaskID),]
## Get BLR
tmp = tmp[which(tmp$Visit %in% 'BLR'),]
dim(tmp)
head(tmp.d[,c(1,which(colnames(tmp.d) %!in% colnames(tmp.o)))])
head(tmp)
tail(tmp)
tail(tmp.d[,c(1,which(colnames(tmp.d) %!in% colnames(tmp.o)))])

tmp.match = data.frame(ID=tmp.d$MaskID,hb=tmp.d$hba1c,hb2=99,hbdiff=99)
#head(tmp.match)
#head(tmp)
i = 1
for (i in 1:nrow(tmp.match)){
  if (length(which(tmp$MaskID %in% tmp.match$ID[i])) > 0){
    tmp.match$hb2[i] = tmp$hba1c[which(tmp$MaskID %in% tmp.match$ID[i])]
    tmp.match$hbdiff[i] = tmp.match$hb[i]-tmp.match$hb2[i]
  }else{
    tmp.match$hb2[i] = NA
  }
}

tmp[]

head(tmp.match)

tmp = read.table('~/data/accord/data/pheno_data/hba1c.txt',skip=1,header=T,sep='\t')
tmp[which(tmp$MaskID %in% tmp.match$ID[which(tmp.match$hbdiff!=0)]),]


tmp.match[which(tmp.match$hbdiff!=0),]

dim(tmp.d)
dim(tmp.o)
head(tmp.d[,c(1,17)])
colnames(tmp.d)
head(tmp)
colnames(tmp)

which(hba1c)
as.vector(outer(as.numeric(tmp.d$hba1c),as.numeric(tmp$hba1c),"="))
#### THERE WAS SOMETHING WRONG IN THIS ONE. 100307 does not match between both files

drop = c(drop,'hba1c')
## From now on drop these
head(tmp.d[,c(1,which(colnames(tmp.d) %!in% c(colnames(tmp.o),drop)))])

#############
## glucose, screat, ualb, uacr
## fpg: Fasting plasma glucose (mg/dL)
## ualb: Urinary albumin (mg/dL)
## uacr: Urinary albumin to creatinine ratio (g/mg)
## screat: Serum creatinine (mg/dL)
tmp = read.table('~/data/accord/data/pheno_data/otherlabs.txt',skip=9,header=T,sep='\t')
tmp = tmp[which(tmp$MaskID %in% tmp.d$MaskID),]
tmp = tmp[order(tmp$MaskID),]
## Get BLR
tmp = tmp[which(tmp$Visit %in% 'BLR'),]
dim(tmp)
dim(tmp)


head(tmp.d[,c(1,which(colnames(tmp.d) %in% c('ckd_gfr')))])
head(tmp[,c(1,which(colnames(tmp) %in% c('gfr')))])

tail(tmp[,c(1,which(colnames(tmp) %in% c('total_dsc','total_mmse','stroop','ravlt')))])
tail(tmp.d[,c(1,which(colnames(tmp.d) %in% c('base_dsst','base_mmse','base_stroop','base_ravlt')))])



colnames(tmp.d)[which(colnames(tmp.d) %in% 'glucose')] = 'fpg' 
drop = c(drop,'fpg','ualb','uacr','screat')
head(tmp.d[,c(1,which(colnames(tmp.d) %!in% c(colnames(tmp.o),drop)))])


#############
## ldl, hdl, trig, chol
##
## chol: Total Cholesterol (mg/dL)
## trig: Triglycerides (mg/dL)
## vldl: Very low density lipoprotein (mg/dL)
## ldl: Low density lipoprotein (mg/dL)
## hdl: High density lipoprotein (mg/dL)
tmp = read.table('~/data/accord/data/pheno_data/lipids.txt',skip=5,header=T,sep='\t')
tmp = tmp[which(tmp$MaskID %in% tmp.d$MaskID),]
tmp = tmp[order(tmp$MaskID),]
## Get BLR
tmp = tmp[which(tmp$Visit %in% 'BLR'),]
dim(tmp)

head(tmp.d[,c(1,which(colnames(tmp.d) %!in% c(colnames(tmp.o),drop)))])
head(tmp)

tail(tmp)
tail(tmp.d[,c(1,which(colnames(tmp.d) %!in% c(colnames(tmp.o),drop)))])

#colnames(tmp.d)[which(colnames(tmp.d) %in% 'glucose')] = 'fpg' 
drop = c(drop,'ldl','hdl','trig','chol')
head(tmp.d[,c(1,which(colnames(tmp.d) %!in% c(colnames(tmp.o),drop)))])


#############
## statins
## statin: HMG CoA reducatce inhibitors (statins)                                                                                                   $
tmp = read.table('~/data/accord/data/pheno_data/concomitantmeds.txt',skip=57,header=T,sep='\t')
tmp = tmp[which(tmp$MaskID %in% tmp.d$MaskID),]
tmp = tmp[order(tmp$MaskID),]
## Get BLR
tmp = tmp[which(tmp$Visit %in% 'BLR'),]
dim(tmp)

head(tmp.d[,c(1,which(colnames(tmp.d) %in% c('statins','ace','arb','insu')))])
head(tmp[,c(1,which(colnames(tmp) %in% c('statin','acei','a2rb','reg_insulin')))])

tail(tmp[,c(1,which(colnames(tmp) %in% c('statin','acei','a2rb','reg_insulin')))])
tail(tmp.d[,c(1,which(colnames(tmp.d) %in% c('statins','ace','arb','insu')))])

colnames(tmp.d)[which(colnames(tmp.d) %in% 'statins')] = 'statin' 
colnames(tmp.d)[which(colnames(tmp.d) %in% 'ace')] = 'ace' 
colnames(tmp.d)[which(colnames(tmp.d) %in% 'arb')] = 'a2rb' 
drop = c(drop,'statin','ace','a2rb')
head(tmp.d[,c(1,which(colnames(tmp.d) %!in% c(colnames(tmp.o),drop)))])

##!!!!!!!!!!!!!!!!!!! ace_arb is probably ace + arb
##!!!!!!!!!!!!!!!!!!! insu does not match reg_insulin

##drop = unique(drop)
#############
## base_dsst, base_stroop, base_mmse, base_ravlt
# total_mmse: MMSE Total Score
# total_dsc: DSST Total Score
# stroop: STROOP Interference Score
# ravlt: RAVLT Score

tmp = read.table('~/data/accord/data/pheno_data/mind.txt',skip=4,header=T,sep='\t')
tmp = tmp[which(tmp$MaskID %in% tmp.d$MaskID),]
tmp = tmp[order(tmp$MaskID),]
## Get BLR
tmp = tmp[which(tmp$Visit %in% 'BLR'),]
dim(tmp)

head(tmp.d[,c(1,which(colnames(tmp.d) %in% c('base_dsst','base_mmse','base_stroop','base_ravlt')))])
head(tmp[,c(1,which(colnames(tmp) %in% c('total_dsc','total_mmse','stroop','ravlt')))])

tail(tmp[,c(1,which(colnames(tmp) %in% c('total_dsc','total_mmse','stroop','ravlt')))])
tail(tmp.d[,c(1,which(colnames(tmp.d) %in% c('base_dsst','base_mmse','base_stroop','base_ravlt')))])

colnames(tmp.d)[which(colnames(tmp.d) %in% 'base_dsst')] = 'total_dsc' 
colnames(tmp.d)[which(colnames(tmp.d) %in% 'base_mmse')] = 'total_mmse' 
colnames(tmp.d)[which(colnames(tmp.d) %in% 'base_stroop')] = 'stroop' 
colnames(tmp.d)[which(colnames(tmp.d) %in% 'base_ravlt')] = 'ravlt' 
drop = c(drop,'total_dsc','total_mmse','stroop','ravlt')
head(tmp.d[,c(1,which(colnames(tmp.d) %!in% c(colnames(tmp.o),drop)))])



#############
## ticv, wm, gm, wmlv, 
# Brain_Volume_total_sum: Brain Volume, Total (cc, summed)
# White_Matter_total_sum: White Matter, Total (cc, summed)
# Gray_Matter_total_sum: Gray Matter, Total (cc, summed)
# White_Matter_abnormal_sum: White Matter, Abnormal (cc, summed)
tmp = read.table('~/data/accord/data/pheno_data/mind_mri.txt',skip=4,header=T,sep='\t')
tmp = tmp[which(tmp$MaskID %in% tmp.d$MaskID),]
tmp = tmp[order(tmp$MaskID),]
## Get BLR
tmp = tmp[which(tmp$Visit %in% 'BLR'),]
dim(tmp)

head(tmp.d[-which(is.na(tmp.d$ticv)),c(1,which(colnames(tmp.d) %in% c('ticv','wm','gm','wmlv')))])
head(tmp[,c(1,which(colnames(tmp) %in% c('Brain_Volume_total_sum','White_Matter_total_sum','Gray_Matter_total_sum','White_Matter_abnormal_sum')))])


tail(tmp[,c(1,which(colnames(tmp) %in% c('Brain_Volume_total_sum','White_Matter_total_sum','Gray_Matter_total_sum','White_Matter_abnormal_sum')))])
tail(tmp.d[-which(is.na(tmp.d$ticv)),c(1,which(colnames(tmp.d) %in% c('ticv','wm','gm','wmlv')))])

colnames(tmp.d)[which(colnames(tmp.d) %in% 'ticv')] = 'Brain_Volume_total_sum' 
colnames(tmp.d)[which(colnames(tmp.d) %in% 'wm')] = 'White_Matter_total_sum' 
colnames(tmp.d)[which(colnames(tmp.d) %in% 'gm')] = 'Gray_Matter_total_sum' 
colnames(tmp.d)[which(colnames(tmp.d) %in% 'wmlv')] = 'White_Matter_abnormal_sum' 
drop = c(drop,'Brain_Volume_total_sum','White_Matter_total_sum','Gray_Matter_total_sum','White_Matter_abnormal_sum')
head(tmp.d[,c(1,which(colnames(tmp.d) %!in% c(colnames(tmp.o),drop)))])

colnames(tmp.d[,c(1,which(colnames(tmp.d) %!in% c(colnames(tmp.o),drop)))])

# !!!!!!!!!!!!!!! ALL rounded down: gm = 1/100th, wm = 1/100th, ticv = 1/100th, wmlv = 1 


### lts hs techschool colgrad
## 1, 2, 3, 4
tmp = read.table('~/data/accord/data/pheno_data/f07_baselinehistoryphysicalexam.txt',header=T,sep='\t')
tmp = tmp[which(tmp$MaskID %in% tmp.d$MaskID),]
tmp = tmp[order(tmp$MaskID),]
## Get BLR
tmp = tmp[which(tmp$Visit %in% 'BLR'),]
dim(tmp)
head(tmp)

head(tmp.d)

head(tmp.d[-which(is.na(tmp.d$ticv)),c(1,which(colnames(tmp.d) %in% c('ticv','wm','gm','wmlv')))])
head(tmp[,c(1,which(colnames(tmp) %in% c('edu')))])

head(tmp.d[,c(1,which(colnames(tmp.d) %!in% c(colnames(tmp.o),drop)))])

## prev_mi
# x2mi
tmp = read.table('~/data/accord/data/pheno_data/f01_inclusionexclusionsummary.txt',header=T,sep='\t')




## END MIND COG Jasmin
############################################################

### John House - Resistant Hypertension (rhtn) ###
### We started R from Terminal SSH from the home/accord/data  folder
### R script that created the cases/controls and corresponding .HTML are located in 
### pheno_data/rhtn_phenos!

library(tidyverse)
# Read in covariates file
# Female =1 in the gender column! Look at the data dictionary from Jon Leirer
cov.data <- read.table(file.path("pheno_data","original_pheno_data_backup_v2.txt"),stringsAsFactors=FALSE,sep="\t",header=TRUE)
# read in my case/control definition
read.table(file.path("pheno_data","rhtn_phenos","RHTN_Case_Control_Genotyped.11.1.2019.csv"), 
                     sep = ",", 
                     header = T, 
                    stringsAsFactors = F) %>% 
  dplyr::select(MaskID, RHTN_status) %>%  
  dplyr::filter(RHTN_status != "REMOVE") %>% 
  mutate(MaskID = as.character(MaskID)) %>% 
  mutate(RHTN = ifelse(RHTN_status == "CONTROL", 1, 
                       ifelse(RHTN_status == "CASE", 2, 0))) %>% 
  left_join(., cov.data, by = "MaskID") %>% 
  dplyr::select(MaskID, LabID, everything()) -> o.data

## remove missing forced covariates
o.data %>% dplyr::filter(!(is.na(gender) | is.na(bmi))) -> o.data

  # dplyr::filter(gfr >= 30) -> o.data
# ARMs for our dataset ### this will need to be a covariate
# 4 = Intensive Glycemia/Standard BP
# 2 = Standard Glycemia/Standard BP
# Note: Exclusionary Criteria has been applied here as well
write.table(o.data,file.path("analysis","pheno_data_rhtn.txt"),row.names=FALSE,sep="\t",quote=FALSE)
### END RHTN jsh ###




