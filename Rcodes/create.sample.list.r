
################ BEGIN ##################################
##### RUN THIS SECTION OF CODE EVERY TIME ###############
rm(list=ls())
#setwd("Y://NCSU_Projects//Accord//Data//analysis")
#setwd("/Volumes/Kirk 1/NCSU_Projects/Accord/Data/analysis/")
#setwd('~/data/accord/data/analysis')
#setwd('~/NCSU_Projects/Accord/Data/analysis')
#setwd('~/NCSU_Projects/Accord/accord_data/analysis')
setwd('/home/accord/data/analysis')
data<-read.table("pheno_data_hba1c.txt",header=T,stringsAsFactors=F,sep='\t')
# data <- read.table("pheno_data_tzd_weight_gain.txt",header=T,stringsAsFactors=F,sep='\t')
#lipidmeds<-read.table("../pheno_data/f16_lipidmedicationsmanagement.txt",sep="\t",header=T,stringsAsFactors=F)
################## END ##################################
####### RUN THIS SECTION OF CODE EVERY TIME #############


# x<-data[data$arm%in%c(5,6,7,8) & data$ethnicity%in%"White",]
# f<-"cholesterol/blr_lipid_white/starting_samples.txt"
# write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)

# x<-data[data$arm%in%c(5,6,7,8),]
# f<-"cholesterol/d4m_lipid/starting_samples.txt"
# write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)

# x<-data[data$arm%in%c(5,6,7,8) & data$ethnicity%in%"White",]
# f<-"cholesterol/d4m_lipid_white/starting_samples.txt"
# write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)

# x<-data[data$arm%in%c(5,6,7,8) & data$ethnicity%in%"Black",]
# f<-"cholesterol/d4m_lipid_black/starting_samples.txt"
# write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)

# x<-data[data$arm%in%c(5,6,7,8) & data$statin%in%0 & data$f04_zadhere%in%1,]
# f<-"cholesterol/d4m_lipid_arm_no_statin_blr/starting_samples.txt"
# write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)

#x<-data[data$arm%in%c(5,7) & data$fibrate%in%0 & data$f04_fadhere%in%1,]
#f<-"cholesterol/d4m_lipid_arm_no_fibrate_blr/starting_samples.txt"
#write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)

#pull white subjects for change in ldl, hdl, chol, trig after fibrate trtmnt
x<-data[data$arm%in%c(5,7) & data$fibrate%in%0 & data$f04_fadhere%in%1,]
x <- x[which(x[,"ethnicity"]=="Black"),]
f<-"Fibrate/GOLDN_black_log_d4m_lipid_arm_no_fibrate_blr/starting_samples.txt"
write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)


#create starting sample list for Metformin based on 18 month scores and hba1c_time_adjusted vals
x<-data[data$metforminscore_18mo==3 & !is.na(data$prt_adj_hba1c) & !is.na(data$min_hba1c_timeadj),]
f.list<-c("Metformin_MetGen/metformin18mo_glymed3_sulfonyl34_update/starting_samples.txt",
     "Metformin_MetGen/metformin18mo_glymed3_sulfonyl34_update_forcetime/starting_samples.txt",
     "Metformin_MetGen/metformin18mo_minhba1cCov_glymed3_sulfonyl34_update/starting_samples.txt",
     "Metformin_MetGen/metformin18mo_minhba1cCov_glymed3_sulfonyl34_update_forcetime/starting_samples.txt")
for(f in f.list){
  write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)
}


#create starting sample list for Metformin based on 18 month scores and hba1c_time_adjusted vals only STD ARM!!
x<-data[data$arm%in%c(1,2,5,6) & data$metforminscore_18mo==3 & !is.na(data$prt_adj_hba1c) & !is.na(data$min_hba1c_timeadj),]
f.list<-c("Metformin_MetGen/metformin18mo_glymed3_sulfonyl34_forcetime_stdArm/starting_samples.txt",
          "Metformin_MetGen/metformin18mo_minhba1cCov_glymed3_sulfonyl34_forcetime_stdArm/starting_samples.txt")
for(f in f.list){
  write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)
}

#create starting sample list for Metformin based on 6-18 month scores and hba1c_time_adjusted vals
x<-data[data$metforminscore_6to18mo==3 & !is.na(data$prt_adj_hba1c) & !is.na(data$min_hba1c_timeadj),]
f.list<-c("Metformin_MetGen/metformin8to16mo_minhba1c_glymed3/starting_samples.txt",
          "Metformin_MetGen/metformin8to16mo_minhba1c_glymed4/starting_samples.txt",
          "Metformin_MetGen/metformin8to16mo_minhba1c_glymed34/starting_samples.txt")
for(f in f.list){
  write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)
}


#create starting sample list for Metformin based on 18 month scores and hba1c_time_adjusted vals and must be white--Std Arm
x<-data[data$arm%in%c(1,2,5,6) & data$metforminscore_18mo==3 & !is.na(data$min_hba1c_timeadj),]
x <- x[which(x[,"ethnicity"]=="White"),]
f <- "Metformin_MetGen/metformin18mo_minhba1cCov_stdArm_10list/starting_samples.txt"
write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)
#create starting sample list for Metformin based on 18 month scores and hba1c_time_adjusted vals and must be white--ALL arms
x<-data[data$metforminscore_18mo==3 & !is.na(data$min_hba1c_timeadj),]
x <- x[which(x[,"ethnicity"]=="White"),]
f <- "Metformin_MetGen/metformin18mo_minhba1cCov_ALLArm_10list/starting_samples.txt"
write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)



#create starting sample list for Metformin based on 18 month scores and hba1c_time_adjusted vals no other drugs!!
tzd34 <- read.table(file.path("../pheno_data/tzdscore_3and4.txt"),header=TRUE,sep="\t") 
sulf34 <-read.table(file.path("../pheno_data/sulfonylureascore_3and4.txt"),header=TRUE,sep="\t") 
meg <- read.table(file.path("../pheno_data/meglitinidescore_3and4.txt"),header=TRUE,sep="\t") 
ins <-  read.table(file.path("../pheno_data/g2avetid_time_adjusted_score.txt"),header=TRUE,sep="\t") 
drug.data <- cbind(tzd34,sulf34[,2],meg[,2],ins[,4])  
drug.data[is.na(drug.data[,5]),5] <- 0
d.free.samps <- drug.data[which(rowSums(drug.data[,2:ncol(drug.data)])==0),1]

x <-data[data$metforminscore_18mo==3 & !is.na(data$prt_adj_hba1c) & !is.na(data$min_hba1c_timeadj),]
x.samps <- intersect(x[,"MaskID"],d.free.samps)
x <-x[which(x[,"MaskID"]%in%x.samps),]
f.list<-c("Metformin_MetGen/metformin18mo_glymed3_noOtherDrug/starting_samples.txt",
          "Metformin_MetGen/metformin18mo_minhba1cCov_glymed3_noOtherDrug/starting_samples.txt")
for(f in f.list){
  write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)
}


#split above into whites and non-whites
x <-data[data$metforminscore_18mo==3 & !is.na(data$prt_adj_hba1c) & !is.na(data$min_hba1c_timeadj) & data$ethnicity=="White",]
x.samps <- intersect(x[,"MaskID"],d.free.samps)
x <-x[which(x[,"MaskID"]%in%x.samps),]
f.list<-c("Metformin_MetGen/metformin18mo_glymed3_noOtherDrug_whiteOnly/starting_samples.txt",
          "Metformin_MetGen/metformin18mo_minhba1cCov_glymed3_noOtherDrug_whiteOnly/starting_samples.txt")
for(f in f.list){
  write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)
}

x <-data[data$metforminscore_18mo==3 & !is.na(data$prt_adj_hba1c) & !is.na(data$min_hba1c_timeadj) & data$ethnicity!="White",]
x.samps <- intersect(x[,"MaskID"],d.free.samps)
x <-x[which(x[,"MaskID"]%in%x.samps),]
f.list<-c("Metformin_MetGen/metformin18mo_glymed3_noOtherDrug_noWhite/starting_samples.txt",
          "Metformin_MetGen/metformin18mo_minhba1cCov_glymed3_noOtherDrug_noWhite/starting_samples.txt")
for(f in f.list){
  write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)
}

### Now for no insulin, 0 or 4 for other drugs, and glymed 3

tzd4 <- read.table(file.path("../pheno_data/tzdscore_4.txt"),header=TRUE,sep="\t") 
sulf4 <-read.table(file.path("../pheno_data/sulfonylureascore_4.txt"),header=TRUE,sep="\t") 
meg4 <- read.table(file.path("../pheno_data/meglitinidescore_4.txt"),header=TRUE,sep="\t") 
ins <-  read.table(file.path("../pheno_data/g2avetid_time_adjusted_score.txt"),header=TRUE,sep="\t") 
drug.data4 <- cbind(tzd4,sulf4[,2],meg4[,2],ins[,4])  
drug.data4[is.na(drug.data4[,5]),5] <- 0

tzd0 <- read.table(file.path("../pheno_data/tzdscore_0.txt"),header=TRUE,sep="\t") 
sulf0 <-read.table(file.path("../pheno_data/sulfonylureascore_0.txt"),header=TRUE,sep="\t") 
meg0 <- read.table(file.path("../pheno_data/meglitinidescore_0.txt"),header=TRUE,sep="\t") 
drug.data0 <- cbind(tzd0,sulf0[,2],meg0[,2])  

samps.0 <- drug.data0[which(rowSums(drug.data0[,2:ncol(drug.data0)])==3),1]

drug.data4 <- drug.data[which(drug.data4[,5]==0),]
samps.4 <- drug.data4[which(rowSums(drug.data4[,2:4])>0),1]
d.samps <- c(samps.0,samps.4)

x <-data[data$metforminscore_18mo==3 & !is.na(data$prt_adj_hba1c) & !is.na(data$min_hba1c_timeadj),]
x.samps <- intersect(x[,"MaskID"],d.samps)
x <-x[which(x[,"MaskID"]%in%x.samps),]
f.list<-c("Metformin_MetGen/metformin18mo_glymed3_stableOtherDrug/starting_samples.txt",
          "Metformin_MetGen/metformin18mo_minhba1cCov_glymed3_stableOtherDrug/starting_samples.txt")
for(f in f.list){
  write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)
}


#split above into whites and non-whites
x <-data[data$metforminscore_18mo==3 & !is.na(data$prt_adj_hba1c) & !is.na(data$min_hba1c_timeadj) & data$ethnicity=="White",]
x.samps <- intersect(x[,"MaskID"],d.samps)
x <-x[which(x[,"MaskID"]%in%x.samps),]
f.list<-c("Metformin_MetGen/metformin18mo_glymed3_stableOtherDrug_whiteOnly/starting_samples.txt",
          "Metformin_MetGen/metformin18mo_minhba1cCov_glymed3_stableOtherDrug_whiteOnly/starting_samples.txt")
for(f in f.list){
  write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)
}

x <-data[data$metforminscore_18mo==3 & !is.na(data$prt_adj_hba1c) & !is.na(data$min_hba1c_timeadj) & data$ethnicity!="White",]
x.samps <- intersect(x[,"MaskID"],d.samps)
x <-x[which(x[,"MaskID"]%in%x.samps),]
f.list<-c("Metformin_MetGen/metformin18mo_glymed3_stableOtherDrug_noWhite/starting_samples.txt",
          "Metformin_MetGen/metformin18mo_minhba1cCov_glymed3_stableOtherDrug_noWhite/starting_samples.txt")
for(f in f.list){
  write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)
}




### Now for 0 or 4 for other drugs, insulin OK, and glymed 3

tzd4 <- read.table(file.path("../pheno_data/tzdscore_4.txt"),header=TRUE,sep="\t") 
sulf4 <-read.table(file.path("../pheno_data/sulfonylureascore_4.txt"),header=TRUE,sep="\t") 
meg4 <- read.table(file.path("../pheno_data/meglitinidescore_4.txt"),header=TRUE,sep="\t") 
#ins <-  read.table(file.path("../pheno_data/g2avetid_time_adjusted_score.txt"),header=TRUE,sep="\t") 
drug.data4 <- cbind(tzd4,sulf4[,2],meg4[,2])  
#drug.data4[is.na(drug.data4[,5]),5] <- 0

tzd0 <- read.table(file.path("../pheno_data/tzdscore_0.txt"),header=TRUE,sep="\t") 
sulf0 <-read.table(file.path("../pheno_data/sulfonylureascore_0.txt"),header=TRUE,sep="\t") 
meg0 <- read.table(file.path("../pheno_data/meglitinidescore_0.txt"),header=TRUE,sep="\t") 
drug.data0 <- cbind(tzd0,sulf0[,2],meg0[,2])  

samps.0 <- drug.data0[which(rowSums(drug.data0[,2:ncol(drug.data0)])==3),1]

#drug.data4 <- drug.data[which(drug.data4[,5]==0),]
samps.4 <- drug.data4[which(rowSums(drug.data4[,2:4])>0),1]
d.samps <- c(samps.0,samps.4)

x <-data[data$metforminscore_18mo==3 & !is.na(data$prt_adj_hba1c) & !is.na(data$min_hba1c_timeadj),]
x.samps <- intersect(x[,"MaskID"],d.samps)
x <-x[which(x[,"MaskID"]%in%x.samps),]
f.list<-c("Metformin/metformin18mo_glymed4/starting_samples.txt",
          "Metformin/metformin18mo_minhba1cCov_glymed4/starting_samples.txt")
for(f in f.list){
  write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)
}


#split above into whites and non-whites
x <-data[data$metforminscore_18mo==3 & !is.na(data$prt_adj_hba1c) & !is.na(data$min_hba1c_timeadj) & data$ethnicity=="White",]
x.samps <- intersect(x[,"MaskID"],d.samps)
x <-x[which(x[,"MaskID"]%in%x.samps),]
f.list<-c("Metformin/metformin18mo_glymed4_whiteOnly/starting_samples.txt",
          "Metformin/metformin18mo_minhba1cCov_glymed4_whiteOnly/starting_samples.txt")
for(f in f.list){
  write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)
}

x <-data[data$metforminscore_18mo==3 & !is.na(data$prt_adj_hba1c) & !is.na(data$min_hba1c_timeadj) & data$ethnicity!="White",]
x.samps <- intersect(x[,"MaskID"],d.samps)
x <-x[which(x[,"MaskID"]%in%x.samps),]
f.list<-c("Metformin/metformin18mo_glymed4_noWhite/starting_samples.txt",
          "Metformin/metformin18mo_minhba1cCov_glymed4_noWhite/starting_samples.txt")
for(f in f.list){
  write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)
}



##################################
#pull in fibrate samples-using 4month fibrate score

fib.score <- read.table(file.path("../pheno_data/fibrate_phenos/fibrate_score-4mo.txt"),header=TRUE,sep="\t") 
#pull all subjects for change in ldl, hdl, chol, trig after fibrate trtmnt
d.samps <- fib.score[which(fib.score[,"Score"]==3),"MaskID"]
x<-data[data$arm%in%c(5,7),]
x.samps <- intersect(x[,"MaskID"],d.samps)
x <-x[which(x[,"MaskID"]%in%x.samps),]
f<-"fibrate/d4m_fibrate_wdrugScores/starting_samples.txt"
write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)


#pull white subjects for change in ldl, hdl, chol, trig after fibrate trtmnt
d.samps <- fib.score[which(fib.score[,"Score"]==3),"MaskID"]
x<-data[data$arm%in%c(5,7),]
x.samps <- intersect(x[,"MaskID"],d.samps)
x <-x[which(x[,"MaskID"]%in%x.samps),]
x <- x[which(x[,"ethnicity"]=="White"),]
f<-"fibrate/d4m_fibrate_wdrugScores_white/starting_samples.txt"
write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)


#pull black subjects for change in ldl, hdl, chol, trig after fibrate trtmnt
d.samps <- fib.score[which(fib.score[,"Score"]==3),"MaskID"]
x<-data[data$arm%in%c(5,7),]
x.samps <- intersect(x[,"MaskID"],d.samps)
x <-x[which(x[,"MaskID"]%in%x.samps),]
x <- x[which(x[,"ethnicity"]=="Black"),]
f<-"fibrate/d4m_fibrate_wdrugScores_black/starting_samples.txt"
write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)


######################################################################
# Nov 2, 2014

# create list of fibrate subjects 3-4mo trt window and 0/4 for other drugs
x <-data[data$fibrateScore_3==1 & 
           !is.na(data$ldl_lr) &
           !is.na(data$hdl_lr) & 
           !is.na(data$trig_lr) & 
           !is.na(data$chol_lr) & 
           (data$fib_tzdScore_0==1 | data$fib_tzdScore_4==1) & 
           (data$fib_statinScore_0==1 | data$fib_statinScore_34==1) &
           data$arm%in%c(5,7)
         ,]
f<-"Fibrate/fibrate3-4mo_lipids_med4/starting_samples.txt"
write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)

# create list of fibrate subjects 3-4mo trt window and 0 for other drugs
x <-data[data$fibrateScore_3==1 & 
           !is.na(data$ldl_lr) &
           !is.na(data$hdl_lr) & 
           !is.na(data$trig_lr) & 
           !is.na(data$chol_lr) & 
           data$fib_tzdScore_0==1 & 
           (data$fib_statinScore_0==1 | data$fib_statinScore_4==1) &
           data$fib_other_lipidmedScore==0 &
           data$fib_niacinScore==0 &
           data$fib_cholest_abiScore==0 &
           data$fib_bile_sequestrantScore==0 &
           data$arm%in%c(5,7)
         ,]
f<-"Fibrate/fibrate3-4mo_lipids_med0/starting_samples.txt"
write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)


# create list of placebo subjects 3-4mo trt window and 0/4 for other drugs
x <-data[data$fibrateScore_3==1 & 
           !is.na(data$ldl_lr) &
           !is.na(data$hdl_lr) & 
           !is.na(data$trig_lr) & 
           !is.na(data$chol_lr) & 
           (data$fib_tzdScore_0==1 | data$fib_tzdScore_4==1) & 
           (data$fib_statinScore_0==1 | data$fib_statinScore_34==1) &
           data$arm%in%c(8,6)
         ,]
f<-"Fibrate/placebo3-4mo_lipids_med4/starting_samples.txt"
write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)

# create list of placebo subjects 3-4mo trt window and 0 for other drugs
x <-data[data$fibrateScore_3==1 & 
           !is.na(data$ldl_lr) &
           !is.na(data$hdl_lr) & 
           !is.na(data$trig_lr) & 
           !is.na(data$chol_lr) & 
           data$fib_tzdScore_0==1 & 
           (data$fib_statinScore_0==1 | data$fib_statinScore_4==1) &
           data$fib_other_lipidmedScore==0 &
           data$fib_niacinScore==0 &
           data$fib_cholest_abiScore==0 &
           data$fib_bile_sequestrantScore==0 &
           data$arm%in%c(8,6)
         ,]
f<-"Fibrate/placebo3-4mo_lipids_med0/starting_samples.txt"
write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)

# create list of metformin subjects 3-9mo trt window and 0/4 for other drugs
x <-data[data$metforminScore_3==1 & 
           !is.na(data$prt_adj_hba1c) &
           (data$met_ACE_inhibitorsScore_0==1 | data$met_ACE_inhibitorsScore_4==1) & 
           (data$met_carvedilolScore_0==1 | data$met_carvedilolScore_4==1) & 
           (data$met_terazosinScore_0==1 | data$met_terazosinScore_4==1) & 
           (data$met_acarboseScore_0==1 | data$met_acarboseScore_4==1) & 
           (data$met_pramlintideScore_0==1 | data$met_pramlintideScore_4==1) & 
           (data$met_ARBScore_0==1 | data$met_ARBScore_4==1) & 
           (data$met_beta_blockerScore_0==1 | data$met_beta_blockerScore_4==1) & 
           (data$met_CCBScore_0==1 | data$met_CCBScore_4==1) & 
           (data$met_sitagliptinScore_0==1 | data$met_sitagliptinScore_4==1) & 
           (data$met_exanatideScore_0==1 | data$met_exanatideScore_4==1) &
           (data$met_triamtereneScore_0==1 | data$met_triamtereneScore_4==1) &
           (data$met_furosemideScore_0==1 | data$met_furosemideScore_4==1) &
           (data$met_meglitinideScore_0==1 | data$met_meglitinideScore_4==1) &
           (data$met_reserpineScore_0==1 | data$met_reserpineScore_4==1) &
           (data$met_statinScore_0==1 | data$met_statinScore_4==1) &
           (data$met_sulfonylureaScore_0==1 | data$met_sulfonylureaScore_4==1) &
           (data$met_thiazide_diureticScore_0==1 | data$met_thiazide_diureticScore_4==1) &
           (data$met_tzdScore_0==1 | data$met_tzdScore_4==1) &
           (data$met_hydralazineScore_0==1 | data$met_hydralazineScore_4==1)
           ,]
f<-"Metformin/metformin3-9mo_hba1c_med4/starting_samples.txt"
write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)

# create list of metformin subjects 3-9mo trt window and 0 for other drugs
x <-data[data$metforminScore_3==1 & 
           !is.na(data$prt_adj_hba1c) &
           data$met_ACE_inhibitorsScore_0==1 & 
           data$met_carvedilolScore_0==1 & 
           data$met_terazosinScore_0==1 & 
           data$met_acarboseScore_0==1 & 
           data$met_pramlintideScore_0==1 & 
           data$met_ARBScore_0==1 & 
           data$met_beta_blockerScore_0==1 & 
           data$met_CCBScore_0==1 & 
           data$met_sitagliptinScore_0==1 & 
           data$met_exanatideScore_0==1 &
           data$met_triamtereneScore_0==1 &
           data$met_furosemideScore_0==1 &
           data$met_meglitinideScore_0==1 &
           data$met_reserpineScore_0==1 &
           data$met_statinScore_0==1 &
           data$met_sulfonylureaScore_0==1 &
           data$met_thiazide_diureticScore_0==1 &
           data$met_tzdScore_0==1 &
           data$met_hydralazineScore_0==1
         ,]
f<-"Metformin/metformin3-9mo_hba1c_med0/starting_samples.txt"
write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)




# Statin NA, look at f16_lipidmedicationsmanagement.txt for MaskIDs below
#x<-data[data$arm%in%c(5,6,7,8) & data$statin%in%0,]
#na<-which(is.na(x$f04_zadhere))
#ii<-which(lipidmeds$MaskID%in%x$MaskID[na] & lipidmeds$Visit=="F04")
# > length(na)
# [1] 197
# > length(ii)
# [1] 143
# > lipidmeds$MaskID[ii[1]]
# [1] 100031
# > lipidmeds$MaskID[ii[2]]
# [1] 100080
# > lipidmeds$MaskID[ii[3]]
# [1] 100087

# table(x$f04_zadhere)
# for (i in 1:dim(x)[1]) {
#     if (is.na(x$f04_zadhere[i])) {
#         ii<-which(lipidmeds$MaskID==x$MaskID[i] & lipidmeds$Visit=="F04")
#         if (lipidmeds$zadhere[ii-1]%in%c(2,3))
#         x$f04_zadhere[i]<-1
#     }
# }
# table(x$f04_zadhere)

# Fibrate NA
x<-data[data$arm%in%c(5,7) & data$fibrate%in%0,]
na<-which(is.na(x$f04_fadhere))
ii<-which(lipidmeds$MaskID%in%x$MaskID[na] & lipidmeds$Visit=="F04")
# > length(na)
# [1] 113
# > length(ii)
# [1] 47
# > lipidmeds$MaskID[ii[1]]
# [1] 100157
# > lipidmeds$MaskID[ii[2]]
# [1] 100225
# > lipidmeds$MaskID[ii[3]]
# [1] 100389



######################################################################
#Mar 18, 2015

##################################
#pull in fibrate samples-using 90d fibrate score

#fib.score <- read.table(file.path("../pheno_data/fibrate_phenos/90d/fibrate_score-4mo.txt"),header=TRUE,sep="\t") 
#pull all subjects for change in ldl, hdl, chol, trig after fibrate trtmnt
d.samps <- data[which(data[,"fibrate_score"]==3),"MaskID"]
x<-data[data$arm%in%c(5,7),]
x.samps <- intersect(x[,"MaskID"],d.samps)
#x <-x[which(x[,"MaskID"]%in%x.samps),]
f<-"Fibrate/fibrate_90d_lipids_alldrugs_allraces3/starting_samples.txt"
write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)


#pull white subjects for change in ldl, hdl, chol, trig after fibrate trtmnt
d.samps <- fib.score[which(fib.score[,"Score"]==3),"MaskID"]
x<-data[data$arm%in%c(5,7),]
x.samps <- intersect(x[,"MaskID"],d.samps)
x <-x[which(x[,"MaskID"]%in%x.samps),]
x <- x[which(x[,"ethnicity"]=="White"),]
f<-"Fibrate/fibrate_90d_lipids_alldrugs_white3/starting_samples.txt"
write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)


#pull black subjects for change in ldl, hdl, chol, trig after fibrate trtmnt
d.samps <- fib.score[which(fib.score[,"Score"]==3),"MaskID"]
x<-data[data$arm%in%c(5,7),]
x.samps <- intersect(x[,"MaskID"],d.samps)
x <-x[which(x[,"MaskID"]%in%x.samps),]
x <- x[which(x[,"ethnicity"]=="Black"),]
f<-"Fibrate/fibrate_90d_lipids_alldrugs_black3/starting_samples.txt"
write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)

#pull hispanic subjects for change in ldl, hdl, chol, trig after fibrate trtmnt
d.samps <- data[which(data[,"fibrate_score"]==3),"MaskID"]#d.samps <- fib.score[which(fib.score[,"Score"]==3),"MaskID"]
x<-data[data$arm%in%c(5,7),]
x.samps <- intersect(x[,"MaskID"],d.samps)
x <-x[which(x[,"MaskID"]%in%x.samps),]
x <- x[which(x[,"ethnicity"]=="Hispanic"),]
f<-"Fibrate/fibrate_90d_lipids_alldrugs_hispanic3/starting_samples.txt"
write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)



##placebo

#pull all subjects for change in ldl, hdl, chol, trig after fibrate trtmnt
d.samps <- data[which(data[,"fibrate_score"]==3),"MaskID"]
x<-data[data$arm%in%c(6,8),]
x.samps <- intersect(x[,"MaskID"],d.samps)
#x <-x[which(x[,"MaskID"]%in%x.samps),]
f<-"Fibrate/placebo_90d_lipids_alldrugs_allraces3/starting_samples.txt"
write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)


#pull white subjects for change in ldl, hdl, chol, trig after fibrate trtmnt
d.samps <- fib.score[which(fib.score[,"Score"]==3),"MaskID"]
x<-data[data$arm%in%c(6,8),]
x.samps <- intersect(x[,"MaskID"],d.samps)
x <-x[which(x[,"MaskID"]%in%x.samps),]
x <- x[which(x[,"ethnicity"]=="White"),]
f<-"Fibrate/placebo_90d_lipids_alldrugs_white3/starting_samples.txt"
write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)


#pull black subjects for change in ldl, hdl, chol, trig after fibrate trtmnt
d.samps <- fib.score[which(fib.score[,"Score"]==3),"MaskID"]
x<-data[data$arm%in%c(6,8),]
x.samps <- intersect(x[,"MaskID"],d.samps)
x <-x[which(x[,"MaskID"]%in%x.samps),]
x <- x[which(x[,"ethnicity"]=="Black"),]
f<-"Fibrate/placebo_90d_lipids_alldrugs_black3/starting_samples.txt"
write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)


#pull hispanic subjects for change in ldl, hdl, chol, trig after fibrate trtmnt
d.samps <- data[which(data[,"fibrate_score"]==3),"MaskID"]#d.samps <- fib.score[which(fib.score[,"Score"]==3),"MaskID"]
x<-data[data$arm%in%c(6,8),]
x.samps <- intersect(x[,"MaskID"],d.samps)
x <-x[which(x[,"MaskID"]%in%x.samps),]
x <- x[which(x[,"ethnicity"]=="Hispanic"),]
f<-"Fibrate/placebo_90d_lipids_alldrugs_hispanic3/starting_samples.txt"
write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)







### Feb 16, 2015
# pull samples for fibrate-tzd interaction analysis

d.samps <- data[which(data[,"fibrate_score"]==3),"MaskID"]
x<-data[data$arm%in%c(5,7),]
x.samps <- intersect(x[,"MaskID"],d.samps)
x <-x[which(x[,"MaskID"]%in%x.samps),]
f<-"starting_samples.txt"
write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)





######################################################################
#July 13, 2015

##################################
#pull in TZD samples-using 90-270d TZD score

fib.score <- read.csv(file.path("../pheno_data/tzd_phenos/90d/tzdscores.csv"),header=TRUE) 
#pull all subjects for change in ldl, hdl, chol, trig after fibrate trtmnt
d.samps <- data[which(data[,"tzd_score"]==3),"MaskID"]
#x<-data[data$arm%in%c(5,7),]
x.samps <- intersect(data[,"MaskID"],d.samps)
x <-data[which(data[,"MaskID"]%in%x.samps),]
f<-"TZD/tzd_90d_allraces/starting_samples.txt"
write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)


#pull white subjects for change in ldl, hdl, chol, trig after fibrate trtmnt
d.samps <- fib.score[which(fib.score[,"Score"]==3),"MaskID"]
#x<-data[data$arm%in%c(5,7),]
x.samps <- intersect(data[,"MaskID"],d.samps)
x <-x[which(x[,"MaskID"]%in%x.samps),]
x <- x[which(x[,"ethnicity"]=="White"),]
f<-"TZD/tzd_90d_white/starting_samples.txt"
write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)


#pull black subjects for change in ldl, hdl, chol, trig after fibrate trtmnt
d.samps <- fib.score[which(fib.score[,"Score"]==3),"MaskID"]
x<-data#[data$arm%in%c(5,7),]
x.samps <- intersect(x[,"MaskID"],d.samps)
x <-x[which(x[,"MaskID"]%in%x.samps),]
x <- x[which(x[,"ethnicity"]=="Black"),]
f<-"TZD/tzd_90d_black/starting_samples.txt"
write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)

#pull hispanic subjects for change in ldl, hdl, chol, trig after fibrate trtmnt
d.samps <- fib.score[which(fib.score[,"Score"]==3),"MaskID"]
x<-data#x<-data[data$arm%in%c(5,7),]
x.samps <- intersect(x[,"MaskID"],d.samps)
x <-x[which(x[,"MaskID"]%in%x.samps),]
x <- x[which(x[,"ethnicity"]=="Hispanic"),]
f<-"TZD/tzd_90d_hispanic/starting_samples.txt"
write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)


######################################################################
#Aug 19, 2015

##################################
#pull in Metformin samples-using 90-270d Metformin score

met.score <- read.csv(file.path("../pheno_data/metformin_phenos/90d/metforminscores.csv"),header=TRUE) 
#pull all subjects for change in hba1c after metformin trtmnt
d.samps <- data[which(data[,"metformin_score"]==3),"MaskID"]
#x<-data[data$arm%in%c(5,7),]
x.samps <- intersect(data[,"MaskID"],d.samps)
x <-data[which(data[,"MaskID"]%in%x.samps),]
f<-"Metformin/met_90d_allraces/starting_samples.txt"
write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)


#pull white subjects for change in hba1c after metformin trtmnt
d.samps <- met.score[which(met.score[,"Score"]==3),"MaskID"]
#x<-data[data$arm%in%c(5,7),]
x.samps <- intersect(data[,"MaskID"],d.samps)
x <-x[which(x[,"MaskID"]%in%x.samps),]
x <- x[which(x[,"ethnicity"]=="White"),]
f<-"Metformin/met_90d_white/starting_samples.txt"
write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)


#pull black subjects for change in hba1c after metformin trtmnt
d.samps <- met.score[which(met.score[,"Score"]==3),"MaskID"]
x<-data#[data$arm%in%c(5,7),]
x.samps <- intersect(x[,"MaskID"],d.samps)
x <-x[which(x[,"MaskID"]%in%x.samps),]
x <- x[which(x[,"ethnicity"]=="Black"),]
f<-"Metformin/met_90d_black/starting_samples.txt"
write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)

#pull hispanic subjects for change in hba1c after metformin trtmnt
d.samps <- met.score[which(met.score[,"Score"]==3),"MaskID"]
x<-data#x<-data[data$arm%in%c(5,7),]
x.samps <- intersect(x[,"MaskID"],d.samps)
x <-x[which(x[,"MaskID"]%in%x.samps),]
x <- x[which(x[,"ethnicity"]=="Hispanic"),]
f<-"Metformin/met_90d_hispanic/starting_samples.txt"
write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)




######################################################################
#Nov 20, 2015

##################################
#pull in STATIN samples-using 120d Statin score and from placebo arm

#stat.score <- read.csv(file.path("../pheno_data/statin_phenos/120d/statinscores.txt"),header=TRUE,stringsAsFactors=FALSE,sep="\t") 
#pull all subjects for change in LDL etc after statin trtmnt
d.samps <- data[intersect(which(data[,"statin_score"]==3),which(data[,"arm"]%in%c(6,8))),"MaskID"]
#x<-data[data$arm%in%c(5,7),]
x.samps <- intersect(data[,"MaskID"],d.samps)
x <-data[which(data[,"MaskID"]%in%x.samps),]
f<-"Statin/statin_120d_lipids_alldrugs_allraces/starting_samples.txt"
write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)


#pull white subjects for change in LDL etc after statin trtmnt
d.samps <- data[which(data[,"statin_score"]==3 & data[,"arm"]%in%c(6,8) & data[,"ethnicity"]=="White"),"MaskID"]
#x<-data[data$arm%in%c(5,7),]
x.samps <- intersect(data[,"MaskID"],d.samps)
x <-x[which(x[,"MaskID"]%in%x.samps),]
#x <- x[which(x[,"ethnicity"]=="White"),]
f<-"Statin/statin_120d_lipids_alldrugs_white/starting_samples.txt"
write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)



######################################################################
#April 11, 2016

##################################
#pull in Sulfonylurea samples-using 90-270d Sulfonylurea score

sulf.score <- read.csv(file.path("../pheno_data/sulfonylurea_phenos/90d/sulfscores.csv"),header=TRUE) 
#pull all subjects for change in hba1c after metformin trtmnt
d.samps <- data[which(data[,"sulf_score"]==3),"MaskID"]
#x<-data[data$arm%in%c(5,7),]
x.samps <- intersect(data[,"MaskID"],d.samps)
x <-data[which(data[,"MaskID"]%in%x.samps),]
f<-"Sulfonylurea/sulf_90d_allraces/starting_samples.txt"
write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)


#pull white subjects for change in hba1c after metformin trtmnt
d.samps <- sulf.score[which(sulf.score[,"Score"]==3),"MaskID"]
#x<-data[data$arm%in%c(5,7),]
x.samps <- intersect(data[,"MaskID"],d.samps)
x <-x[which(x[,"MaskID"]%in%x.samps),]
x <- x[which(x[,"ethnicity"]=="White"),]
f<-"Sulfonylurea/sulf_90d_white/starting_samples.txt"
write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)


#pull black subjects for change in hba1c after metformin trtmnt
d.samps <- sulf.score[which(sulf.score[,"Score"]==3),"MaskID"]
x<-data#[data$arm%in%c(5,7),]
x.samps <- intersect(x[,"MaskID"],d.samps)
x <-x[which(x[,"MaskID"]%in%x.samps),]
x <- x[which(x[,"ethnicity"]=="Black"),]
f<-"Sulfonylurea/sulf_90d_black/starting_samples.txt"
write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)




######################################################################
# May 13 2016

##################################
#pull in baseline HGI samples
f <- "hgi_baseline/hgi_allraces/starting_samples.txt"
write.table(cbind(data$MaskID,data$LabID),f,quote=F,row.names=F,col.names=F)


#pull white subjects for change in hba1c after metformin trtmnt
#d.samps <- sulf.score[which(sulf.score[,"Score"]==3),"MaskID"]
#x<-data[data$arm%in%c(5,7),]
#x.samps <- intersect(data[,"MaskID"],d.samps)
#x <-x[which(x[,"MaskID"]%in%x.samps),]
x <- data[which(data[,"ethnicity"]=="White"),]
#head(x)
f<-"hgi_baseline/hgi_white/starting_samples.txt"
write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)

x <- data[which(data[,"ethnicity"]=="Black"),]
#head(x)
f <- "hgi_baseline/hgi_black/starting_samples.txt"
write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)

x <- data[which(data[,"ethnicity"]=="Hispanic"),]
#head(x)
f <- "hgi_baseline/hgi_hispanic/starting_samples.txt"
write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)



######################################################################
# July 5 2016

##################################
#pull in samples for statin only test-Drug Scoring Paper

stat.score <- read.csv(file.path("../pheno_data/statin_phenos/120d_drugScorePaper/statinscores.csv"),header=TRUE) 
#pull all subjects for change in hba1c after metformin trtmnt
d.samps <- data[which(data[,"statin_score"]==3),"MaskID"]
#x<-data[data$arm%in%c(5,7),]
x.samps <- intersect(data[,"MaskID"],d.samps)
x <-data[which(data[,"MaskID"]%in%x.samps),]
f<-"HillaryTestData/statin_allraces/starting_samples.txt"
write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)



######################################################################
# Aug 7 2016

##################################
#pull in samples for diet associations with dietary info
diet.score <- read.csv("../pheno_data/dietary_scores.csv") 
#pull all subjects
d.samps <- diet.score$MaskID
x.samps <- intersect(data[,"MaskID"],d.samps)
x <-data[which(data[,"MaskID"]%in%x.samps),]
f<-"diet/diet_allraces/starting_samples.txt"
write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)
# add phenos here
pheno.data <- read.table("pheno_data_diet.txt",header=TRUE,nrows=1)
pheno.list <-colnames(pheno.data)[25:ncol(pheno.data)]
f<-"diet/diet_allraces/phenotypes.txt"
write.table(pheno.list,f,quote=F,row.names=F,col.names=F)



######################################################################
# Nov 7 2016

##################################
#pull in baseline HGI samples
f <- "GV/GV_allraces/starting_samples.txt"
length(data$MaskID)
write.table(cbind(data$MaskID,data$LabID),f,quote=F,row.names=F,col.names=F)

## stratify all samples by race 
x <- data[which(data[,"ethnicity"]=="White"),]
#head(x)
f<-"GV/GV_white/starting_samples.txt"
write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)
#x <- data[which(data[,"ethnicity"]=="Black"),]
#head(x)
#f <- "GV/GV_black/starting_samples.txt"
#write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)
#x <- data[which(data[,"ethnicity"]=="Hispanic"),]
#head(x)
#f <- "GV/GV_hispanic/starting_samples.txt"
#write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)

## Stratify all samples by arm
table(data$arm)

i <- names(table(data$arm))[1]
for (i in names(table(data$arm))){
    f <- paste0("GV/GV_arm",i,"_allraces/starting_samples.txt")
    x <- data[which(data[,"arm"]==i),]
    write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)
    
    x <- data[which((data[,"ethnicity"]=="White") &(data[,"arm"]==i) ),]
    f<-paste0("GV/GV_arm",i,"_white/starting_samples.txt")
    write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)
}
    #x <- data[which(data[,"ethnicity"]=="Black"),]
    #f <- "GV/GV_black/starting_samples.txt"
    #write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)
    
    #x <- data[which(data[,"ethnicity"]=="Hispanic"),]    
    #f <- "GV/GV_hispanic/starting_samples.txt"
    #write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)
    

## Dec 6th ICAPS Study
##################################

dim(data)
## Drop all individuals not on BP subtrial
data <- data[which(!is.na(data$Other)),]
dim(data) ## 3667


colnames(data)
data$tmp <- 0
## Drop everyone based on drug exclusion criteria (<2 visits on drug--see analysis plan and implementation by Michael
i <- 30
for (i in 1:nrow(data)){
    data$tmp[i] <- sum(data[i,grep('exclusion',colnames(data))])
}

dim(data)
data <- data[-which(data$tmp>0),]
dim(data)

## Thiazide
f <- "icaps/icaps.thiazide.all_races/starting_samples.txt"
x <- data[which(data[,"Thiazide"]==1),] ## Individual is on at least Thiazide
head(x)
dim(x) ##2618
write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)

## CCB
f <- "icaps/icaps.ccb.all_races/starting_samples.txt"
x <- data[which(data[,"CCB"]==1),] ## Individual is on at least Thiazide
head(x)
dim(x) ##1653
write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)

f <- "icaps/icaps.bb.all_races/starting_samples.txt"
x <- data[which(data[,"BB"]==1),] ## Individual is on at least Thiazide
head(x)
dim(x) ##2023
write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)

f <- "icaps/icaps.aceia.arb.all_races/starting_samples.txt"
x <- data[which(data[,"ACEI.ARB"]==1),] ## Individual is on at least Thiazide
head(x)
dim(x) ##3356
write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)

'''
f <- "icaps/icaps.thiazide.white/staring_samples.txt"
xx <- x[which(x[,"ethnicity"]=="White" ),]
dim(xx) ## 1573
head(xx)
write.table(cbind(xx$MaskID,xx$LabID),f,quote=F,row.names=F,col.names=F)

f <- "icaps/icaps.thiazide.black/staring_samples.txt"
xx <- x[which(x[,"ethnicity"]=="Black"),]
dim(xx) ##617
head(xx)
write.table(cbind(xx$MaskID,xx$LabID),f,quote=F,row.names=F,col.names=F)

f <- "icaps/icaps.thiazide.hispanic/staring_samples.txt"
xx <- x[which(x[,"ethnicity"]=="Hispanic"),]
dim(xx) ## 154
head(xx)
write.table(cbind(xx$MaskID,xx$LabID),f,quote=F,row.names=F,col.names=F)

f <- "icaps/icaps.ccb.white/staring_samples.txt"
f <- "icaps/icaps.ccb.black/staring_samples.txt"
f <- "icaps/icaps.ccb.hispanic/staring_samples.txt"

f <- "icaps/icaps.bb.white/staring_samples.txt"
f <- "icaps/icaps.bb.black/staring_samples.txt"
f <- "icaps/icaps.bb.hispanic/staring_samples.txt"

f <- "icaps/icaps.aceia.arb.white/staring_samples.txt"
f <- "icaps/icaps.aceia.arb.black/staring_samples.txt"
f <- "icaps/icaps.aceia.arb.hispanic/staring_samples.txt"

write.table(cbind(data$MaskID,data$LabID),f,quote=F,row.names=F,col.names=F)
#pull white subjects for change in hba1c after metformin trtmnt
#d.samps <- sulf.score[which(sulf.score[,"Score"]==3),"MaskID"]
#x<-data[data$arm%in%c(5,7),]
#x.samps <- intersect(data[,"MaskID"],d.samps)
#x <-x[which(x[,"MaskID"]%in%x.samps),]

#head(x)
f <- "icaps/icaps_black/starting_samples.txt"
write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)

x <- data[which(data[,"ethnicity"]=="Hispanic"),]
#head(x)
f <- "icaps/icaps_hispanic/starting_samples.txt"
write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)
'''

## 4/6/2017 HbA1c Baseline
##################################

f<-"hba1c_blr/hba1c_blr_allraces/starting_samples.txt"
write.table(cbind(data$MaskID,data$LabID),f,quote=F,row.names=F,col.names=F)

x <- data[which(data[,"ethnicity"]=="White"),]
#head(x)
f<-"hba1c_blr/hba1c_blr_white/starting_samples.txt"
write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)

x <- data[which(data[,"ethnicity"]=="Black"),]
#head(x)
f<-"hba1c_blr/hba1c_blr_black/starting_samples.txt"
write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)

x <- data[which(data[,"ethnicity"]=="Hispanic"),]
#head(x)
f<-"hba1c_blr/hba1c_blr_hispanic/starting_samples.txt"
write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)



## 6/19/2017 Metformin HbA1c Baseline
##################################

f<-"met_hba1c_blr/met_hba1c_blr_allraces/starting_samples.txt"
met.data <- data[which(data$biguanide==1),]
write.table(cbind(met.data$MaskID,met.data$LabID),f,quote=F,row.names=F,col.names=F)

x <- met.data[which(met.data[,"ethnicity"]=="White"),]
#head(x)
f<-"met_hba1c_blr/met_hba1c_blr_white/starting_samples.txt"
write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)

x <- met.data[which(met.data[,"ethnicity"]=="Black"),]
#head(x)
f<-"met_hba1c_blr/met_hba1c_blr_black/starting_samples.txt"
write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)

x <- met.data[which(met.data[,"ethnicity"]=="Hispanic"),]
#head(x)
f<-"met_hba1c_blr/met_hba1c_blr_hispanic/starting_samples.txt"
write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)


## 9/04/2017 Metformin HbA1c Baseline
##################################

f<-"no_met_hba1c_blr/no_met_hba1c_blr_allraces/starting_samples.txt"
met.data <- data[which(data$biguanide!=1),]
write.table(cbind(met.data$MaskID,met.data$LabID),f,quote=F,row.names=F,col.names=F)

x <- met.data[which(met.data[,"ethnicity"]=="White"),]
#head(x)
f<-"no_met_hba1c_blr/no_met_hba1c_blr_white/starting_samples.txt"
write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)

x <- met.data[which(met.data[,"ethnicity"]=="Black"),]
#head(x)
f<-"no_met_hba1c_blr/no_met_hba1c_blr_black/starting_samples.txt"
write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)

x <- met.data[which(met.data[,"ethnicity"]=="Hispanic"),]
#head(x)
f<-"no_met_hba1c_blr/no_met_hba1c_blr_hispanic/starting_samples.txt"
write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)

## 11/20/2017 Hypoglycemia
##################################
f<-"no_met_hba1c_blr/no_met_hba1c_blr_allraces/starting_samples.txt"
met.data <- data[which(data$biguanide!=1),]
write.table(cbind(met.data$MaskID,met.data$LabID),f,quote=F,row.names=F,col.names=F)


## 3/23/2018
## mindcog replication effort for Jasmin
## sample_list was created in create.pheno_data.r
###################################################



######################################################################
# October 04, 2018

##################################
#pull in TZD samples-using 12mo TZD score
phenopath <- file.path("/home/accord/data/pheno_data/tzd_phenos/12mo")
analysispath <- file.path("/home/accord/data/analysis/tzd_weight_gain")

tzd.score <- read.csv(file.path(phenopath, "tzdscores_weight_gain.csv"),header = TRUE) 

#pull all subjects
d.samps <- data[which(data[,"tzd_score_weight_gain"]==3),"MaskID"]
#x<-data[data$arm%in%c(5,7),]
x.samps <- intersect(data[,"MaskID"],d.samps)
x <-data[which(data[,"MaskID"]%in%x.samps),]
f <- file.path(analysispath, "tzd_allraces", "starting_samples.txt")
write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)


#pull white subjects
d.samps <- tzd.score[which(tzd.score[,"Score"]==3),"MaskID"]
#x<-data[data$arm%in%c(5,7),]
x.samps <- intersect(data[,"MaskID"],d.samps)
x <-x[which(x[,"MaskID"]%in%x.samps),]
x <- x[which(x[,"ethnicity"]=="White"),]
f <- file.path(analysispath, "tzd_white", "starting_samples.txt")
write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)


#pull black subjects
d.samps <- tzd.score[which(tzd.score[,"Score"]==3),"MaskID"]
x<-data#[data$arm%in%c(5,7),]
x.samps <- intersect(x[,"MaskID"],d.samps)
x <-x[which(x[,"MaskID"]%in%x.samps),]
x <- x[which(x[,"ethnicity"]=="Black"),]
f <- file.path(analysispath, "tzd_black", "starting_samples.txt")
write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)

#pull hispanic subjects
d.samps <- tzd.score[which(tzd.score[,"Score"]==3),"MaskID"]
x<-data#x<-data[data$arm%in%c(5,7),]
x.samps <- intersect(x[,"MaskID"],d.samps)
x <-x[which(x[,"MaskID"]%in%x.samps),]
x <- x[which(x[,"ethnicity"]=="Hispanic"),]
f<- file.path(analysispath, "tzd_hispanic", "starting_samples.txt")
write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)


######################################################################
# August 30, 2019

# Lipids Sample List

analysispath <- file.path("/home/accord/data/analysis")

data<-read.table("pheno_data_simvastatin_targeted.txt",header=T,stringsAsFactors=F,sep='\t')
#pull all subjects
d.samps <- data[,"MaskID"]
#x<-data[data$arm%in%c(5,7),]
x.samps <- intersect(data[,"MaskID"],d.samps)
x <-data[which(data[,"MaskID"]%in%x.samps),]
f <- file.path(analysispath, "simvastatin_targeted_allraces", "starting_samples.txt")
write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)

######################################################################
# October 11, 2019

# Lipids Sample List - On Simvastatin Only
analysispath <- file.path("/home/accord/data/analysis")

data<-read.table(file.path(analysispath, "pheno_data_simvastatin_targeted.txt"),header=T,stringsAsFactors=F,sep='\t')
#pull all subjects
d.samps <- data[,"MaskID"]
#x<-data[data$arm%in%c(5,7),]
x.samps <- intersect(data[,"MaskID"],d.samps)
x <-data[which(data[,"MaskID"]%in%x.samps),]
x <- x[x$on_simvistatin == 1, ]
f <- file.path(analysispath, "simvastatin_targeted", "simvastatin_targeted_onSimvaOnly", "starting_samples.txt")
write.table(cbind(x$MaskID,x$LabID),f,quote=F,row.names=F,col.names=F)


######### John House RHTN 11/7/19 #######################
library(tidyverse)
analysispath <- file.path("/home/accord/data/analysis")
data<-read.table(file.path(analysispath, "pheno_data_rhtn.txt"),header=T,stringsAsFactors=F,sep='\t')

### Exclusion Criteria = gfr < 30 and bmi >= 40
### We had seven missing gfr, but NONE with gfr <30, so kept missing ###
#data %>% dplyr::filter(gfr >= 30, bmi < 40 ) -> excluded.data

data %>% dplyr::select(MaskID, LabID) %>% 
  write.table(., file = file.path(analysispath,"rhtn","rhtn_combined","starting_samples.txt"), quote=F,row.names=F,col.names=F)

data %>% dplyr::filter(gender == 0) %>% 
  dplyr::select(MaskID, LabID) %>% 
  write.table(., file = file.path(analysispath,"rhtn","rhtn_male","starting_samples.txt"), quote=F,row.names=F,col.names=F)

data %>% dplyr::filter(gender == 1) %>% 
  dplyr::select(MaskID, LabID) %>% 
  write.table(., file = file.path(analysispath,"rhtn","rhtn_female","starting_samples.txt"), quote=F,row.names=F,col.names=F)

##############################################






