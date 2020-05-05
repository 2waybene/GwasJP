#########################################
# ACCORD Statin Scores
# Author: Jon Leirer (built on legacy code)
# Date: 08/14/2019
# Description: Generate Statin Drug Scores
# Changelog:
# 08/14/2019 - Added filepath variables for flexibility
# 08/15/2019 - Add in scoring parameters for readability
#            - 
# Input Files:
#  f16_lipidmedicationsmanagement.txt
#  concomitantmeds.txt
#  activitystatus.txt
#  lipids.txt



#----------------
# Load Packages
#----------------
library(zoo)
library(doBy)
library(plyr)
options(warn = -1)

#-------------
# Filepaths
#-------------
inputpath <- "/home/accord/data"
outputpath <- "/home/jleirer/Small_Projects/Lipids/data"

#-------------
# Load Data
#-------------
# read in all necessary data files:
dat.statinLog <- read.table(file.path(inputpath, "pheno_data/f16_lipidmedicationsmanagement.txt"),
                            header = T, sep = "\t")
dat.conMeds <- read.table(file.path(inputpath,"pheno_data/concomitantmeds.txt"), 
                          header = T, sep = "\t")
dat.actStat <- read.table(file.path(inputpath, "pheno_data/activitystatus.txt"), 
                          header = T, sep = "\t")
dat.ldl <- read.table(file.path(inputpath, 'pheno_data/lipids.txt'), 
                      header = T, sep = '\t')

#-------------
# Merge Data
#-------------
# merge files to contain all relevant data:
# length(unique(dat.complete$MaskID))
dat.complete <- merge(dat.statinLog, dat.conMeds[,c("MaskID","Visit","statin")], by = c("MaskID","Visit"), all = T)
dat.complete <- merge(dat.complete, dat.actStat[,c('MaskID','Visit','days_from_baseline')], by = c('MaskID','Visit'), all.x = T)
dat.complete <- merge(dat.complete, dat.ldl, by = c('MaskID','Visit'), all = T)

# Change Class of variables
dat.complete$MaskID <- factor(dat.complete$MaskID)
dat.complete$ldl <- as.numeric(dat.complete$ldl)
dat.complete$days_from_baseline <- as.numeric(as.character(dat.complete$days_from_baseline))
dat.complete$statin <- as.numeric(as.character(dat.complete$statin))

# Remove observations missing "days from baseline"
dat.complete <- dat.complete[!is.na(dat.complete$days_from_baseline),]
dat.complete$comp <- as.numeric(as.character(dat.complete$zadhere))

#dat.complete$action <- as.numeric(as.character(dat.complete$action))
dat.complete <- merge(dat.complete, dat.actStat[,c("MaskID","Visit","exit_mapsto")], by = c("MaskID","Visit"), all.x = T)
exits <- subset(dat.complete, Visit == 'EXIT')
exits$Visit <- exits$exit_mapsto
dat.complete <- rbind(dat.complete, exits)
dat.complete <- dat.complete[dat.complete$Visit != 'EXIT',]
dat.complete <- dat.complete[dat.complete$Visit != '>7yrs',]
dat.complete$exit_mapsto <- NULL



# Score 0 subset 
# Create indicator for metformin present
temp <- dat.complete
temp$zadhere[which(temp$zadhere=='.')] <- 0
temp$statInd[temp$zadhere %in% c(1,4)] <- 1
temp$statInd[!temp$zadhere %in% c(1,4)] <- 0
temp$statConcInd[temp$statin!=1] <- 0
temp$statConcInd[temp$statin==1] <- 1
temp$statConcInd[is.na(temp$statConcInd)] <- 0

#temp$biguanide[is.na(temp$biguanide)] <- 0
# collapse data to identify score 0 patients
temp <- summaryBy(statInd + statConcInd ~ MaskID, FUN=sum, data=temp)
# if metInd = 0 and biguanide = 0 then we have score = 0 (patient was never on metformin)
dat.score0 <- subset(temp, statInd.sum == 0 & statConcInd.sum == 0, select = 'MaskID')
# include the rest of the score file variables
dat.score0$drug <- 'statin'
dat.score0$Time_Frame <- NA
dat.score0$Washout_Length_Months <- NA
dat.score0$Score <- 0


# DEAL WITH THE OTHER CASES STILL (3959 patients with no zadhere record of statin but record on conMed--MOSTLY DUE TO TRIAL ARM FOR STATIN)
# if have a metformin record on conMed not on baseline or exit visit, we toss out the patient completely (9/10251 patients)
# NOTE: I did this by hand: may want to come up with some code in the future
temp1 <- subset(temp, statInd.sum == 0 & statConcInd.sum != 0)
view <- merge(temp1, dat.complete, by = "MaskID", all.x = T)
view <- view[order(view$MaskID, view$days_from_baseline),]  
#write.csv(view, file = 'Z:/accord/bin/Hillary/Metformin/test/inconsistentFiles.csv', row.names = FALSE)
dat.complete <- dat.complete[!dat.complete$MaskID %in% temp1$MaskID,]



# Other subsets:
# filter down to patients who were on statin at any point in time (not score = 0 or tossed out patients)
dat.statPatients <- merge(dat.complete, dat.score0[,c("MaskID",'Score')], all = T)
dat.statPatients <- dat.statPatients[is.na(dat.statPatients$Score),]
#dat.metPatients$comp[dat.metPatients$oral_gmed != 'metformin']  <- NA
#dat.metPatients$action[dat.metPatients$oral_gmed != 'metformin']  <- NA
#dat.metPatients$oral_gmed[(dat.metPatients$oral_gmed != 'metformin') | dat.metPatients$oral_gmed == '']  <- NA
dat.statPatients$comp[(dat.statPatients$comp == 1 | dat.statPatients$comp == 4)] <- 1
dat.statPatients$comp[(dat.statPatients$comp == 2)] <- 0


# impute missing compliance values and classify patients
dat.statPatients <- dat.statPatients[order(dat.statPatients$MaskID, dat.statPatients$days_from_baseline),]
dat.patients <- c()
patient <- data.frame(matrix(ncol = 0, nrow = 1))

for (i in 1:length(unique(dat.statPatients$MaskID))) {
  print(paste(i,"of",length(unique(dat.statPatients$MaskID))))
  # pull out one patients observations
  temp <- subset(dat.statPatients, MaskID==unique(dat.statPatients$MaskID)[i])
  
  # Create Action Code indicator
#   temp$actStop <- 0
#   temp$actStop[temp$action == 2] <- 1
#   temp$actStart <- 0
#   temp$actStart[temp$action == 4] <- 1
#  


  temp$ind <- 0
  temp$ind[(temp$statin == 1 | temp$comp == 1)] <- 1
  
  # collapse down to one record per visit
  #temp <- ddply(temp, .(MaskID, Visit, days_from_baseline), summarize, ind = max(ind), actStart = max(actStart), actStop = max(actStop), met = max(biguanide, na.rm = TRUE), comp = max(comp, na.rm = TRUE), hba1c = max(hba1c, na.rm = TRUE))
  temp <- ddply(temp, .(MaskID, Visit, days_from_baseline), summarize, ind = max(ind),statin = max(statin, na.rm = TRUE), comp = max(comp, na.rm = TRUE), ldl = max(ldl, na.rm = TRUE))
  temp$comp[is.infinite(temp$comp)] <- NA
  temp$statin[is.infinite(temp$statin)] <- NA
  temp$ldl[is.infinite(temp$ldl)] <- NA
  
  # Create days variable and locate window 
  temp <- temp[order(temp$days_from_baseline),]
  for (i in 1:nrow(temp)){
    temp$first[i] <- sum(temp$ind[1:i])
  } 
  for (i in 1:nrow(temp)){
    temp$first2[i] <- sum(temp$first[1:i])
  } 
  temp <- temp[order(temp$days_from_baseline),]
  for (i in 1:nrow(temp)){
    temp$days[i] <- temp$days_from_baseline[i] - temp$days_from_baseline[temp$first2 == 1]
  }
  temp <- temp[(temp$days <= 270 & temp$days >= -30), ]
  temp$last.ind <- 0
  temp$last.ind[temp$days >= 90 & !is.na(temp$hba1c)] <- 1
  temp <- temp[order(temp$days),]
  for (i in 1:nrow(temp)){
    temp$l.window[i] <- sum(temp$last.ind[1:i])
  }
  for (i in 1:nrow(temp)){
    temp$l.window2[i] <- sum(temp$l.window[1:i])
  }
  temp <- temp[temp$l.window2 <= 1,] 
  
  temp$first.ind <- 0
  temp$first.ind[temp$days <= 0 & !is.na(temp$ldl)] <- 1
  temp <- temp[order(temp$days),]
  for (i in 1:nrow(temp)){
    temp$f.window[i] <- sum(temp$first.ind[i:nrow(temp)])
  }
  for (i in 1:nrow(temp)){
    temp$f.window2[i] <- sum(temp$f.window[i:nrow(temp)])
  }
  temp <- temp[temp$f.window2 <= 1,] 
  
  # Grab Start window value and then delete anything prior to when patient started statin
  patient$MaskID <- temp$MaskID[1]
  patient$Start <- temp$Visit[1]
  patient$Startldl <- temp$ldl[1]
  temp <- temp[temp$days >= 0,]
  
  # Identify those who were already on statin prior to window
  temp <- temp[order(temp$days),]
  temp$startPreTrial <- 0
  temp$startPreTrial[temp$statin[1] == 1] <- 1
  
  # impute missing compliance records
  temp$comp.fill <- na.locf(temp$comp, na.rm = FALSE, fromLast = TRUE)
  temp$stat.fill <- na.locf(temp$statin, na.rm = FALSE, fromLast = TRUE)   
  temp$comp.fill[is.na(temp$comp.fill)] <- 0
  temp$stat.fill[is.na(temp$stat.fill)] <- 0
  temp$notCompliant <- 0
  temp$notCompliant[(temp$comp.fill == 3)] <- 1
  
  # Determine whether the patient ever stopped taking metformin in the window
  temp$continuous <- 1
  #temp$continuous[(temp$actStart != temp$actStop)] <- 0
  temp$continuous[1] <- 1
  temp$continuous[nrow(temp)] <- 1
  
  # get relevant patient information in one row
  patient$Time_Frame <- paste0(patient$Start,':',temp$Visit[nrow(temp)])
  patient$Endldl <- temp$ldl[nrow(temp)]
  patient$trialLength <- temp$days[nrow(temp)]
  patient$notCompliant <- max(temp$notCompliant)
  patient$complianceAv <- mean(temp$comp.fill[2:nrow(temp)])
  patient$complianceAv[is.na(patient$complianceAv)] <- 1 
  patient$startPreTrial <- max(temp$startPreTrial)
  patient$continuous <- min(temp$continuous)
  
  dat.patients <- rbind(dat.patients,patient)
}


dat.throwout <- subset(dat.patients, notCompliant == 1 | continuous == 0 | complianceAv < .8 | is.na(Startldl) | is.na(Endldl))

dat.score1 <- subset(dat.patients, notCompliant == 0 & continuous == 1 & startPreTrial == 1 & trialLength < 90 & complianceAv >= .8 & !is.na(Startldl) & !is.na(Endldl))
dat.score1$drug <- 'statin'
dat.score1$Washout_Length_Months <- NA
dat.score1$Score <- 1
dat.score1 <- dat.score1[,c('MaskID','drug','Time_Frame','Washout_Length_Months','Score')]

dat.score2 <- subset(dat.patients, notCompliant == 0 & continuous == 1 & startPreTrial == 0 & trialLength < 90 & complianceAv >= .8 & !is.na(Startldl) & !is.na(Endldl))
dat.score2$drug <- 'statin'
dat.score2$Washout_Length_Months <- NA
dat.score2$Score <- 2
dat.score2 <- dat.score2[,c('MaskID','drug','Time_Frame','Washout_Length_Months','Score')]

dat.score3 <- subset(dat.patients, notCompliant == 0 & continuous == 1 & startPreTrial == 0 & trialLength >= 90 & complianceAv >= .8 & !is.na(Startldl) & !is.na(Endldl))
dat.score3$drug <- 'statin'
dat.score3$Washout_Length_Months <- NA
dat.score3$Score <- 3
dat.score3 <- dat.score3[,c('MaskID','drug','Time_Frame','Washout_Length_Months','Score')]

dat.score4 <- subset(dat.patients, notCompliant == 0 & continuous == 1 & startPreTrial == 1 & trialLength >= 90 & complianceAv >= .8 & !is.na(Startldl) & !is.na(Endldl))
dat.score4$drug <- 'statin'
dat.score4$Washout_Length_Months <- NA
dat.score4$Score <- 4
dat.score4 <- dat.score4[,c('MaskID','drug','Time_Frame','Washout_Length_Months','Score')]

statinscores <- rbind(dat.score0, dat.score1, dat.score2, dat.score3, dat.score4)
# write.csv(statinscores, file = 'pheno_data/statin_phenos/120d_drugScorePaper/statinscores.csv', row.names = FALSE)
write.csv(statinscores, file = file.path(outputpath, 'statinscores.csv'), row.names = FALSE)

