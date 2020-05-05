conmed.score <- function(drug.list,drug.file,out.dir){
  
  cat("Reading in Concomitant Medications File...")
  conmeds.initial <- read.table('~/NCSU_Projects/Accord/Data/pheno_data/concomitantmeds.txt',sep="\t",header=TRUE)
  act.status <- read.table('~/NCSU_Projects/Accord/Data/pheno_data/activitystatus.txt', sep = '\t', header = TRUE)
  conmeds <- merge(conmeds.initial, act.status[,c("MaskID","Visit","exit_mapsto")], by = c("MaskID","Visit"), all.x = T)
  exits <- subset(conmeds, Visit == 'EXIT')
  exits$Visit <- exits$exit_mapsto
  conmeds <- rbind(conmeds, exits)
  conmeds <- conmeds[conmeds$Visit != 'EXIT',]
  conmeds <- conmeds[conmeds$Visit != '>7yrs',]
  conmeds$exit_mapsto <- NULL
  cat("Done.\n")
  cat("Reading in Medications File...")
  meds.file <- read.csv(drug.file,header=TRUE)
  cat("Done.\n")
  
  #tf.label <- unlist(strsplit(out.dir,split="/"))
  #tf.label <- unlist(strsplit(tf.label,split="\\",fixed=T))
  #tf.label <- tf.label[length(tf.label)]
  
  ind.list <- unique(meds.file[,"MaskID"])
  for(drug in drug.list){
    output <- NULL
    for(i in ind.list){
      cat("Calculating ",drug," For Subject:",i,grep(i,ind.list),"/",length(ind.list),"\n")
      target.time.frame <- as.character(meds.file[which(meds.file[,"MaskID"]==i),"Time_Frame"])
      time.frame <- unlist(strsplit(target.time.frame,split=":"))
      #drug.score <- as.character(meds.file[which(meds.file[,"MaskID"]==i),"Score"])
      if(length(time.frame)>1){
        t1 <- time.frame[1]
        t2 <- time.frame[2]
        ind.data <- conmeds[which(conmeds[,"MaskID"]==i),]
        ## update logic
        if(is.null(dim(ind.data))){
          response.time.frame <- NA
          score <- NA
        }else{
          un.times <- as.character(sort(unique(ind.data[,"Visit"])))
          if(length(grep("EXIT",un.times))>0){
            un.times <- un.times[-c(which(un.times=="EXIT"))]
            un.times <-c(un.times, "100")
          } 
          if(length(grep("BLR",un.times))>0){
            un.times <- un.times[-c(which(un.times=="BLR"))]
            un.times <-c("0",un.times)
          } 
          
          if(t1=="BLR")t1 <- "0"
          if(t2=="BLR")t2 <- "0"
          if(t1=="EXIT")t1 <- "1000"
          if(t2=="EXIT")t2 <- "1000"
          t1 <- as.numeric(gsub("F","",t1))
          t2 <- as.numeric(gsub("F","",t2))
          un.times <- as.numeric(gsub("F","",un.times))
          
          start.time <- max(which(un.times<=t1))
          end.time <- which(un.times<=t2)
          con.time <- unique(start.time,end.time)
          score <- as.numeric(as.character(ind.data[max(con.time),drug]))
          response.time.frame <- paste(ind.data[max(con.time),"Visit"],ind.data[max(con.time),"Visit"],sep=":") 
        
        }
        
      }else{
        response.time.frame <- NA
        score<- NA
      }
      output.temp <- c(i,target.time.frame,response.time.frame,score)
      output <- rbind(output,output.temp)
    }
    colnames(output) <- c("MaskID","Target_Time_Frame","Available_Time_Frame",paste0(drug,"Score"))
    write.table(output,file.path(out.dir,paste0(drug,"scores",".txt")),row.names=FALSE,quote=FALSE,sep="\t")
  }
  cat("***FINISHED***\n")
}


drug.list <- c("biguanide","aspirin","acei","tzd","vitamin","other_med","cholest_abi","beta_blocker","thiazide","sulfonylurea")#c('acei','meglitinide','ag_inhibitor','other_diabmed','loop','a2rb','cholest_abi','nitrate')
drug.file <- '~/NCSU_Projects/Accord/Data/pheno_data/statin_phenos/120d_drugScorePaper/statinscores.csv'#'Z:/accord/bin/Hillary/tzdScoreFiles/tzdscores.csv'
out.dir <- '~/NCSU_Projects/Accord/Data/pheno_data/statin_phenos/120d_drugScorePaper'#file.path('Z:','accord','bin','Hillary','tzdScoreFiles')
conmed.score(drug.list,drug.file,out.dir)
