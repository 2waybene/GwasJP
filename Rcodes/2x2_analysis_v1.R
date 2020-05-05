

#create table
txt <- matrix(NA,ncol=2,nrow=2) 

#identify classes
class <- c("BCLL","TCLL")
#assign row and column names
rownames(txt) <- c(paste0(class,"_predicted"))
colnames(txt) <- c(paste0(class,"_truth"))
print(txt)
#assign values
txt[1,1] <- 22
txt[2,1] <- 3
txt[1,2] <- 6
txt[2,2] <- 40

print(txt)

## Summary ##
tot.num <- sum(txt)

### Calculates Sens & Spec for each combination###

output <- NULL
for(i in unique(class))
{
  if(which(class==i)==1)
  {
    tp <- txt[1,1]
    fn <- txt[2,1]
    tn <- txt[2,2]
    fp <- txt[1,2]
    
  }else{
    tp <- txt[2,2]
    fn <- txt[1,2]
    tn <- txt[1,1]
    fp <- txt[2,1]
  }
  
  sens <- round(tp/(tp+fn),digits = 3)
  spec <- round(tn/(tn+fp), digits = 3)
  ba <- mean(c(sens,spec))
  OR <- (tp/fp)/(fn/tn)
  RR <- (tp/(tp+fp)) / (fn/(tn+fn)) 
  TPR <- sens
  FPR <- 1-spec
  PPV <- tp/(tp+fp)
  NPV <- tn/(tn+fn)
  p.val <- fisher.test(txt)[[1]]
  
  output.temp <- c(i,tp,fp,fn,tn,tot.num,sens,spec,ba,OR,RR,TPR,FPR,PPV,NPV,p.val) 
  output <- rbind(output,output.temp)
}

all.data <- output[,7:ncol(output)]
class(all.data) <- "numeric"

overall <- colMeans(all.data)
output <- rbind(output,c("Average",rep(NA,5),overall))
colnames(output) <- c("Class","True.Positive","False.Positive","False.Negative","True.Negative","Total","Sens","Spec","BA","OR",
                      "RR","TPR","FPR","PPV","NPV","Fisher.p.val")


output.path <- "Y:\\NCSU_Projects\\Canine_Cancer\\leukemia\\Data\\Predictive_Model\\DecisionTree\\CLL-T_CLL-B"
output.name <- "2x2_SummaryStats_decision_tree_3mad_2.5-regions_CLL-T_CLL-B_03Apr2014.csv"
write.csv(output,file.path(output.path,output.name),row.names=FALSE)