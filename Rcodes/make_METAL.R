#####################
# Convert allChrom_[var]_gc_adj.csv to  allChrom_[var]_gc_adj_METAL.tsv
# Author: Jon Leirer
# Date: 04/19/2019
# 

# Something like ~/association_cv/
datapath <- getwd()

# Read in csv file.
df <- read.csv(file.path(datapath, "allChrom_tzd_hba1c_gc_adj.csv"))
head(df)
df2 <- df[,c("SNP", "P_unadj", "P_adj")]

allChrom_tzd_ch_wt_kg_gc_adj_METAL.tsv

write.table(df2, file.path(datapath, "allChrom_tzd_hba1c_gc_adj_METAL.tsv"), sep = "\t")
