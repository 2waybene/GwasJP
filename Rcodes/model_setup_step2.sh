#!/bin/bash
p=$1

echo "Remove designated covars and related individuals. Add first 10 PCs..."
# Remove selected covars and related individuals. Add first 10 PCs.
R --slave --vanilla --file=bin/pheno_data_step2.r --args $p

# Perform log transformation on pheno_data_step2.txt. Creates histograms and replaces vals in d4m cols
#R --slave --vanilla --file=bin/rotroff_scripts/log_transform_and_hist_v1.R --args $p

echo "Create modeltypes.txt. If only unique(phenotype values)=2, then logistic model is chosen..."
## Create a file called modeltypes.txt which explains the model type for each line of phenotypes.txt (lm or glm models)
R --slave --vanilla --file=bin/create.model.types.r --args $p

echo "Perform backwards selection on covars..."
# Backwards select non-forced covars. Create pheno files for R script, PLINK, and GCTA
R --slave --vanilla --file=bin/covar_backwards_selection_BIC.r --args $p

echo "Create samplelist.txt and frequency file..."
# Create sample list and frequency file
cut -f1,2 <(tail -n +2 $p/pheno_data/pheno_data_step2.txt) > $p/pheno_data/sample_list.txt
plink --bfile /home/accord/data/geno_data/post_qc.unc.uva.merged --keep $p/pheno_data/sample_list.txt --silent --freq --out $p/association_cv/plink
