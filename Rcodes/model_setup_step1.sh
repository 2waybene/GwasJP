#!/bin/bash

# Make sure dir path is input
if [ $# -ne 2 ];then echo "Must input path to script";exit 0;fi

# Make sure path exists
p=$1
if [ ! -d $p ];then echo "Path \"$p\" does not exist";exit 0;fi

# Format time, print date
TIMEFORMAT='    %lR';
date

# Create subdirs if they haven't been setup
if [ ! -d $p/pheno_data ];then
	mkdir $p/association_cv
	mkdir $p/association_cv/imputed_chunks
	mkdir $p/association_cv/imputed_chunks/imputed_chunks_forMeta
	mkdir $p/association_rv
	mkdir $p/cluster_plots
	mkdir $p/gcta
	mkdir $p/outputs
	mkdir $p/outputs/gc
	mkdir $p/pca
	mkdir $p/peak_data
	mkdir $p/pheno_data
	mkdir $p/relatedness
	mkdir $p/sbatch_logs
	mkdir $p/reg_plots
	while read pheno;do
		if [ -n "${pheno}" ];then 
			mkdir $p/reg_plots/${pheno}_call
			mkdir $p/reg_plots/${pheno}_call_bar
			mkdir $p/reg_plots/${pheno}_dosage
			mkdir $p/reg_plots/${pheno}_dosage_bar
		fi
	done < $p/phenotypes.txt
fi

echo;echo "Create complete cases phenotype data (bin/pheno_data_step1.r)"
R --slave --vanilla --file=bin/pheno_data_step1.r --args $p $2

echo;echo "Compute relatedness (bin/relatedness.sh)"
time ./bin/relatedness.sh $p

echo;echo "Compute PCs (bin/pca.sh)"
time ./bin/pca.sh $p $2
