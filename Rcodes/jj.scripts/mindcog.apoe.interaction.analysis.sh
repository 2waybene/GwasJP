#!/bin/bash
p=$1
pheno=$2
model=$3
snplist=$4

#echo $snplist
#echo $model
## Assume there is a snplist
plink --bfile /home/accord/data/geno_data/post_qc.unc.uva.merged \
	    --extract $p/snp_list.txt \
	    --linear --vif 1000 \
	    --maf 0.000001 \
	    --ci .95 \
	    --pheno $p/pheno_data/pheno_$pheno.txt \
	    --covar $p/apoe/covar_$pheno.txt \
	    --interaction \
	    --parameters 1,2,3,4,5,6,7,8,9,10,11,12,13,14,26,27 \
	    --noweb \
	    --out $p/apoe/association_cv/chr0.$pheno \