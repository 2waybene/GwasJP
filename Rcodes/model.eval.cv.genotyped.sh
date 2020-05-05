#!/bin/bash
p=$1
pheno=$2
model=$3
snplist=$4

#echo $snplist
#echo $model

if [ "$snplist" == "False" ] ## Run for all snps
then
    if [ "$model" == "linear" ]
    then
	echo "$pheno\n"
	echo '1'
	plink --bfile /home/accord/data/geno_data/post_qc.unc.uva.merged \
	    --linear --vif 1000 \
	    --maf 0.000001 \
	    --pheno $p/pheno_data/pheno_$pheno.txt \
	    --covar $p/pheno_data/covar_$pheno.txt \
	    --hide-covar \
	    --silent --noweb --out $p/association_cv/chr0.$pheno
    else ##  logistic
	echo "$pheno\n"
	echo '2'
	plink --bfile /home/accord/data/geno_data/post_qc.unc.uva.merged \
	    --logistic --vif 1000 \
	    --maf 0.000001 \
	    --1 \
	    --ci .95 \
	    --pheno $p/pheno_data/pheno_$pheno.txt \
	    --covar $p/pheno_data/covar_$pheno.txt \
	    --hide-covar \
	    --silent --noweb --out $p/association_cv/chr0.$pheno
    fi
else ## run only subset of snps 
    if [ "$model" == "linear" ]
    then
	echo "$pheno\n"
	echo '3'
	plink --bfile /home/accord/data/geno_data/post_qc.unc.uva.merged \
	    --extract $p/snp_list.txt \
	    --linear --vif 1000 \
	    --maf 0.000001 \
	    --pheno $p/pheno_data/pheno_$pheno.txt \
	    --covar $p/pheno_data/covar_$pheno.txt \
	    --hide-covar \
	    --silent --noweb --out $p/association_cv/chr0.$pheno
    else ## logistic
	echo "$pheno\n"
	echo '4'
	plink --bfile /home/accord/data/geno_data/post_qc.unc.uva.merged \
	    --extract $p/snp_list.txt \
	    --logistic --vif 1000 \
	    --maf 0.000001 \
	    --1 \
	    --ci .95 \
	    --pheno $p/pheno_data/pheno_$pheno.txt \
	    --covar $p/pheno_data/covar_$pheno.txt \
	    --hide-covar \
	    --silent --noweb --out $p/association_cv/chr0.$pheno
    fi
fi
