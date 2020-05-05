#!/bin/bash

p=$1
pheno=$2
model=$3


if [ "$model" == "linear" ]
    then
	plink --meta-analysis $p/association_cv/allChrImputed_forMetaAnalysis.$pheno.assoc $p/association_cv/chr0.$pheno.assoc.$model + qt --silent --noweb --out $p/association_cv/plink_meta_$pheno
echo 'plink --meta-analysis '$p'/association_cv/allChrImputed_forMetaAnalysis.'$pheno'.assoc '$p'/association_cv/chr0.'$pheno'.assoc.'$model' + qt --silent --noweb --out '$p'/association_cv/plink_meta_'$pheno


else ##  logistic
	plink --meta-analysis $p/association_cv/allChrImputed_forMetaAnalysis.$pheno.assoc $p/association_cv/chr0.$pheno.assoc.$model --silent --noweb --out $p/association_cv/plink_meta_$pheno
echo 'plink --meta-analysis '$p'/association_cv/allChrImputed_forMetaAnalysis.'$pheno'.assoc '$p'/association_cv/chr0.'$pheno'.assoc.'$model' + logscale --silent --noweb --out '$p'/association_cv/plink_meta_'$pheno

fi
# plink --meta-analysis allChrImputed_forMetaAnalysis.met_hba1c.assoc chr0.met_hba1c.assoc.linear + qt\
#      --silent --noweb --out association_cv/plink_meta_met_hba1c
