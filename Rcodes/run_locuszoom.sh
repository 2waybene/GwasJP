#!/bin/bash
p=$1
pheno=$2
snp=$3

mf=$p/association_cv/allChrom_tzd_ch_wt_kg_gc_adj_METAL.tsv
locuszoom --metal=$mf --refsnp=$3 --flank 1000kb --build hg19 --pop EUR --build hg19 --source 1000G_March2012 --markercol SNP --pvalcol P_adj

#sbatch -o LZ.log locuszoom --metal allChrom_tzd_ch_wt_kg_gc_adj_METAL.tsv --refsnp rs80332660 --flank 500kb --build hg19 --pop EUR --build hg19 --source 1000G_March2012 --markercol SNP --pvalcol P_adj