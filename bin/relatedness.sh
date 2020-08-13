#!/bin/bash
p=$1

# Keep only those samples with phenotype data
cut -f 1-2 <(tail -n +2 $p/pheno_data/pheno_data_step1.txt) > $p/relatedness/keep.txt
plink --bfile ../geno_data/unc.jj/post_qc.v3 \
      --keep $p/relatedness/keep.txt \
      --silent --noweb --recode --make-bed --out $p/relatedness/data

# Compute and plot relatedness
bin/king -b $p/relatedness/data.bed \
     --kinship --related --degree 5 \
     --prefix $p/relatedness/king > $p/relatedness/king.log
R --slave --vanilla --file=bin/relatedness_plot.r --args $p

# Create discard list
R --slave --vanilla --file=bin/relatedness_discard.r --args $p
