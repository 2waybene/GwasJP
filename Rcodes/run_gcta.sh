#!/bin/bash
p=$1
gcta64 --bfile /home/accord/data/geno_data/post_qc.unc.uva.merged \
       --keep $p/pheno_data/sample_list.txt \
       --autosome --make-grm --out $p/gcta/out
while read pheno;do if [ -n "${pheno}" ];then
  gcta64 --reml --grm $p/gcta/out \
	 --thread-num 8 \
         --pheno $p/gcta/pheno_${pheno}.txt \
         --covar $p/gcta/dcovar_${pheno}.txt \
         --qcovar $p/gcta/qcovar_${pheno}.txt \
         --out $p/gcta/out_${pheno}
fi;done < $p/phenotypes.txt

