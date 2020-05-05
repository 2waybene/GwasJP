#!/bin/bash
p=$1

# Filter SNPs in LD
plink --bfile $p/relatedness/data \
      --remove $p/relatedness/discard.txt \
      --maf 0.01 --indep 50 5 1.5 --silent --noweb --out $p/pca/data_maf_r2
plink --bfile $p/relatedness/data \
      --remove $p/relatedness/discard.txt \
      --extract $p/pca/data_maf_r2.prune.in \
      --recode12 --transpose --silent --noweb --out $p/pca/data_pruned

# Create ind, snp, and geno files
R --slave --vanilla --file=bin/pca_ind.r --args $p $2
awk '{print $2"\t"$1"\t0.0\t"$4}' $p/pca/data_pruned.tped > $p/pca/snp.txt
awk '{for (i=5;i<=NF;i=i+2) {j=i+1;v=$i+$j-2;if (v==-2) printf "%d",9;else printf "%d",v;};printf "\n";}' $p/pca/data_pruned.tped > $p/pca/geno.txt
#rm $p/pca/data_pruned*

# Compute PCs
smartpca.perl \
    -i $p/pca/geno.txt \
    -a $p/pca/snp.txt \
    -b $p/pca/ind.txt \
    -k 10 \
    -o $p/pca/result.pca \
    -p $p/pca/result.plot \
    -e $p/pca/result.eval \
    -l $p/pca/result.log \
    -m 0 \
    -t 5 \
    -s 6.0


#echo "Check this..."
#echo "R --slave --vanilla --file=bin/pca_plot.r --args $p $2"
# Plot PCs
R --slave --vanilla --file=bin/pca_plot.r --args $p
