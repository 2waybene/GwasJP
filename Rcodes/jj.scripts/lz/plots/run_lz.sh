#!/bin/bash

# Plots are created in the current working directory.

for hit_file in /home/jj/data/accord/data/analysis/Sulfonylurea/sulf_90d_allraces2/peak_data/*_hitspec;
do
  echo "hitfile\n"
  echo ${hit_file}

  metal_file="${hit_file/_hitspec/_metal}"
  echo "metal\n"
  echo ${metal_file}

  file_name="${hit_file/*\//}"
  drug="${file_name/_hitspec/}"
  echo "drug\n"
  echo ${drug}

  /usr/local/locuszoom/bin/locuszoom --metal ${metal_file} --hitspec ${hit_file} --plotonly --gwas-cat whole-cat_significant-only --source=1000G_Nov2014 --build=hg19 --pop=EUR --prefix=${drug} showAnnot=T signifLine=6
done
