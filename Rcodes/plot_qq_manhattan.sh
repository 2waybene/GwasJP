#!/bin/bash
p=$1
pheno=$2
maf=$3
#chr=$@
#R --slave --vanilla --file=bin/jj.scripts/plot_manhattan_only.r --args $p $pheno $maf
# R --slave --vanilla --file=bin/plot_qq_manhattan.r --args $p $pheno $maf
R --slave --vanilla --file=bin/jleirer.scripts/plot_qq_manhattan.r --args $p $pheno $maf
