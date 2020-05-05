#!/bin/bash
echo "Make tmp directory on node: /tmp/accord.party/"
mkdir "/tmp/accord.party/"
f_gz="chr8.imputed.4"
echo "gunzip the file on node: /tmp/accord.party/$f_gz.gz"
zcat "../imputation/outputs/$f_gz.gz" > "/tmp/accord.party/$f_gz"
## Set max number threads for GLM (BLAS algorithm). When Ubuntu upgraded OpenBLAS,
## This script started using as many threads as there were cores on any node. That caused a massive slowdown
## 16 threads * 16 jobs = 196 threads and system was overloaded
export OPENBLAS_NUM_THREADS=4
echo "Starting R"
time R --slave --vanilla --file=bin/compute_cv.r --args "test.blas/icaps.ccb.all_races" "8" "4" "/tmp/accord.party/$f_gz"
echo "Remove tmp directory on node: /tmp/accord.party"
rm -r "/tmp/accord.party"
