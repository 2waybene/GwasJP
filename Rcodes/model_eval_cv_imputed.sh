#!/bin/bash
TIMEFORMAT='    %lR';
date

## Set max number threads for GLM (BLAS algorithm). When Ubuntu upgraded OpenBLAS,
## This script started using as many threads as there were cores on any node. That caused a massive slowdown
## 16 threads * 16 jobs = 196 threads and system was overloaded. This forces only 4 threads. Could even bump down to 1
export OPENBLAS_NUM_THREADS=4

p=$1
chr=$2
tmprand=$RANDOM
me=$USER
if [ $chr -ge 100 ]; then
	chr=$((chr-100))
  for f in $(ls ../imputation/outputs/chr${chr}.*gz | sort -nt. -k5);do
    chunk=${f%.gz}
    chunk=${chunk##*.}
    if [[ ( "$chr" -eq 1 && "$chunk" -ge 25 ) || \
          ( "$chr" -eq 2 && "$chunk" -ge 26 ) || \
          ( "$chr" -eq 3 && "$chunk" -ge 20 ) || \
          ( "$chr" -eq 4 && "$chunk" -ge 20 ) || \
          ( "$chr" -eq 5 && "$chunk" -ge 20 ) || \
          ( "$chr" -eq 6 && "$chunk" -ge 16 ) || \
          ( "$chr" -eq 7 && "$chunk" -ge 16 ) || \
          ( "$chr" -eq 8 && "$chunk" -ge 15 ) || \
          ( "$chr" -eq 9 && "$chunk" -ge 15 ) || \
          ( "$chr" -eq 10 && "$chunk" -ge 14 ) || \
          ( "$chr" -eq 11 && "$chunk" -ge 14 ) || \
          ( "$chr" -eq 12 && "$chunk" -ge 14 ) \
        ]]; then

      echo "Make tmp directory on node: /tmp/accord.party.$chr.$chunk.$tmprand"
      mkdir "/tmp/accord.party.$me.$chr.$chunk.$tmprand"
      f_gz="chr$chr.imputed.$chunk"
      echo "gunzip the file on node: /tmp/accord.party.$me.$chr.$chunk.$tmprand/$f_gz"
      zcat "../imputation/outputs/$f_gz.gz" > "/tmp/accord.party.$me.$chr.$chunk.$tmprand/$f_gz"
      time R --slave --vanilla --file=bin/compute_cv.r --args $p $chr $chunk "/tmp/accord.party.$me.$chr.$chunk.$tmprand/$f_gz"
      echo "Remove tmp directory on node: /tmp/accord.party.$me.$chr.$chunk.$tmprand"
      rm -r "/tmp/accord.party.$me.$chr.$chunk.$tmprand"
    fi
  done
else
  for f in $(ls ../imputation/outputs/chr${chr}.*gz | sort -nt. -k5);do
    chunk=${f%.gz}
    chunk=${chunk##*.}
    if [[ ( "$chr" -eq 1 && "$chunk" -lt 25 ) || \
          ( "$chr" -eq 2 && "$chunk" -lt 26 ) || \
          ( "$chr" -eq 3 && "$chunk" -lt 20 ) || \
          ( "$chr" -eq 4 && "$chunk" -lt 20 ) || \
          ( "$chr" -eq 5 && "$chunk" -lt 20 ) || \
          ( "$chr" -eq 6 && "$chunk" -lt 16 ) || \
          ( "$chr" -eq 7 && "$chunk" -lt 16 ) || \
          ( "$chr" -eq 8 && "$chunk" -lt 15 ) || \
          ( "$chr" -eq 9 && "$chunk" -lt 15 ) || \
          ( "$chr" -eq 10 && "$chunk" -lt 14 ) || \
          ( "$chr" -eq 11 && "$chunk" -lt 14 ) || \
          ( "$chr" -eq 12 && "$chunk" -lt 14 ) || \
          ( "$chr" -ge 13 ) \
        ]]; then
      echo "Make tmp directory on node: /tmp/accord.party.$me.$chr.$chunk.$tmprand"
      mkdir "/tmp/accord.party.$me.$chr.$chunk.$tmprand"
      f_gz="chr$chr.imputed.$chunk"
      echo "gunzip the file on node: /tmp/accord.party.$me.$chr.$chunk.$tmprand/$f_gz"
      zcat "../imputation/outputs/$f_gz.gz" > "/tmp/accord.party.$me.$chr.$chunk.$tmprand/$f_gz"
      time R --slave --vanilla --file=bin/compute_cv.r --args $p $chr $chunk "/tmp/accord.party.$me.$chr.$chunk.$tmprand/$f_gz"
      echo "Remove tmp directory on node: /tmp/accord.party.$me.$chr.$chunk.$tmprand"
      rm -r "/tmp/accord.party.$me.$chr.$chunk.$tmprand"
    fi
  done
fi
