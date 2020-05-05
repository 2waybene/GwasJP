#!/bin/bash
TIMEFORMAT='    %lR';
date

p=$1
chr=$2
#f="../imputation/outputs/chr${chr}.4.imputed.gz"
chunk=4

      #echo $f
      echo "Make tmp directory on node: /tmp/accord.party.$chr.$chunk/"
      mkdir "/tmp/accord.party.$chr.$chunk/"
      f_gz="chr$chr.imputed.$chunk"
      echo "gunzip the file on node: /tmp/accord.party.$chr.$chunk/$f_gz"
      zcat "../imputation/outputs/$f_gz.gz" > "/tmp/accord.party.$chr.$chunk/$f_gz"
      #time R --slave --vanilla --file=bin/compute_cv.r --args $p $chr $chunk $f_gz
      #rm -r "/tmp/accord.party.$chr.$chunk/"
      #time R --slave --vanilla --file=bin/compute_cv_imputeForMetaAnalysis.r --args $p $chr $chunk

