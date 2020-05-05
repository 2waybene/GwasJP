#!/bin/bash
  
R --slave --vanilla --file=bin/jj.scripts/mindcog.merge.chunks.r --args $1 $2 $3
