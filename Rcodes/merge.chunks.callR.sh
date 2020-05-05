#!/bin/bash
  
R --slave --vanilla --file=bin/merge.chunks.r --args $1 $2 $3
