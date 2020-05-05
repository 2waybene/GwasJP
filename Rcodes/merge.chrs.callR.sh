#!/bin/bash
  
R --slave --vanilla --file=bin/merge.chrs.r --args $1 $2
