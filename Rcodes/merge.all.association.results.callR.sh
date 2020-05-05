#!/bin/bash
  
R --slave --vanilla --file=bin/merge.all.association.results.r --args $1 $2 $3
