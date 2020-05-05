#!/bin/bash
TIMEFORMAT='    %lR';
date

p=$1
chr=$2
maf=$3

time R --slave --vanilla --file=bin/compute_rv.r --args $p $chr $maf
