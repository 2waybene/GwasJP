#!/bin/bash
TIMEFORMAT='    %lR';
date

p=$1
phen=$2
tmprand=$RANDOM
me=$USER
time R --slave --vanilla --file=bin/jj.scripts/mindcog.interaction.model.setup.r --args $p $phen
