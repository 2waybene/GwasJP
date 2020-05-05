#!/bin/bash
time R --slave --vanilla --file=bin/jleirer.scripts/get_peaks.r --args $1 $2 $3
