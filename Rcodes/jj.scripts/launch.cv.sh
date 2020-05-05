#!/bin/bash
time R --slave --vanilla --file=bin/compute_cv.r --args $1 $2 $3

