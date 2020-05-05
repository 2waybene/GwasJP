#!/bin/bash
p=$1

R --slave --vanilla --file=bin/plot_gc.r --args $p