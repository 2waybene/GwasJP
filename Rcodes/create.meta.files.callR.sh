#!/bin/bash
  
R --slave --vanilla --file=bin/create.meta.files.r --args $1 $2 $3
