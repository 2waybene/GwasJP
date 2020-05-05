#!/bin/bash
time R --slave --vanilla --file=bin/jleirer.scripts/snp_gene_mapping.R --args $1 $2
