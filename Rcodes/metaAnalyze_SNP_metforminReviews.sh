#!/bin/bash

#launch from association_cv directory
#paths run: Metformin/met_90d_allraces4_noGlyArmCorrect/association_cv
#paths run: Metformin/met_90d_white4_noGlyArmCorrect/association_cv


plink --meta-analysis imputed_chunks/imputed_chunks_forMeta/chr8.14.met_hba1c.assoc chr0.met_hba1c.assoc.linear + qt --silent --noweb --out plink_meta_met_hba1c

