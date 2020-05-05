#!/bin/bash
#echo "Hello";
for p in ${pathList[@]};do
  #echo "Yo";
  echo $p;
  while read pheno;do if [ -n "${pheno}" ];then
	#echo "Yo1";
	#echo $p
        for chr in 22 112;do
		echo $chr;
		#echo "LAUNCH CHR:"$chr
		#echo "Yo2";
	        ./bin/jj.scripts/model_eval_cv_imputed.sh $p $chr > $p/sbatch_logs/chr$chr.LAUNCH.CHUNKS.cv.out &
        done;
  fi;done < $p/phenotypes.txt
done;

