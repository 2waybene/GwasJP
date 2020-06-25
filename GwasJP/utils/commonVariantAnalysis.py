import os




##	coverted from model.eval.genotyped.sh
def modelEvalCVGenotyped (path, pheno, model, snpList, genotypeFle ="/ddn/gs1/home/li11/local/accord/data/geno_data/post_qc.unc.uva.merged"):

	phenoFile = path + "/pheno_data/pheno_" +  str(pheno) + ".txt"
	covarFile = path + "/pheno_data/covar_" +  str(pheno) + ".txt"
	outFile   = path + "/association_cv/chr0."  +  str(pheno)
	extractOrNot =  os.path.isfile(snpList)
	extractFile = path + str(snpList)
	cmd = ""
	if (extractOrNot == False):
		print ("NO selected SNP!\n")
		if (model == "liner"):
			print ("Using liner model for genotpye " + str (pheno) + "\n")
			print ("1")
			cmd = "plink --bfile " + genotypeFle + \
				  " --linear --vif 1000 --maf 0.000001 --pheno " + phenoFile + \
				  " --covar  " + covarFile + \
				  " --hide-covar --silent --noweb --out  " + outFile
		else:
		##  logistic
			print ("Using logistic model for genotpye " + str (pheno) + "\n")
			print ("2")
			cmd = "plink --bfile " + genotypeFle + \
				  " --logistic --vif 1000 --maf 0.000001 --1 --ci .95 --pheno " + phenoFile + \
				  " --covar  " + covarFile + \
				  " --hide-covar --silent --noweb --out  " + outFile
	else:
		print ("Here are the selected SNP!\n")
		if (model == "liner"):
			print ("Using liner model for genotpye " + str (pheno) + "\n")
			print ("3")
			cmd = "plink --bfile " + genotypeFle + \
				  " --extract " + extractFile +  \
				  " --linear --vif 1000 --maf 0.000001 --pheno " + phenoFile + \
				  " --covar  " + covarFile + \
				  " --hide-covar --silent --noweb --out  " + outFile
		else:
		##  logistic
			print ("Using logistic model for genotpye " + str (pheno) + "\n")
			print ("4")
			cmd = "plink --bfile " + genotypeFle + \
				  " --extract " + extractFile +\
				  " --logistic  --vif 1000 --maf 0.000001 --1 --ci .95 --pheno " + phenoFile + \
				  " --covar  " + covarFile + \
				  " --hide-covar --silent --noweb --out  " + outFile

	return(cmd)
