

class accordObj:

   'Common base class for an accord trial case'

   dir2make = ["association_cv",
	                    "association_cv/imputed_chunks",
	                    "association_cv/imputed_chunks/imputed_chunks_forMeta",
	                    "association_rv",
	                    "cluster_plots",
	                    "gcta",
	                    "outputs",
	                    "outputs/gc",
	                    "pca",
	                    "peak_data",
	                    "pheno_data",
	                    "relatedness",
	                    "sbatch_logs",
	                    "reg_plots"]
   dir2makeBatch2 = []

   def __init__(self, phenoType, phenoName):
      self.pType = phenoType
      self.pName = phenoName
      accordObj.dir2makeBatch2.append ("reg_plots/" + self.pName + "_call")
      accordObj.dir2makeBatch2.append ("reg_plots/" + self.pName+ "_call_bar")
      accordObj.dir2makeBatch2.append ("reg_plots/" + self.pName+ "_dosage")
      accordObj.dir2makeBatch2.append ("reg_plots/" + self.pName+ "_dosage_bar")


   def displayPhenotypeInfo(self):
     print ("Phenotype information is: " % self.pType)

