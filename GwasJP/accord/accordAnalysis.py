import sys
import shlex
import subprocess as sp
import os

from ..utils import statFittings, createSlurmJob, commonVariantAnalysis
from ..wrappers import gctaCalls,plinkCalls,smartpcaCalls


def modelSetupDirectories (fullPath, prerequisitesdir , projectname):

    print ("Copying presequiste files")
    cmdTmp = "cp " + prerequisitesdir+"/forced_covars.txt " + str(fullPath)
    sp.call(cmdTmp,  shell=True)
    cmdTmp = "cp " + prerequisitesdir+"/starting_covars.txt " + str(fullPath)
    sp.call(cmdTmp,  shell=True)
    cmdTmp = "cp " + prerequisitesdir+"/phenotypes.txt " + str(fullPath)
    sp.call(cmdTmp,  shell=True)
    cmdTmp = "cp " + prerequisitesdir+"/pheno_data_rhtn.txt " + str(fullPath)
    sp.call(cmdTmp,  shell=True)

    print ("Creating directories for " + str(fullPath) + " for project " + str(projectname))

    phenotype = str(fullPath) + "/phenotypes.txt"

    f = open(phenotype, 'r')
    phenoname = f.readline().strip()
    print("phenoname is " + str(phenoname) + "\n")
    #f.close()


    creatingDirs (fullPath, phenoname)

    print ("Finished directory setup, ready for analysis\n")

def modelStep1 (filepath, phenotype , Rdir , bFileInit ):

    print ("****** Begin JOB:' " + str(filepath) + "'")
    print ("****** This is the phenotype data info:' " + str(phenotype) + "'")

    #for path in filepath :
    print ('*************************************')
    print ('This is the working path entered from the user:', str(filepath))

  #  creatingDirs (filepath, phenoname)

    ##============================================================
    ## command 0: replacing model_setup_step1.sh
    ##============================================================

    ## prepare a file pheontypes.txt
   # phenotypes = filepath + "/" + "phenotypes.txt"

   # f = open(phenotype, 'w')
    #f.write(phenoname + "\n")
    #f.close()



    ##  ON Bioinformatic cluster at NIEHS
    ##============================================================
    ## command 1: used to be pheno_data_step1.r
    ##  Now located: /ddn/gs1/home/li11/local/accord/bin/pheno_data_step1.r
    ##============================================================

    outputFile = filepath + "/pheno_data/pheno_data_step1.txt"
    '''
    cmd1 = "R --slave --vanilla --file=/ddn/gs1/home/li11/local/accord/bin/pheno_data_step1.r --args " + filepath + " " + phenotype + " " +  outputFile
    sp.call(cmd1,  shell=True)
    '''
    cmd0 = "R --slave --vanilla --file=" + str(Rdir) + "/pheno_data_step1.r --args " + filepath + " " + phenotype + " " +  outputFile

    ##============================================================
    ## command 2-- : used to be time ./bin/relatedness.sh $p
    # echo;echo "Compute relatedness (bin/relatedness.sh)"
    # Now, re-write the file as command 2--
    ##=============================================================

    keptOut = filepath + "/relatedness/keep.txt"
    cmd1 = "cut -f 1-2 <(tail -n +2 " + outputFile + ") > " + keptOut
    ##=============================================================
    #for plink
    ##=============================================================
   # bFile = "/ddn/gs1/home/li11/local/accord/data/geno_data/unc.jj/post_qc.v3"

    outDir = filepath + "/relatedness/data"
    cmd2 = "plink --bfile " + bFileInit + " --keep " + keptOut+ "  --silent --noweb --recode --make-bed --out  " + outDir


    ##=============================================================
    #for kingship: king
    ##=============================================================
    bedFile = filepath + "/relatedness/data.bed"
    kPrefix = filepath + "/relatedness/king"
    kLog = filepath + "/relatedness/king.log"

    cmd3 = "king  -b " + bedFile + " --kinship --related --degree 5 --prefix " + kPrefix+ " > " + kLog

    ##=============================================================
    #for Compute and plot relatedness
    ##=============================================================
    cmd4 = "R --slave --vanilla --file=" + str(Rdir) + "/relatedness_plot.r  --args "+ filepath
    cmd5 = "R --slave --vanilla --file=" + str(Rdir) + "/relatedness_discard.r  --args " + filepath


    # Filter SNPs in LD
    ##=============================================================
    #for two plink analysis
    ##=============================================================

    bFile = filepath + "/relatedness/data"
    rmFile = filepath + "/relatedness/discard.txt"
    outDir = filepath + "/pca/data_maf_r2"
    cmd6 = "plink --bfile " + bFile + " --remove " + rmFile+ " --maf 0.01 --indep 50 5 1.5 --silent --noweb --out " + outDir

    extractFile = filepath + "/pca/data_maf_r2.prune.in"
    outPruned  = filepath + "/pca/data_pruned"
    cmd7 = "plink --bfile " + bFile + " --remove " + rmFile+  " --extract " + extractFile+ " --recode12 --transpose --silent --noweb --out " + outPruned

    cmd8 = "R --slave --vanilla --file=" + str(Rdir) + "/pca_ind.r  --args " + filepath + " " +  phenotype

    ##======================================================================
    ## two parts of awk script to parse files after pca
    ##======================================================================

    outPrunedTped  = filepath + "/pca/data_pruned.tped"
    snpFile = filepath + "/pca/snp.txt"

    cmd9 = "awk '{print $2\"\\t\"$1\"\\t0.0\\t\"$4}' " +  outPrunedTped + " > " + snpFile
    #  sp.call(cmdTemp,  shell=True, executable="/bin/bash")

    genoFile = filepath + "/pca/geno.txt"

    cmd10 = "awk '{for (i=5;i<=NF;i=i+2) {j=i+1;v=$i+$j-2;if (v==-2) printf \"%d\",9;else printf \"%d\",v;};printf \"\\n\";}' " + outPrunedTped + " > " + genoFile
    #    sp.call(cmdTemp,  shell=True, executable="/bin/bash")

    # Compute PCs
    oPCA =  filepath + "/pca/result.pca"
    pPCA =  filepath + "/pca/result.plot"
    ePCA =  filepath + "/pca/result.eval"
    lPCA =  filepath + "/pca/result.log"

    #cmd11 = "perl /ddn/gs1/home/li11/local/accord/bin/smartpca.perl -i " + genoFile + " -a " + snpFile + " -b " + filepath + "/pca/ind.txt" + " -k 10 -o " + oPCA \
    #        + " -p " + pPCA + " -e " + ePCA + " -l " + lPCA + " -m 0  -t 5   -s 6.0"
    cmd11 = "smartpca.perl -i " + genoFile + " -a " + snpFile + " -b " + filepath + "/pca/ind.txt" + " -k 10 -o " + oPCA \
            + " -p " + pPCA + " -e " + ePCA + " -l " + lPCA + " -m 0  -t 5   -s 6.0"
    # Plot PCs
    cmd12 = "R --slave --vanilla --file=" + str(Rdir) + "/pca_plot.r --args " + filepath


    ##======================================================================
    ## String the commands get a slurm file to submit to bioinfo cluster
    ##======================================================================

    commands = [cmd0,cmd1,cmd2,cmd3,cmd4,cmd5,cmd6,cmd7,cmd8,cmd9,cmd10,cmd11,cmd12]
    jobName = "modelsetupstep1"
    slurmSbatchFile="modelsetupstep1.sh"

    ## create a temporary sbatch file to submit
    (f,d) = createSlurmJob.getASLURMJob (slurmSbatchFile , jobName, commands, filepath)
    print (f)
    print(d)
    cmd = "sbatch --partition=bioinfo --cpus-per-task=8 " + f
    sp.call(cmd,  shell=True)

def creatingDirs (filepath, phenoname):

    dirBatch1 = ["association_cv",
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

    dirBatch2 = [
            "reg_plots/" +phenoname + "_call",
			"reg_plots/" +phenoname + "_call_bar",
			"reg_plots/" +phenoname + "_dosage",
			"reg_plots/" +phenoname + "_dosage_bar"
    ]

    dirs2make = []
    for dir in dirBatch1:
        dirs2make.append(filepath+"/"+dir)
    for dir in dirBatch2:
        dirs2make.append(filepath+"/"+dir)

    for dir in dirs2make:
        if (os.path.isdir(dir) == False):
            try:
                os.mkdir(dir)
                print("Directory '% s' created" % dir)
            except OSError as error:
                print(error)

def checkDirectories (filepath, phenoname):
    dirBatch1 = ["association_cv",
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

    dirBatch2 = [
            "reg_plots/" +phenoname + "_call",
			"reg_plots/" +phenoname + "_call_bar",
			"reg_plots/" +phenoname + "_dosage",
			"reg_plots/" +phenoname + "_dosage_bar"
    ]

    dirs2check = []
    for dir in dirBatch1:
        dirs2check.append(filepath+"/"+dir)
    for dir in dirBatch2:
        dirs2check.append(filepath+"/"+dir)

    for dir in dirs2check:
        if (os.path.isdir(dir) == False):
            print ("this is dir: " + str(dir) + "\n")
            return (1)
    return(0)


def modelStep2 (filepath, Rdir, bFileInit):
    print ("****** Begin JOB:' " + str(filepath) + "'")

    #for path in filepath :
    print ('*************************************')
    print ('This is the working path entered from the user:', str(filepath))

    ## Create system command

    #--file=/ddn/gs1/home/li11/local/accord/bin/pca_plot.r -
    #echo "Remove designated covars and related individuals. Add first 10 PCs..."
    # Remove selected covars and related individuals. Add first 10 PCs.
    cmd1 = "R --slave --vanilla --file=" + str(Rdir) + "/pheno_data_step2.r --args " + filepath

    # Perform log transformation on pheno_data_step2.txt. Creates histograms and replaces vals in d4m cols
    #R --slave --vanilla --file=bin/rotroff_scripts/log_transform_and_hist_v1.R --args $p

    #echo "Create modeltypes.txt. If only unique(phenotype values)=2, then logistic model is chosen..."
    ## Create a file called modeltypes.txt which explains the model type for each line of phenotypes.txt (lm or glm models)
    cmd2 = "R --slave --vanilla --file=" + str(Rdir) + "/create.model.types.r --args " + filepath

    #echo "Perform backwards selection on covars..."
    # Backwards select non-forced covars. Create pheno files for R script, PLINK, and GCTA

    rSourceFile = str(Rdir) + "/load_pheno_data.r"
    cmd3 = "R --slave --vanilla --file=" + str(Rdir) + "/covar_backwards_selection_BIC.r --args " + \
           filepath + "  " + rSourceFile

    #echo "Create samplelist.txt and frequency file..."
    #    Create sample list and frequency file

    outputFile = filepath + "/pheno_data/pheno_data_step2.txt"
    sampleList = filepath + "/pheno_data/sample_list.txt"

    cmd4 = "cut -f1,2 <(tail -n +2 " + outputFile + ") > " + sampleList
    #$p/pheno_data/pheno_data_step2.txt) > $p/pheno_data/sample_list.txt
   # bFile = "/home/accord/data/geno_data/post_qc.unc.uva.merged"

    plinkCV = filepath + "/association_cv/plink"

    cmd5 = "plink --bfile " + bFileInit + " --keep " + sampleList + " --silent --freq --out " + plinkCV

    #$p/association_cv/plink

    ##======================================================================
    ## String the commands get a slurm file to submit to bioinfo cluster
    ##======================================================================

    commands = [cmd1,cmd2,cmd3,cmd4,cmd5]
    jobName = "modelsetupstep2"
    slurmSbatchFile="modelsetupstep2.sh"

    ## create a temporary sbatch file to submit
    (f,d) = createSlurmJob.getASLURMJob (slurmSbatchFile , jobName, commands, filepath)
    print (f)
    print(d)
    cmd = "sbatch --partition=bioinfo --cpus-per-task=8 " + f
    sp.call(cmd,  shell=True)




def heritabilityTest (filepath, sampleList, phenotype,  p , genoTypeData):


    print ("****** Begin JOB:' " + str(filepath) + "'")
    #for path in filepath :
    print ('*************************************')
    print ('This is the working path entered from the user:', str(filepath))

    ## Create system command
	    ## ON NCSU cluter server
   # cmd = 'sbatch -p standard -o '+ genoTypeData +'  --keep " /sbatch_logs/gcta.out ./bin/run_gcta.sh  ' + filepath

    outputPath = filepath + "/gcta/out"
    cmd1 = "gcta64 --bfile " + genoTypeData + " --keep " + sampleList + " --autosome --make-grm --out  " + outputPath

    pheno = filepath + "/gcta/pheno_" + phenotype + ".txt"
    dcov  = filepath + "/gcta/dcovar_" + phenotype + ".txt"
    qcov  = filepath + "/gcta/qcovar_" + phenotype + ".txt"
    outdir = filepath + "/gcta/out_"   + phenotype
    cmd2 = "gcta64 --reml --grm " + outputPath + "  --thread-num " +  str(p)  + " --pheno " + pheno + " --covar " + dcov \
    + " --qcovar  " + qcov + "  --out " + outdir


    commands = [cmd1,cmd2]
    jobName = "heritability"
    slurmSbatchFile="accordHeritability.sh"

    ## create a temporary sbatch file to submit
    (f,d) = createSlurmJob.getASLURMJob (slurmSbatchFile , jobName, commands,filepath)
    print (f)
    print(d)

    cmd = "sbatch --partition=bioinfo --cpus-per-task=8 " + f
    sp.call(cmd,  shell=True)



  ## on Bioinfomatic slurm
    ## cmd = "srun --partition=bioinfo --cpus-per-task=8 -o  " + filepath + "/sbatch_logs/gcta.out  ./bin/run_gcta.sh  " + filepath
 #   print (cmd)
 #   sp.call(cmd,  shell=True)
    print ("Launching launchHeritability step 1 of 3:" + cmd)
    print ("Check the job status with command: squeue ")

def common_variant_analysis_genotyped (filepath, phenosFile, modelsFile, snplistFile , bFile):

 #   accord.common_variant_analysis_genotyped (fullPath, phenoname, modelfile, snplistFile)
    '''
    Current working path: RHTN_testRun/rhtn_combined/
    Launching logistic model for phenotype RHTN:
    sbatch -p bigmem -o RHTN_testRun/rhtn_combined//sbatch_logs/chr0.RHTN.out ./bin/model.eval.cv.genotyped.sh RHTN_testRun/rhtn_combined/ RHTN logistic False
    Submitted batch job 1498222

    '''
    print ("****** Begin JOB: Genotyped Common Variant Analysis  ******")
    print ("Here is the file path:  " + str(filepath) )

  #  phenosFile = filepath + "/" + str(phenosFile)

   # snplistFile = filepath + "/" + str(snplistFile)


    phenos = [line.strip() for line in open(phenosFile, 'r')]
    models = [line.strip() for line in open(modelsFile, 'r')]

   # snplist =  os.path.isfile(snplistFile)

    commands =[]
    clusterJobs = []
    ## For each phenotype/modeltype, launch common variant analysis
    for i,pheno in enumerate(phenos):
            ## modeltype is passed as a parameter to the bash script

            print ("This is this the i:  " + str(i))
            print ("This is the phenotype: " + pheno )
            cmdTemp = [commonVariantAnalysis.modelEvalCVGenotyped (filepath, pheno, models, snplistFile, bFile)]
            ##   using Default genotypeFle ="/home/accord/data/geno_data/post_qc.unc.uva.merged")

            jobName = "GenotypedCommonVariant" + str(i)
            slurmSbatchFile="GenotypedCommonVariant" + str(i) + ".sh"

            ## create a temporary sbatch file to submit
            (f,d) = createSlurmJob.getASLURMJob (slurmSbatchFile , jobName, cmdTemp, filepath )
            print (f)
            print(d)
            clusterJobs.append(f)


    print ("Launching impute common variant analysis  step 3 of 3:")
    print ("Check the job status with command: squeue ")

    for job in clusterJobs:
            cmd = "sbatch --partition=highmem --cpus-per-task=8 " + job
            sp.call(cmd,  shell=True)




def common_variant_analysis_imputed  (filepath, phenosFile, modelsFile, snplistFile, bFile):

   # phenosFile = filepath + "/" + str(phenosFile)
    #modelsFile = filepath + "/" + str(modelsFile)

   # snplistFile = filepath + "/" + str(snplistFile)

    #phenos = [line.strip() for line in open(phenosFile, 'r')]
    #models = [line.strip() for line in open(modelsFile, 'r')]

    phenos = [line.strip() for line in open(phenosFile, 'r')]
    models = [line.strip() for line in open(modelsFile, 'r')]

    chroms = list(range(1, 23, 1)) + list(range(101, 113, 1))
    print (chroms)

    commands =[]
    clusterJobs = []

    '''


    for i,pheno in enumerate(phenos):
            ## modeltype is passed as a parameter to the bash script
            for chrm in chroms:

                print ("This is this the i:  " + str(i))
                print ("This is the phenotype: " + pheno )
                cmdTemp = [commonVariantAnalysis.modelEvalCVGenotyped (filepath, pheno, models, snplistFile)]
            ##   using Default genotypeFle ="/home/accord/data/geno_data/post_qc.unc.uva.merged")
                cmdTemp = ' '.join(('sbatch -x node[1-9] -o '+path+'/sbatch_logs/chr'+str(chrm)+'.cv.for.meta.'+pheno+
                                '.out ./bin/model_eval_cv_imputed.for.meta.sh',
                                path,str(chrm),models[i]))
                jobName = "GenotypedCommonVariant" + str(i)
                slurmSbatchFile="GenotypedCommonVariant" + str(i) + ".sh"

            ## create a temporary sbatch file to submit
                (f,d) = createSlurmJob.getASLURMJob (slurmSbatchFile , jobName, cmdTemp, filepath )
                print (f)
                print(d)
                clusterJobs.append(f)

                cmd = ' '.join(('sbatch -x node[1-11] -o '+path+'/sbatch_logs/chr'+str(chrm)+'.cv.'+pheno+
                                  '.out ./bin/model_eval_cv_imputed.sh',
                                  path,str(chrm),models[i]))
                 print 'Launching full imputed analysis',models[i],'model for phenotype',pheno+':\n',cmd
                 ## Split cmd for shell
                 split_cmd = shlex.split(cmd)


                print ("This is this the i:  " + str(i))
                print ("This is the phenotype: " + pheno )
                cmdTemp = [commonVariantAnalysis.modelEvalCVGenotyped (filepath, pheno, models, snplistFile)]
            ##   using Default genotypeFle ="/home/accord/data/geno_data/post_qc.unc.uva.merged")

                jobName = "GenotypedCommonVariant" + str(i)
                slurmSbatchFile="GenotypedCommonVariant" + str(i) + ".sh"

            ## create a temporary sbatch file to submit
                (f,d) = createSlurmJob.getASLURMJob (slurmSbatchFile , jobName, cmdTemp, filepath )
                print (f)
                print(d)
                clusterJobs.append(f)

		## Impute Analysis for Meta Analysis
		 cmd = ' '.join(('sbatch -x node[1-9] -o '+path+'/sbatch_logs/chr'+str(chrm)+'.cv.for.meta.'+pheno+
                                '.out ./bin/model_eval_cv_imputed.for.meta.sh',
                                path,str(chrm),models[i]))
                 print 'Launching imputed analysis for meta',models[i],'model for phenotype',pheno+':\n',cmd
                 ## Split cmd for shell
                 split_cmd = shlex.split(cmd)
                 ## Launch command
                 sp.call(split_cmd)#,stdout=log_file,stderr=logerr_file)
    '''

    print ("Launching impute common variant analysis  step 3 of 3:")
    print ("Check the job status with command: squeue ")

    for job in clusterJobs:
            cmd = "sbatch --partition=highmem --cpus-per-task=8 " + job
            sp.call(cmd,  shell=True)



    print ("Check the job status with command: squeue ")

def cleanupImpuCommVarData (fullPath, phenotype, modelfile, selectedsnp):

    print ("Error checking and cleaning up data")

def metaAnalysis (filepath, phenosFile, modelsFile, snplistFile = None):

    #FIXME need to test

    phenosFile = filepath + "/" + str(phenosFile)
    modelsFile = filepath + "/" + str(modelsFile)
  #  snplistFile = filepath + "/" + str(snplistFile)
    phenos = [line.strip() for line in open(phenosFile, 'r')]
    models = [line.strip() for line in open(modelsFile, 'r')]


    snplist =  os.path.isfile(filepath+'/snp_list.txt')



        ## For each phenotype/modeltype, launch common variant analysis
    for i,pheno in enumerate(phenos):
            ## modeltype is passed as a parameter to the bash script
        cmd = ' '.join(('sbatch -p standard -o '+filepath+'/sbatch_logs/runMetaAnalysis.'+pheno+
		'.out ./bin/runPlinkMeta.sh',filepath,pheno,models[i],str(snplist)))


        '''
            #this is runPlinkMeta.sh
            if [ "$model" == "linear" ]
            then
	            plink --meta-analysis $p/association_cv/allChrImputed_forMetaAnalysis.$pheno.assoc $p/association_cv/chr0.$pheno.assoc.$model + qt --silent --noweb --out $p/association_cv/plink_meta_$pheno
                echo 'plink --meta-analysis '$p'/association_cv/allChrImputed_forMetaAnalysis.'$pheno'.assoc '$p'/association_cv/chr0.'$pheno'.assoc.'$model' + qt --silent --noweb --out '$p'/association_cv/plink_meta_'$pheno


            else ##  logistic
	            plink --meta-analysis $p/association_cv/allChrImputed_forMetaAnalysis.$pheno.assoc $p/association_cv/chr0.$pheno.assoc.$model --silent --noweb --out $p/association_cv/plink_meta_$pheno
                echo 'plink --meta-analysis '$p'/association_cv/allChrImputed_forMetaAnalysis.'$pheno'.assoc '$p'/association_cv/chr0.'$pheno'.assoc.'$model' + logscale --silent --noweb --out '$p'/association_cv/plink_meta_'$pheno

        '''
         #   print 'Launching meta analysis for phenotype',pheno+':\n',cmd
            ## Split cmd for shell
        split_cmd = shlex.split(cmd)
            ## Launch command
        sp.call(split_cmd)#,stdout=log_file,stderr=logerr_file)

    print ("needs to implement the meta analysis")

def getPlotting (fullPath, phenotype, modelfile, selectedsnp):

    print ("implement all the plottings ")

def rareVariantAnalysis(fullPath, phenotype, modelfile, selectedsnp):

    print ("Rare Variant Analysis")
