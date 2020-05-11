import sys
import shlex
import subprocess as sp
import os

from ..utils import statFittings
from ..wrappers import gctaCalls,plinkCalls,smartpcaCalls


def modelStep1 (filepath, phenotype = "pheno_data_rhtn.txt", phenoname = "RHTN"):

    print ("****** Begin JOB:' " + str(filepath) + "'")
    print ("****** This is the phenotype data info:' " + str(phenotype) + "'")

    #for path in filepath :
    print ('*************************************')
    print ('This is the working path entered from the user:', str(filepath))

    creatingDirs (filepath, phenoname)
    ## replacing model_setup_step1.sh

    ## ON NCSU cluter server
    cmd = "sbatch -p standard -o " + filepath + "/model_setup_step1.out ./bin/model_setup_step1.sh  " + filepath + " " +   str(phenotype)

    ## prepare a file pheontypes.txt
    phenotypes = filepath + "/" + "phenotypes.txt"
    f = open(phenotypes, 'w')
    f.write(phenoname)
    f.close()

    cmd = "king"
    ## ON Bionformatic slurm system
    ## cmd = "srun --partition=bioinfo --cpus-per-task=8 -o  " + filepath + "/model_setup_step1.out ./bin/model_setup_step1.sh  " + filepath +  "  " + str(phenotype)
    print (cmd)
    sp.call(cmd,  shell=True)


    print ("Launching model setup step 1:" +  cmd)
    print ("Check the job status with command: squeue ")

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

def modelStep2 (filepath):

    print ("****** Begin JOB:' " + str(filepath) + "'")

    #for path in filepath :
    print ('*************************************')
    print ('This is the working path entered from the user:', str(filepath))

    ## Create system command


    ## ON NCSU cluter server

    cmd = 'sbatch -p standard -o '+filepath+'/model_setup_step2.out ./bin/model_setup_step2.sh  '  + filepath

    ## ON Bionformatic slurm system
    #3 cmd = "srun --partition=bioinfo --cpus-per-task=8 -o  " + filepath + "/model_setup_step2.out ./bin/model_setup_step2.sh  " + filepath
    print (cmd)
    sp.call(cmd,  shell=True)

    print ("Launching model setup step 2:" +  cmd)
    print ("Check the job status with command: squeue ")


def heritabilityTest (filepath):

    print ("****** Begin JOB:' " + str(filepath) + "'")
    #for path in filepath :
    print ('*************************************')
    print ('This is the working path entered from the user:', str(filepath))

    ## Create system command
	    ## ON NCSU cluter server
    cmd = 'sbatch -p standard -o '+filepath+'/sbatch_logs/gcta.out ./bin/run_gcta.sh  ' + filepath

  ## on Bioinfomatic slurm
    ## cmd = "srun --partition=bioinfo --cpus-per-task=8 -o  " + filepath + "/sbatch_logs/gcta.out  ./bin/run_gcta.sh  " + filepath
    print (cmd)
    sp.call(cmd,  shell=True)
    print ("Launching launchHeritability step 1 of 3:" + cmd)
    print ("Check the job status with command: squeue ")


def genoCommondVarAnalysis (filepath):

    print ("****** Begin JOB:' " + str(filepath) + "'")
    #for path in filepath :
    print ('*************************************')
    print ('This is the working path entered from the user:', str(filepath))

    ## Create system command

    # cmd = 'sbatch -p standard -o '+path+'/sbatch_logs/gcta.out ./bin/run_gcta.sh',path))
    cmd = "place holder"
    print (cmd)
    sp.call(cmd,  shell=True)
    print ("Launching genotype common variant analysis  step 2 of 3:" + cmd)
    print ("Check the job status with command: squeue ")


def imputeCommondVarAnalysis (filepath):

    print ("****** Begin JOB:' " + str(filepath) + "'")
    #for path in filepath :
    print ('*************************************')
    print ('This is the working path entered from the user:', str(filepath))

    ## Create system command

    # cmd = 'sbatch -p standard -o '+path+'/sbatch_logs/gcta.out ./bin/run_gcta.sh',path))
    cmd = "place holder"
    print (cmd)
    sp.call(cmd,  shell=True)
    print ("Launching impute common variant analysis  step 3 of 3:" + cmd)
    print ("Check the job status with command: squeue ")
