import sys
import shlex
import subprocess as sp
import os

from ..utils import statFittings, createSlurmJob
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
   # cmd = "sbatch -p standard -o " + filepath + "/model_setup_step1.out ./bin/model_setup_step1.sh  " + filepath + " " +   str(phenotype)

    ## prepare a file pheontypes.txt
    phenotypes = filepath + "/" + "phenotypes.txt"
    f = open(phenotypes, 'w')
    f.write(phenoname)
    f.close()

    
    # echo;echo "Create complete cases phenotype data (bin/pheno_data_step1.r)"
    # R --slave --vanilla --file=bin/pheno_data_step1.r --args $p $2

    ##  ON Bioinformatic cluster at NIEHS
    # /ddn/gs1/home/li11/local/accord/bin/pheno_data_step1.r
    outputFile = filepath + "/pheno_data/pheno_data_step1.txt"

    cmd1 = "R --slave --vanilla --file=/ddn/gs1/home/li11/local/accord/bin/pheno_data_step1.r --args filepath phenotype outputFile"

    # echo;echo "Compute relatedness (bin/relatedness.sh)"
    # time ./bin/relatedness.sh $p
    #
    # echo;echo "Compute PCs (bin/pca.sh)"
    # time ./bin/pca.sh $p $2

    ##  finished r process and parse the pheno_data_step1.txt for the slurm job

    sp.call(cmd1,  shell=True)
    # Keep only those samples with phenotype data
    keptOut = filepath + "/relatedness/keep.txt"
    #getKeptRelatedness (outputFile, keptOut)

    cmdTemp = "cut -f 1-2 <(tail -n +2 " + outputFile + " ) > " + keptOut
    sp.call(cmdTemp,  shell=True)

    #cut -f 1-2 <(tail -n +2 $p/pheno_data/pheno_data_step1.txt) > $p/relatedness/keep.txt



    #for plink

    bFile = filepath + "/geno_data/unc.jj/post_qc.v3"
    outDir = filepath + "/relatedness/data/"


    cmd2 = "plink --bfile bFile --keep keptOut --silent --noweb --recode --make-bed --out  outDirz"

    bedFile = filepath + "/relatedness/data.bed"
    kPrefix = filepath + "/relatedness/king"
    kLog = filepath + "/relatedness/king.log"

    cmd3 = "king  -b bedFile    --kinship --related --degree 5 --prefix kPrefix > kLog "

# Compute and plot relatedness
    cmd4 = "R --slave --vanilla --file=/ddn/gs1/home/li11/local/accord/bin/relatedness_plot.r  --args filepath"
    cmd5 = "R --slave --vanilla --file=/ddn/gs1/home/li11/local/accord/bin/relatedness_discard.r  --args filepath"


#R --slave --vanilla --file=bin/relatedness_plot.r --args $p

# Create discard list
#R --slave --vanilla --file=bin/relatedness_discard.r --args $p

   # cmd = "king"
    ## ON Bionformatic slurm system
    ## cmd = "srun --partition=bioinfo --cpus-per-task=8 -o  " + filepath + "/model_setup_step1.out ./bin/model_setup_step1.sh  " + filepath +  "  " + str(phenotype)

    commands = [cmd2,cmd3,cmd4,cmd5]
    jobName = "modelsetupstep1"
    slurmSbatchFile="modelsetupstep1.sh"


    ## create a temporary sbatch file to submit

    (f,d) = createSlurmJob (slurmSbatchFile , jobName, commands)

    cmd =  "sbatch -p standard " + f
  #  print (cmd)
    sp.call(cmd,  shell=True)


  #  print ("Launching model setup step 1:" +  cmd)
  #  print ("Check the job status with command: squeue ")

# def getKeptRelatedness (outputFile, keptOut):



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
