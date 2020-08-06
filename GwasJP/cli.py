import click
import sys
import os
#from . import  analysisPipeline
from . import accord
"""Console script for accordJP."""


@click.group()
def main(args=None):
    pass

@main.command()
@click.option('--shout/--no-shout', default=False)
def systemInfo(shout):
    '''
    Running accordJP pipeline, reporting the operating system.

    Default setting is --no-shout, gives lower case information.
    When --shout is provide, upper case information will be reported

    '''

    rv = sys.platform
    if shout:
        rv = rv.upper() + '!!!!'
    click.echo("This is your operating system: " + rv + ". Please pay enough attention!")

@main.command()
@click.argument('rootdir', type=click.Path(exists=True))
@click.argument('prerequisitesdir', type=click.Path(exists=True))
@click.argument('projectname', type=str)

def accordWorkingDirSetup(rootdir, prerequisitesdir, projectname):
    '''
    Running accordJP pipeline, this is step 1:

    rootdir is any directory that a user want to start
    projectname will be used as the "working directory" i.e. rhtn_combined etc.
    prerequeistes is a directory that contains all the prerequisite files,
    currently prepared outside (rhtn(s) are prepared by John House) the pipeline.

    /rundir/rhtn/pheno_data_rhtn.txt
    /rundir/rhtn/rhtn_combined/forced_covars.txt
    /rundir/rhtn/rhtn_combined/starting_covars.txt
    /rundir/rhtn/rhtn_combined/phenotypes.txt

    As of this moment, JYL -- FIXME

    '''

    click.echo(click.format_filename(rootdir))
    click.echo(click.format_filename(prerequisitesdir))
    click.echo(projectname)

    file1 = prerequisitesdir+"/forced_covars.txt"
    file2 = prerequisitesdir+"/starting_covars.txt"
    file3 = prerequisitesdir+"/phenotypes.txt"
    file4 = prerequisitesdir+"/pheno_data_rhtn.txt"

    if (os.path.exists(file1) and os.path.getsize(file1) > 0)  and \
            (os.path.exists(file2) and os.path.getsize(file2) > 0) and \
            (os.path.exists(file3) and os.path.getsize(file3) > 0) and \
            (os.path.exists(file4) and os.path.getsize(file4) > 0):
        print ("All prerequisites checked and passed!\n")
    else:
        print ("You don't seem to have the prequisite files, please check all the prerequisite files needed for the analysis!\n")
        exit(1)


    projectdir = rootdir + "/" + projectname
    if (os.path.isdir(projectdir) == False):
            try:
                os.mkdir(projectdir)
                print("Directory '% s' created" % projectdir)
            except OSError as error:
                print(error)

    fullPath = os.path.abspath(projectdir)
    accord.modelSetupDirectories (fullPath, prerequisitesdir, projectname)


@main.command()
@click.argument('runningdir', type=click.Path(exists=True))
#@click.argument('inputdir',  type=str)
#@click.argument('phenodata', type=str)
@click.argument('phenoname', type=str)

def accordModelStep1(runningdir,  phenoname = None):

    '''
    Running accordJP pipeline, step 1:

    launch model setup step 1

    As of this moment, JYL -- FIXME

    Instruction from Dr. John House:

    Make sure your phenotype data file is called '/home/accord/data/analysis/pheno_data.txt' before launching step 1
    Update: File should be named "pheno_data_*****.txt" Step 1 will prompt user for pheno file name
    Make sure you are in ~/accord/data/analysis when you start running this code

    Print INPUTDIR if the directory exists.
    Print PHENODATA from the input or use defaulty: pheno_data.txt

    '''

    click.echo(click.format_filename(runningdir))
  #  click.echo(click.format_filename(inputdir))
    click.echo(phenoname)

    file1 = runningdir+"/forced_covars.txt"
    file2 = runningdir+"/starting_covars.txt"
    phenotypes = runningdir+"/phenotypes.txt"
    phenodata  = runningdir+"/pheno_data_rhtn.txt"

    if (os.path.exists(file1) and os.path.getsize(file1) > 0)  and \
            (os.path.exists(file2) and os.path.getsize(file2) > 0) and \
            (os.path.exists(phenotypes) and os.path.getsize(phenotypes) > 0) and \
            (os.path.exists(phenodata) and os.path.getsize(phenodata) > 0):
        print ("All prerequisites checked and passed!\n")
    else:
        print ("You don't seem to have the prequisite files, please consider running accordWorkingDirSetup first or \
        manually set up directories and prepare all the prequisite files!\n")
        exit(1)

    if (accord.checkDirectories (runningdir, phenoname) == 1 ):
        print ("Subdirectory missing! Please check manuaal\n")

        exit(1)
  #  inputdir = rootdir + "/" + inputdir
    fullPath = os.path.abspath(runningdir)
    print ("This is the full path:  " + fullPath)

   # phenodata = rootdir + "/" + phenodata

    ## Kick off main analysis
    accord.modelStep1(fullPath, phenodata, phenoname, )


@main.command()
@click.argument('rootdir', type=click.Path(exists=True))
@click.argument('inputdir',  type=str)
def accordModelStep2(rootdir, inputdir=None):

    """
    Running accordJP pipeline, step 2:
    launch model setup step 2

    User needs to provide the directory as the argument to kick off the analysis.
    Print INPUTDIR if the directory exists.

    As of this moment, JYL -- FIXME
    """
    click.echo(click.format_filename(rootdir))
    click.echo(click.format_filename(inputdir))

    inputdir = rootdir + "/" + inputdir
    fullPath = os.path.abspath(inputdir)
    print ("This is the full path:  " + fullPath)

    accord.modelStep2(fullPath)

@main.command()
@click.argument('rootdir', type=click.Path(exists=True))
@click.argument('inputdir',  type=str)
@click.argument('phenoname', type=str)
@click.option('--samplelist', default= "sample_list.txt", type=str,
              help='The sample_list.txt is needed for the analysis.')
@click.option('--thread', default= 8, type=int,
              help='The defaulty thread is 8')
def accordHeritability(rootdir, samplelist, phenoname, thread, inputdir=None):

    """
    Run accordJP pipeline,
    starting heritability analysis.

    As of this moment, JYL -- FIXME

    Print INPUTDIR if the directory exists.

    """
    click.echo(click.format_filename(rootdir))
    click.echo(click.format_filename(inputdir))
    click.echo(click.format_filename(phenoname))

    inputdir = rootdir + "/" + inputdir
    fullPath = os.path.abspath(inputdir)
    print ("This is the full path:  " + fullPath)

    sampleList = fullPath + "/pheno_data/" + samplelist
    if (os.path.isfile(sampleList)):
        accord.heritabilityTest(fullPath, sampleList, phenoname, thread)
    else:
        print ("Sample list file does no exist!\n")
        exit(1)


@main.command()
@click.argument('rootdir', type=click.Path(exists=True))
@click.argument('inputdir',  type=str)
@click.argument('phenotype', type=str)
@click.argument('modelfile', type=str)
@click.argument('selectedsnp', type=str)

def accordGenoCommVar(rootdir, phenotype, modelfile, inputdir, selectedsnp=None):

    """

    Run accordJP pipeline, "genotyped common variant analysis"

    As of this moment, JYL -- FIXME

    Print INPUTDIR if the directory exists.

    """

    click.echo(click.format_filename(rootdir))
    click.echo(click.format_filename(inputdir))
    click.echo(click.format_filename(phenotype))
    click.echo(click.format_filename(modelfile))

    inputdir = rootdir + "/" + inputdir
    fullPath = os.path.abspath(inputdir)
    print ("This is the full path:  " + fullPath)
    click.echo(click.format_filename(inputdir))

    accord.common_variant_analysis_genotyped (fullPath, phenotype, modelfile, selectedsnp)


@main.command()
@click.argument('rootdir', type=click.Path(exists=True))
@click.argument('inputdir',  type=str)
@click.argument('phenotype', type=str)
@click.argument('modelfile', type=str)
@click.argument('selectedsnp', type=str)

def accordImpuCommVar(rootdir, phenotype, modelfile, inputdir, selectedsnp=None):

    """

    Run accordJP pipeline, "genotyped common variant analysis"

    As of this moment, JYL -- FIXME

    Print INPUTDIR if the directory exists.

    """

    click.echo(click.format_filename(rootdir))
    click.echo(click.format_filename(inputdir))
    click.echo(click.format_filename(phenotype))
    click.echo(click.format_filename(modelfile))

    inputdir = rootdir + "/" + inputdir
    fullPath = os.path.abspath(inputdir)
    print ("This is the full path:  " + fullPath)
    click.echo(click.format_filename(inputdir))

    accord.common_variant_analysis_imputed (fullPath, phenotype, modelfile, selectedsnp)


@main.command()
@click.argument('rootdir', type=click.Path(exists=True))
@click.argument('inputdir',  type=str)
@click.argument('phenotype', type=str)
@click.argument('modelfile', type=str)
@click.argument('selectedsnp', type=str)

def accordCleanupImpuCommVarData (rootdir, phenotype, modelfile, inputdir, selectedsnp=None):

    """

    Run accordJP pipeline, "imputed common variant analysis"

    As of this moment, JYL -- FIXME

    Print INPUTDIR if the directory exists.

    """

    click.echo(click.format_filename(rootdir))
    click.echo(click.format_filename(inputdir))
    click.echo(click.format_filename(phenotype))
    click.echo(click.format_filename(modelfile))

    inputdir = rootdir + "/" + inputdir
    fullPath = os.path.abspath(inputdir)
    print ("This is the full path:  " + fullPath)
    click.echo(click.format_filename(inputdir))

    accord.cleanupImpuCommVarData (fullPath, phenotype, modelfile, selectedsnp)


@main.command()
@click.argument('rootdir', type=click.Path(exists=True))
@click.argument('inputdir',  type=str)
@click.argument('phenotype', type=str)
@click.argument('modelfile', type=str)
@click.argument('selectedsnp', type=str)

def accordMetaAnalysis (rootdir, phenotype, modelfile, inputdir, selectedsnp=None):

    """

    Run accordJP pipeline, "meta analysis"

    As of this moment, JYL -- FIXME

    Print INPUTDIR if the directory exists.

    """

    click.echo(click.format_filename(rootdir))
    click.echo(click.format_filename(inputdir))
    click.echo(click.format_filename(phenotype))
    click.echo(click.format_filename(modelfile))

    inputdir = rootdir + "/" + inputdir
    fullPath = os.path.abspath(inputdir)
    print ("This is the full path:  " + fullPath)
    click.echo(click.format_filename(inputdir))

    accord.metaAnalysis (fullPath, phenotype, modelfile, selectedsnp)


@main.command()
@click.argument('rootdir', type=click.Path(exists=True))
@click.argument('inputdir',  type=str)
@click.argument('phenotype', type=str)
@click.argument('modelfile', type=str)
@click.argument('selectedsnp', type=str)

def accordDoThePlottings (rootdir, phenotype, modelfile, inputdir, selectedsnp=None):

    """

    Run accordJP pipeline, "doing plot"

    As of this moment, JYL -- FIXME

    Print INPUTDIR if the directory exists.

    """

    click.echo(click.format_filename(rootdir))
    click.echo(click.format_filename(inputdir))
    click.echo(click.format_filename(phenotype))
    click.echo(click.format_filename(modelfile))

    inputdir = rootdir + "/" + inputdir
    fullPath = os.path.abspath(inputdir)
    print ("This is the full path:  " + fullPath)
    click.echo(click.format_filename(inputdir))

    accord.getPlotting(fullPath, phenotype, modelfile, selectedsnp)


@main.command()
@click.argument('rootdir', type=click.Path(exists=True))
@click.argument('inputdir',  type=str)
@click.argument('phenotype', type=str)
@click.argument('modelfile', type=str)
@click.argument('selectedsnp', type=str)

def accordRareVariantAnalysis (rootdir, phenotype, modelfile, inputdir, selectedsnp=None):

    """

    Run accordJP pipeline, "rare variant analysis"

    As of this moment, JYL -- FIXME

    Print INPUTDIR if the directory exists.

    """

    click.echo(click.format_filename(rootdir))
    click.echo(click.format_filename(inputdir))
    click.echo(click.format_filename(phenotype))
    click.echo(click.format_filename(modelfile))

    inputdir = rootdir + "/" + inputdir
    fullPath = os.path.abspath(inputdir)
    print ("This is the full path:  " + fullPath)
    click.echo(click.format_filename(inputdir))

    accord.rareVariantAnalysis(fullPath, phenotype, modelfile, selectedsnp)


if __name__ == "__main__":
    main()
