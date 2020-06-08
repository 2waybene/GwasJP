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
@click.argument('inputdir',  type=str)
@click.argument('phenodata', type=str)
@click.argument('phenoname', type=str)
def accordModelStep1(rootdir, inputdir=None, phenodata=None, phenoname = None):

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

    click.echo(click.format_filename(rootdir))
    click.echo(click.format_filename(inputdir))
    click.echo(phenodata)


    inputdir = rootdir + "/" + inputdir
    fullPath = os.path.abspath(inputdir)
    print ("This is the full path:  " + fullPath)

    phenodata = rootdir + "/" + phenodata

    accord.modelStep1(fullPath, phenodata, phenoname)

@main.command()
#@click.argument('inputdir', type=click.Path(exists=True))
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
#@click.argument('inputdir', type=click.Path(exists=True))
@click.argument('rootdir', type=click.Path(exists=True))
@click.argument('inputdir',  type=str)
@click.option('--samplelist', default= "sample_list.txt", type=str,
              help='The sample_list.txt is needed for the analysis.')
@click.option('--thread', default= 8, type=int,
              help='The defaulty thread is 8')
def accordHeritability(rootdir, samplelist, thread,   inputdir=None):

    """
    Run accordJP pipeline,
    starting heritability analysis.

    As of this moment, JYL -- FIXME

    Print INPUTDIR if the directory exists.

    """
    click.echo(click.format_filename(rootdir))
    click.echo(click.format_filename(inputdir))

    inputdir = rootdir + "/" + inputdir
    fullPath = os.path.abspath(inputdir)
    print ("This is the full path:  " + fullPath)

    sampleList = fullPath + "/pheno_data/" + samplelist
    if (os.path.isfile(sampleList)):
        accord.heritabilityTest(fullPath, sampleList, thread)
    else:
        print ("Sample list file does no exist!\n")
        exit(1)

@main.command()
@click.argument('inputdir', type=click.Path(exists=True))
def accordGenoCommVar(inputdir):

    """

    Run accordJP pipeline, starting launch heritability test

    As of this moment, JYL -- FIXME

    Print INPUTDIR if the directory exists.

    """
    click.echo(click.format_filename(inputdir))
    accord.genoCommondVarAnalysis (inputdir)



if __name__ == "__main__":
    main()
