import click
import sys
#from . import  analysisPipeline
from . import accord
"""Console script for accordJP."""


@click.group()
def main(args=None):
    pass


@main.command()
@click.argument('inputdir', type=click.Path(exists=True))
@click.argument('phenodata', type=str)
def accordModelStep1(inputdir, phenodata=None):

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


    click.echo(click.format_filename(inputdir))
    click.echo(phenodata)
    accord.modelStep1(inputdir, phenodata)

@main.command()
@click.argument('inputdir', type=click.Path(exists=True))
def accordModelStep2(inputdir):

    """
    Running accordJP pipeline, step 2:
    launch model setup step 2

    User needs to provide the directory as the argument to kick off the analysis.
    Print INPUTDIR if the directory exists.

    As of this moment, JYL -- FIXME
    """
    click.echo(click.format_filename(inputdir))
    accord.modelStep1(inputdir)


@main.command()
@click.argument('inputdir', type=click.Path(exists=True))
def accordHeritability(inputdir):

    """
    Run accordJP pipeline,
    starting heritability analysis.



    As of this moment, JYL -- FIXME

    Print INPUTDIR if the directory exists.

    """
    click.echo(click.format_filename(inputdir))
    accord.heritabilityTest(inputdir)

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
