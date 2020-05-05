import os

from . import __init__ as inW


def align(fastq_1, ref_fn, output_fn, fastq_2=None, p=1):
    '''
    Align reads using Bowtie2.
    '''
    assert os.path.exists(fastq_1)
    if fastq_2:
        assert os.path.exists(fastq_2)

    if fastq_2:
        inW.quiet_call([
            inW.PATHS['bowtie2'],
            '-q',
            '--phred33',
            '-p', str(p),
            '-I', '0',
            '-X', '1000',
            '--fr',
            '--local',
            '--sensitive-local',
            '-S', output_fn,
            '-x', ref_fn,
            '-1', fastq_1,
            '-2', fastq_2,
        ])

    else:
        inW.quiet_call([
            inW.PATHS['bowtie2'],
            '-q',
            '--phred33',
            '-p', str(p),
            '-I', '0',
            '-X', '1000',
            '--local',
            '--sensitive-local',
            '-S', output_fn,
            '-x', ref_fn,
            '-U', fastq_1,
        ])
