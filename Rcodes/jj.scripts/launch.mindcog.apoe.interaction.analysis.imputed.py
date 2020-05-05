import sys
import shlex
import subprocess as sp

## n.b.,  Chromosomes 1-12 are split into 2 sections for imputed data,
##        chr=101-112 corresponds to second half of chromosomes 1-12.

def common_variant_analysis_imputed():
    #print 'Number of arguments', len(sys.argv), 'arguments.'
    #print 'Argument List:', str(sys.argv)
    print '\n************************************'
    print '******* Begin JOB:',sys.argv[0]

    for path in sys.argv[1:]:
        print '*************************************'
        print 'Current working path:',path

        ## Read in phenotypes.txt and modeltypes.txt
        phenos = [line.strip() for line in open(path+'/phenotypes.txt', 'r')]
        models = [line.strip() for line in open(path+'/modeltypes.txt', 'r')]

        ## Make list of chromosome chunks
        ## the list(range(101, 112, 1)) corresponds to XXY where XX is chunk and Y is chr
	chroms = list(range(1, 23, 1)) + list(range(101, 113, 1))
	#chroms = list(range(8, 10, 1))
	print chroms
        ## For each phenotype/modeltype, launch common variant analysis
        for i,pheno in enumerate(phenos):
            for chrm in chroms:
		## Impute Analysis
                 cmd = ' '.join(('sbatch -o '+path+'/sbatch_logs/mindcog.apor.chr'+str(chrm)+'.cv.'+pheno+
                                  '.out ./bin/jj.scripts/mindcog.apoe.interaction.analysis.imputed.sh',
                                  path,str(chrm),models[i]))
                 print 'Launching full imputed analysis',models[i],'model for phenotype',pheno+':\n',cmd
                 ## Split cmd for shell
                 split_cmd = shlex.split(cmd)
                 ## Launch command
                 sp.call(split_cmd)#,stdout=log_file,stderr=logerr_file)

		## Impute Analysis for Meta Analysis
		# cmd = ' '.join(('sbatch -p bigmem -x node[1-6,11,13,15] -o '+path+'/sbatch_logs/chr'+str(chrm)+'.cv.for.meta.'+pheno+
        #                        '.out ./bin/model_eval_cv_imputed.for.meta.sh',
        #                        path,str(chrm),models[i]))
        #         print 'Launching imputed analysis for meta',models[i],'model for phenotype',pheno+':\n',cmd
        #         ## Split cmd for shell
        #         split_cmd = shlex.split(cmd)
        #         ## Launch command
        #         sp.call(split_cmd)#,stdout=log_file,stderr=logerr_file)

    print '*************************************'
    print '******* End JOB:',sys.argv[0]
    print '*************************************\n'

if __name__ == '__main__':
    common_variant_analysis_imputed()


