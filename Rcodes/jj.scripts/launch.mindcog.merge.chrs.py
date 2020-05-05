import sys
import shlex
import subprocess as sp

def collapse_all_imputed_results():
    #print 'Number of arguments', len(sys.argv), 'arguments.'
    #print 'Argument List:', str(sys.argv)
    print '\n************************************'        
    print '******* Begin JOB:',sys.argv[0]

    for path in sys.argv[1:]:
        print '*************************************'
        print 'Current working path:',path

        ## Read in phenotypes.txt and modeltypes.txt
        phenos = [line.strip() for line in open(path+'/phenotypes.txt', 'r')]
        #models = [line.strip() for line in open(path+'/modeltypes.txt', 'r')]
        
        ## For each phenotype/modeltype, launch common variant analysis
        for i,pheno in enumerate(phenos):
		cmd = ' '.join(('sbatch -p long -o '+path+'/sbatch_logs/'+
				'mindcog.merge.chrs.'+pheno+
	                        '.out ./bin/jj.scripts/mindcog.merge.chrs.callR.sh',path,pheno))
                print 'Collapsing chromosomes into single file for phenotype',pheno+':\n',cmd         
		## Split cmd for shell
                split_cmd = shlex.split(cmd)
                ## Launch command
		sp.call(split_cmd)#,stdout=log_file,stderr=logerr_file)

    print '*************************************'        
    print '******* End JOB:',sys.argv[0]
    print '*************************************\n'
    
if __name__ == '__main__':
	collapse_all_imputed_results()


