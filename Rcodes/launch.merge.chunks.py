import sys
import shlex
import subprocess as sp

def collapse_imputed_chunks():
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
        
        ## Make list of chromosome chunks
        ## the list(range(101, 112, 1)) corresponds to XXY where XX is chunk and Y is chr
        chroms = list(range(1, 23, 1))
	print chroms
        ## For each phenotype/modeltype, launch common variant analysis
        for i,pheno in enumerate(phenos):
            for chrm in chroms:                
                #cmd = ' '.join(('sbatch -p standard -o '+path+'/sbatch_logs/chr0.'+pheno+'.%j.out ./bin/model.eval.cv.genotyped.sh',
                #path,pheno,models[i]))
                cmd = ' '.join(('sbatch -p standard -x node[1-17] -o '+path+'/sbatch_logs/chr'+str(chrm)+'.merge.chunks.'+pheno+
                                '.out ./bin/merge.chunks.callR.sh',
                                path,pheno,str(chrm)))
                print 'Collapsing chromosome',str(chrm),'for phenotype',pheno+':\n',cmd            		
                ## Split cmd for shell
                split_cmd = shlex.split(cmd)
                ## Launch command
		sp.call(split_cmd)#,stdout=log_file,stderr=logerr_file)

    print '*************************************'        
    print '******* End JOB:',sys.argv[0]
    print '*************************************\n'
    
if __name__ == '__main__':
	collapse_imputed_chunks()


