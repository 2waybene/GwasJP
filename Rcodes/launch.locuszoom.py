import sys
import shlex
import subprocess as sp

def model_setup_step1():
    #print 'Number of arguments', len(sys.argv), 'arguments.'
    #print 'Argument List:', str(sys.argv)
    print '\n************************************'        
    print '******* Begin JOB:',sys.argv[0]
    print "************************************"
    print 'LocusZoom plots:'
    for path in sys.argv[1:]:
	print '*************************************'
        print 'Current working path:',path
        #print '*************************************'
        ## Read in phenotypes.txt and modeltypes.txt
        phenos = [line.strip() for line in open(path+'/phenotypes.txt', 'r')]
        models = [line.strip() for line in open(path+'/modeltypes.txt', 'r')]
        #snplist =  os.path.isfile(path+'/snp_list.txt')
        #snp = sys.argv[2]
        #snp = "rs62373187"
        snp = "rs80332660"
        #snp = "rs12043219"

	for i,pheno in enumerate(phenos):
	        ## Create system command
            cmd = ' '.join(('sbatch -p bigmem -o '+path+'/sbatch_logs/qq.locuszoom.out ./bin/jleirer.scripts/run_locuszoom.sh',path,pheno,snp))
            print cmd
	        ## Split cmd for shell
            split_cmd = shlex.split(cmd)
	        ## Launch command
            sp.call(split_cmd)#,stdout=log_file,stderr=logerr_file)

    print '*************************************'        
    print '******* End JOB:',sys.argv[0]
    print '*************************************\n'
    
if __name__ == '__main__':
    model_setup_step1()