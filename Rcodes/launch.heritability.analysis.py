import sys
import shlex
import subprocess as sp

def estimate_heritability():
    #print 'Number of arguments', len(sys.argv), 'arguments.'
    #print 'Argument List:', str(sys.argv)
    print '\n************************************'        
    print '******* Begin JOB:',sys.argv[0]

    for path in sys.argv[1:]:
        print '*************************************'
        print 'Current working path:',path
        ## Create system command
        cmd = ' '.join(('sbatch -p standard -o '+path+'/sbatch_logs/gcta.out ./bin/run_gcta.sh',path))
	print 'Launching heritability analysis:\n',cmd            
        ## Split cmd for shell
        split_cmd = shlex.split(cmd)
        ## Launch command
        sp.call(split_cmd)#,stdout=log_file,stderr=logerr_file)

    print '*************************************'        
    print '******* End JOB:',sys.argv[0]
    print '*************************************\n'
    
if __name__ == '__main__':
    estimate_heritability()


