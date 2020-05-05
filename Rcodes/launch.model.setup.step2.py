import sys
import shlex
import subprocess as sp

def model_setup_step2():
    #print 'Number of arguments', len(sys.argv), 'arguments.'
    #print 'Argument List:', str(sys.argv)
    print '\n************************************'        
    print '******* Begin JOB:',sys.argv[0]

    for path in sys.argv[1:]:
        print '*************************************'
        print 'Current working path:',path
        ## Create system command
        cmd = ' '.join(('sbatch -p standard -o '+path+'/model_setup_step2.out ./bin/model_setup_step2.sh',path))
	print 'Launching model setup step 2:\n',cmd            
        ## Split cmd for shell
        split_cmd = shlex.split(cmd)
        ## Launch command
        sp.call(split_cmd)#,stdout=log_file,stderr=logerr_file)

    print '*************************************'        
    print '******* End JOB:',sys.argv[0]
    print '*************************************\n'
    
if __name__ == '__main__':
    model_setup_step2()


