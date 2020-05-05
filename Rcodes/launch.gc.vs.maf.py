import sys
import shlex
import subprocess as sp

def model_setup_step1():
    #print 'Number of arguments', len(sys.argv), 'arguments.'
    #print 'Argument List:', str(sys.argv)
    print '\n************************************'        
    print '******* Begin JOB:',sys.argv[0]
    print "************************************"
    print 'Plotting GC vs MAF:'
    for path in sys.argv[1:]:
        ## Create system command
        cmd = ' '.join(('sbatch -p bigmem -o '+path+'/sbatch_logs/plot.gc.out ./bin/plot_gc_callR.sh',path))
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


