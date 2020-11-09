import sys
import shlex
import subprocess as sp
import os

def create_files_for_meta():
    #print 'Number of arguments', len(sys.argv), 'arguments.'
    #print 'Argument List:', str(sys.argv)
    print '\n************************************'        
    print '******* Begin JOB:',sys.argv[0]

    for path in sys.argv[1:]:
        print '*************************************'
        print 'Current working path:',path
        #print '*************************************'
        ## Read in phenotypes.txt and modeltypes.txt
        phenos = [line.strip() for line in open(path+'/phenotypes.txt', 'r')]
        models = [line.strip() for line in open(path+'/modeltypes.txt', 'r')]
        #snplist =  os.path.isfile(path+'/snp_list.txt')
        
        #print phenos
        #print models
        #print snplist

        ## For each phenotype/modeltype, launch common variant analysis
        for i,pheno in enumerate(phenos):
            ## modeltype is passed as a parameter to the bash script
            cmd = ' '.join(('sbatch -p standard -o '+path+'/sbatch_logs/create.meta.files.'+
			pheno+'.out ./bin/create.meta.files.callR.sh',
                            path,pheno,models[i]))
            print 'Creating files for meta analysis for phenotype, ',pheno+':\n',cmd            
            ## Split cmd for shell
            split_cmd = shlex.split(cmd)
            ## Launch command
            sp.call(split_cmd)#,stdout=log_file,stderr=logerr_file)

    print '*************************************'        
    print '******* End JOB:',sys.argv[0]
    print '*************************************\n'
    
if __name__ == '__main__':
    create_files_for_meta()

