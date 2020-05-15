import tempfile, os


slurmInfo = []
slurmInfo.append("#!/bin/bash")
slurmInfo.append("#SBATCH --job-name=TRM000_phonetic.job")
slurmInfo.append("#SBATCH --output=TRM000_phonetic.job.out")
slurmInfo.append("#SBATCH --error=TRM000_phonetic.job.err")
slurmInfo.append("#SBATCH --time=2-00:00")
slurmInfo.append("#SBATCH --mem=12000")
slurmInfo.append("#SBATCH --qos=normal")
slurmInfo.append("#SBATCH --mail-type=END,FAIL")
slurmInfo.append("#SBATCH --mail-user=li11@niehs.nih.gov")

def getASLURMJob (jobFile, jobName, cmds, memory = 12000, runTimeallowed = "2-00:00" ):
    job    = jobName + ".job"
    output = jobName + ".out"
    errout = jobName + ".err"

    dir = tempfile.mkdtemp()
    file = dir+"/" + jobFile
    with open (file, 'w') as f:
        f.write("#!/bin/bash\n")
        f.write("#SBATCH --job-name=" + job + "\n")
        f.write("#SBATCH --output=" + output+ "\n")
        f.write("#SBATCH --error=" + errout+ "\n")
        f.write("#SBATCH --time=" + runTimeallowed + "\n")
        f.write("#SBATCH --mem="  + str(memory) + "\n")

       # f.write("#SBATCH --time=2-00:00\n")
       # f.write("#SBATCH --mem=12000\n")

        f.write("#SBATCH --qos=normal\n")
        f.write("#SBATCH --mail-type=END,FAIL\n")
        f.write("#SBATCH --mail-user=li11@niehs.nih.gov\n")
        for cmd in cmds:
            f.write( cmd + "\n")
    f.close()
    return(file,dir)

if __name__ == '__main__':
    (f,d) = getASLURMJob( "firsJob.sh", "firstJob", ["hello1", "hello2"])
    print(f)
    print(d)
