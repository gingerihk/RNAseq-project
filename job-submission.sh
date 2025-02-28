# Single job 

#### job-submission.sh START ####
#!/bin/bash
#$ -cwd   ###ensures that the job is run in the directory from which the job was submitted
# error = Merged with joblog
#$ -o joblog.$JOB_ID.$TASK_ID
#$ -j y
#$ -l h_rt=24:00:00,h_data=8G
#$ -pe shared 1


# Load the module system
. /u/local/Modules/default/init/modules.sh
module load 


## substitute the command to run your code
## in the two lines below:
echo '/usr/bin/time -v ./convert_fastq.sh'
path-to-file/convert_fastq.sh

#### job-submission.sh STOP ####


# Multiple job (Kallisto)

#### job-submission.sh START ####
#!/bin/bash
#$ -cwd   ###ensures that the job is run in the directory from which the job was submitted\
# error = Merged with joblog
#$ -o $SCRATCH/joblogs-kallisto/joblog.$JOB_ID.$TASK_ID
#$ -j y
#$ -l h_rt=01:00:00,h_data=4G
#$ -pe shared 1
#$ -t 1-80:1 ###array job with 80 tasks, tasks will be executed sequentially

# echo job info on joblog:
echo "Job $JOB_ID.$SGE_TASK_ID started on:   " `hostname -s`
echo "Job $JOB_ID.$SGE_TASK_ID started on:   " `date `
echo " "

# Load the module system
. /u/local/Modules/default/init/modules.sh
module load 


## substitute the command to run your code
## in the two lines below:
echo '/usr/bin/time -v ./kallisto.sh'
path-to-file/kallisto.sh

# echo job info on joblog:
echo "Job $JOB_ID.$SGE_TASK_ID ended on:   " `hostname -s`
echo "Job $JOB_ID.$SGE_TASK_ID ended on:   " `date `
echo " "
#### job-submission.sh STOP ####

