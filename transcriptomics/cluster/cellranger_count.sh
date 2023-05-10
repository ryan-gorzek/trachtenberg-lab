#### cellranger_count.sh START ####
#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblog.$JOB_ID
#$ -j y
## Edit the line below as needed:
#$ -l h_rt=12:00:00,h_data=32G,exclusive
## Modify the parallel environment
## and the number of cores as needed:
## -pe shared 1
# Email address to notify
#$ -M $USER@g.ucla.edu
# Notify when
#$ -m bea

# echo job info on joblog:
echo "Job $JOB_ID started on:   " `hostname -s`
echo "Job $JOB_ID started on:   " `date `
echo " "

# load the job environment:
. /u/local/Modules/default/init/modules.sh
# run the code
module load cellranger/7.1.0

cellranger count \
    --id=NW_TX0077-7_S01_L003-001 \
    --transcriptome=$HOME/genomes/Monodelphis_domestica_genome \
    --fastqs=$SCRATCH/fastqs/NW_TX0077-7_S01_L003-001/ \
    --include-introns=true

# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "
#### cellranger_count.sh STOP ####