#### cellranger_mkref.sh START ####
#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblog.$JOB_ID
#$ -j y
## Edit the line below as needed:
#$ -l h_rt=24:00:00,h_data=32G,exclusive
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

cellranger mkgtf \
    $HOME/genomes/Monodelphis_domestica/Monodelphis_domestica.ASM229v1.109.gtf \
\
    $SCRATCH/Monodelphis_domestica.ASM229v1.109.filtered.gtf \
    --attribute=gene_biotype:protein_coding

cellranger mkref \
    --genome=Monodelphis_domestica_genome \
    --fasta=$HOME/genomes/Monodelphis_domestica/Monodelphis_domestica.ASM229v1.\
dna.toplevel.fa \
    --genes=$SCRATCH/Monodelphis_domestica.ASM229v1.109.filtered.gtf

# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "
#### cellranger_mkref.sh STOP ####