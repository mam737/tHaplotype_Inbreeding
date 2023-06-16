#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mmunasin@umn.edu
#SBATCH -p small,amdsmall
#SBATCH -t 96:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=40gb
#SBATCH -J Model1
#SBATCH --array=1-847

# set start time
d1=$(date +%s)

# set local directory name to unique job ID (unique for each task within a job)
newdir=${SLURM_JOB_ID}

# get parameter values for this task
myRow=`sed -n "${SLURM_ARRAY_TASK_ID}p" /path/to/Model1_params.txt`
mys=`sed -n "${SLURM_ARRAY_TASK_ID}p" /path/to/Model1_params.txt | cut -f 1`
myt=`sed -n "${SLURM_ARRAY_TASK_ID}p" /path/to/Model1_params.txt | cut -f 2`
myk=`sed -n "${SLURM_ARRAY_TASK_ID}p" /path/to/Model1_params.txt | cut -f 3`

# print task info to stdout (.out) file
echo host node
echo $HOSTNAME
echo task index
echo $SLURM_ARRAY_TASK_ID
echo job ID
echo $newdir
echo parameter set
echo $myRow

mkdir -p /path/to/scratch/$newdir
cp /path/to/Model1.R /path/to/scratch/$newdir
cd /path/to/scratch/$newdir
module load R/4.0.4

Rscript Model1.R ${mys} ${myt} ${myk}

# clean up
cd ..
rm -r ./$newdir

# print elapsed time
d2=$(date +%s)
sec=$(( ( $d2 - $d1 ) ))
hour=$(echo - | awk '{ print '$sec'/3600}')
echo Runtime: $hour hours \($sec\s\)
