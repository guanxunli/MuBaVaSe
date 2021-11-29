#!/bin/sh
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH -t 12:00:00     # 12 hours
#SBATCH -p sched_mit_hill     # partition name
#SBATCH -J subset  # sensible name for the job

. /etc/profile.d/modules.sh
module add /cm/shared/modulefiles/engaging/gmp/6.1.1
R=/home/yuhaow/software/R/bin/R

