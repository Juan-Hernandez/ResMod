#!/bin/sh
#PBS -N calsobolPBS
#PBS -l nodes=5:ppn=10
#PBS -l walltime=03:00:00:00

 ./etc/profile.d/modules.sh
module load Julia

julia calsobReentry.jl > hawkout.txt



