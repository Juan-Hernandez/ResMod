#!/bin/sh

#PBS -N calsobolPBS
#PBS -l nodes=1:ppn=12
#PBS -l walltime=03:00:00:00

#$ -V

#$ -cwd

./etc/profile.d/modules.sh 
module load julia

julia calusGRy100.jl > hwkCALsobolOUT.txt



