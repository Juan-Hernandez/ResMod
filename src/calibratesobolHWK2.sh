#!/bin/sh

#PBS -N calsobolPBS
#PBS -l nodes=1:ppn=12
#PBS -l walltime=03:00:00:00

#$ -pe threaded 12

#$ -V

#$ -cwd

./etc/profile.d/modules.sh 
module load julia

julia -p 12 calUSnomGR2.jl > hwkCALsobolOUT.txt



