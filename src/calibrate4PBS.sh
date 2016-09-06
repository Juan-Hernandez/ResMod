#!/bin/sh
#PBS -N calibrate4PBS
#PBS -l nodes=1:ppn=26
#PBS -l walltime=07:00:00:00

cd /data/scratch/m1jmh09/Julia/ResMod/src

/opt/julia/0.4.2/bin/julia calibrate4NM.jl > NelderMead4-2out.txt



