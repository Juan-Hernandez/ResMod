#!/bin/sh
#PBS -N calsobolPBS
#PBS -l nodes=1:ppn=26
#PBS -l walltime=07:00:00:00

cd /data/scratch/m1dsw00/JuanStuff/ReservesMod/src

/opt/julia/0.4.2/bin/julia calsobRho.jl > calRhoOut.txt

/opt/julia/0.4.2/bin/julia calsobLambda.jl > calLambdaOut.txt

/opt/julia/0.4.2/bin/julia calsobReentry.jl > calReentryOut.txt



