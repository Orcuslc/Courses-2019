#!/bin/sh
##Job Submission script

#$ -q UI
#$ -l h_rt=00:10:00
#$ -pe 32cpn 32
#$ -N sample1
#$ -o sample1.log
#$ -j y
#$ -V
#$ -cwd

echo "Job begin:"`date`
echo "Run sample1 on Argon, 32 proc"
mpirun -n 32 hydro.e sample1.in
echo "Job   end:"`date`




