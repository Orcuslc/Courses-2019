#!/bin/sh
##Job Submission script

#$ -q all.q
#$ -l h_rt=00:05:00
#$ -pe orte 1
#$ -N hw6
#$ -o hw6.log
#$ -j y
#$ -V
#$ -cwd

echo "Job begin:"`date`
echo "Run hw6 on Argon, 1 proc"
./hw6.e
echo "Job   end:"`date`




