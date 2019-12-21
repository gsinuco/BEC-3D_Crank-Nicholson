##!/bin/sh 
#$ -N CN-3D
#$ -S /bin/sh 
#$ -cwd 
# -M g.sinuco@sussex.ac.uk
#$ -m bea
#$ -t 1-6:5
module  add sge shared gsl/gcc blas/gcc lapack/gcc gcc
#echo Starting at: `date`
./zeeman_Floquet $SGE_TASK_ID
#echo Now it is: `date` 