##!/bin/sh 
#$ -N CN-3D
#$ -S /bin/sh 
#$ -cwd 
#$ -m bea
#$ -t 1-1:1
module  add sge shared gsl/gcc blas/gcc lapack/gcc gcc
#echo Starting at: `date`
./LS 
#echo Now it is: `date` 