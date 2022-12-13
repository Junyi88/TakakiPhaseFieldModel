#!/bin/sh
#PBS -lwalltime=24:00:00
#PBS -lselect=1:ncpus=16:mem=2gb

module load gcc
module load mpi

cp -r $WORK/TakakiPhaseFieldModel/* $TMPDIR
cd $TMPDIR/TakakiModel
make && cd $TMPDIR/build
#./Test1Con.out ../jobs/InitialInputs.inp
cp $TMPDIR/TakakiModel/exec/*.out $TMPDIR/build
mpiexec -n 2 ./Test1Con.out test2 ./InitialInputs.inp ./InputParameters.inp

mkdir $WORK/$PBS_JOBID
cp -r * $WORK/$PBS_JOBID
