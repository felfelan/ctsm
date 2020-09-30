#!/bin/bash
### Job Name
#PBS -N mkmapdata
### Project code
#PBS -A UMSU0010
#PBS -l walltime=08:30:00
#PBS -q economy
### Merge output and error files
#PBS -j oe
#PBS -k eod
### Select 2 nodes with 36 CPUs each for a total of 72 MPI processes
#PBS -l select=1:ncpus=36:mpiprocs=36
### Send email on abort, begin and end
#PBS -m ae
### Specify mail recipient
#PBS -M felfelan@egr.msu.edu


export TMPDIR=/glade/scratch/$USER/temp
mkdir -p $TMPDIR

export REGRID_PROC=36

### Run program
mpiexec_mpt ./mkmapdata.sh -f ../mkmapgrids/SCRIPgrid_2minGLakes_nomask_modifiedNOfilled_c200925.nc -r 2minGlakes -t regional -v 