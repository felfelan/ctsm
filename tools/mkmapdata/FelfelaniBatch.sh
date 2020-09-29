#!/bin/bash -l
#SBATCH --job-name=mkmapdata
#SBATCH --account=UMSU0010
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=08:30:00
#SBATCH --partition=dav
#SBATCH --output=mkmapdata.out.%j

export TMPDIR=/glade/scratch/$USER/temp
mkdir -p $TMPDIR

export REGRID_PROC=1

esmfvers=7.1.0r
intelvers=17.0.1
module purge
module load intel/$intelvers
module load ncl
module load nco
module load netcdf
module load ncarcompilers
module load esmflibs/$esmfvers

module load esmf-${esmfvers}-ncdfio-uni-O

### Run program
srun ./mkmapdata.sh -f ../mkmapgrids/SCRIPgrid_2minGLakes_nomask_modifiedNOfilled_c200925.nc -r 2minGlakes -t regional -v