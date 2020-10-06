#!/bin/bash -l
#SBATCH --job-name=mkmapdata
#SBATCH --account=UMSU0010
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --partition=dav
#SBATCH --output=mkmapdata.out.%j

export TMPDIR=/glade/scratch/$USER/temp
mkdir -p $TMPDIR

export REGRID_PROC=1

# esmfvers=7.1.0r
# intelvers=17.0.1
# module purge
# module load intel/$intelvers
# module load ncl
# module load nco
# module load netcdf
# module load ncarcompilers
# module load esmflibs/$esmfvers

# module load esmf-${esmfvers}-ncdfio-uni-O

### Run program
# srun ./mksurfdata.pl -res usrspec -usr_gname 2minGLakes -usr_gdate 200928 -y 1850-2000 -dinlc /glade/p/cesm/cseg/inputdata/ -crop -hirespft -ssp_rcp hist -usr_mapdir /glade/work/felfelan/CTSM/clm5N/tools/mkmapdata/GLakes_200928 -inlandwet -dynpft ./landuse_timeseries_hist_78pfts_simyr1850-2015.txt
# srun ./mksurfdata.pl -res usrspec -usr_gname 2minGLakes -usr_gdate 200928 -y 1850-2000 -dinlc /glade/p/cesm/cseg/inputdata/ -crop -hirespft -ssp_rcp hist -usr_mapdir /glade/work/felfelan/CTSM/clm5N/tools/mkmapdata/GLakes_200928 -inlandwet

srun ./mksurfdata_map < ./GLakes_200930_y18502000_hirespft_ssp_rcphist_inlandwet/surfdata_2minGLakes_hist_78pfts_CMIP6_simyr1850_c200930.namelist.modnomask