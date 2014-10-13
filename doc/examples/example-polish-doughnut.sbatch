#!/bin/bash
#SBATCH --job-name=Gyoto
#SBATCH --nodes=4 --ntasks-per-node=4

set -e

SCRATCH=/scratch/$USER/run.${SLURM_JOBID}

XMLFILE=$HOME/git/Gyoto/doc/examples/example-polish-doughnut.xml
FITSFILE=\!toto.fits

export LD_LIBRARY_PATH=/usr/local/lib:$HOME/usr/lib:$LD_LIBRARY_PATH
export PATH=$HOME/usr/bin:$PATH:$HOME/usr/lib/gyoto/3.0.1

mpirun -n 1 mkdir -p $SCRATCH
mpirun -n 1 cp $XMLFILE $SCRATCH/

mpirun -n 1 -wdir $SCRATCH ${HOME}/usr/bin/gyoto \
    --nprocesses=15 `basename $XMLFILE` $FITSFILE

mpirun -n 1 cp -a $SCRATCH $HOME/ && rm -Rf $SCRATCH

exit 0
