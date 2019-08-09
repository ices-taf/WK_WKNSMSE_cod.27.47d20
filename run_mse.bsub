#!/bin/sh
## queue
#BSUB -q cefas-ib
## project code
#BSUB -P MA011A
## number of CPU cores
#BSUB -n 28
## request exclusive access to allocated compute nodes
#BSUB -x
## use only 28-core nodes
#BSUB -R 'span[ptile=28]'
## job name for array job [scenarios]%maximum running simultaneously
#BSUB -J WKNSMSE_A[453-456]%6
## where to save standard output and error streams
###BSUB -oo reports/%J.out
###BSUB -eo reports/%J.err
## save for array jobs
##BSUB -oo /gpfs/afmcefas/simonf/reports/%J.%I.out
##BSUB -eo /gpfs/afmcefas/simonf/reports/%J.%I.err
#BSUB -oo reports/%J.%I.out
#BSUB -eo reports/%J.%I.err
## send email when job starts and finishes
#BSUB -B
#BSUB -N
### start job after previous job finished
##BSUB -w "numended(248796,*)"

. /etc/profile

## remove any loaded modules
module purge
## load R and MPI environment
module load mpi/openmpi/3.1.3/gcc/mellanox  R/3.5.0/gcc-mellanox

echo "got $LSB_DJOB_NUMPROC slots"

### print details about job
echo "This is job $LSB_JOBID and index $LSB_JOBINDEX"
echo "The following ressources have been allocated"
echo $LSB_MCPU_HOSTS

### set working directory
cd $HOME/git/WK_WKNSMSE_cod.27.47d20

# install_local(build=TRUE, path = "mse")

echo "starting the simulations..."
### run MPI job
# mpiexec R CMD BATCH --vanilla --quiet '--args iters=1 years=20 nblocks=1 par_env=2 n_workers=1 HCRoption=1 HCR_comb=6 TAC_constraint=0 BB=0' $HOME/git/WK_WKNSMSE_cod.27.47d20/run_mse.R $HOME/reports/$LSB_JOBID.Rout
### run single job without MPI
#R CMD BATCH --vanilla --quiet '--args iters=1 years=20 nblocks=1 par_env=2 n_workers=1 HCRoption=1 HCR_comb=6 TAC_constraint=0 BB=0' $HOME/git/WK_WKNSMSE_cod.27.47d20/run_mse.R $HOME/reports/$LSB_JOBID.Rout
### run array job
R CMD BATCH --vanilla --quiet "--args iters=1000 years=20 nblocks=200 par_env=2 n_workers=28 HCRoption=1 HCR_comb=$LSB_JOBINDEX TAC_constraint=0 BB=0 OM=0" $HOME/git/WK_WKNSMSE_cod.27.47d20/run_mse.R $HOME/reports/$LSB_JOBID.$LSB_JOBINDEX.Rout

echo "done!"
