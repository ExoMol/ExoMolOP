#!/bin/csh
#!

setenv wdir $pwd

module load intel/compilers/18
module load python/intel/3.6

echo $name
echo $nproc
#!
setenv OMP_NUM_THREADS $nproc
setenv MKL_NUM_THREADS $nproc
#!
setenv KMP_LIBRARY turnaround
setenv KMP_STACKSIZE 1gb
setenv OMP_NESTED FALSE
#setenv ppn $(uniq -c "$PBS_NODEFILE" | head --lines=1 | sed -e 's/^ *\([0-9]\+\) .*$/\1/g')
#!
hostname
#!
limit
limit datasize unlimited
limit
#!
cd   $wdir
echo -e "Changed directory to `pwd`.\n"

setenv LAUNCH time  ###"dplace -x2"
setenv TMPDIR $wdir
#!
echo "TMPDIR = " $TMPDIR
echo "USER = " $USER
echo "OMP_NUM_THREADS = " $OMP_NUM_THREADS
echo "wdir" $wdir
#!
##setenv wdir  $PBS_O_WORKDIR
echo "wdir" $wdir
echo "OMP_NUM_THREADS=" $OMP_NUM_THREADS
#!
cd $wdir
#!
echo $wdir
#!
if (-e $name.out) then
    /bin/rm $name.out
endif


setenv LAUNCH "time numactl --interleave=all"

setenv TMPDIR $wdir
echo "TMPDIR = " $TMPDIR
echo "USER = " $USER
echo "OMP_NUM_THREADS = " $OMP_NUM_THREADS
echo "OMP_NUM_THREADS = " $MKL_NUM_THREADS
echo "PBS_O_WORKDIR" $PBS_O_WORKDIR


setenv JOBID $PBS_JOBID

echo -e "JobID: $JOBID\n======"
echo "Time: `date`"
echo "Running on master node: `hostname`"
echo "Current directory: `pwd`"
#!
$LAUNCH /cm/shared/apps/python/intelpython3/bin/python3 $pwd/$exec --dict $pwd/$name --wl 0.3,50 --resolution 1000 --ngauss 20 --ncores 36 --key_iso_ll rep1 --mname rep2 --mol_mass rep3
#!
echo "DONE"
