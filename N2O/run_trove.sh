#!/bin/bash -l
#!

export wdir=$pwd

echo $name
echo $nproc
#!
export OMP_NUM_THREADS=$nproc
export MKL_NUM_THREADS=$nproc
#!
export KMP_LIBRARY=turnaround
export KMP_AFFINITY=disabled

export KMP_STACKSIZE=1gb
export OMP_MAX_ACTIVE_LEVELS=1
#!
hostname
#!
#limit
#limit datasize unlimited
#limit
#!
cd   $wdir
echo -e "Changed directory to `pwd`.\n"

export LAUNCH=time  ###"dplace -x2"
export TMPDIR=$wdir
#!
echo "TMPDIR = " $TMPDIR
echo "USER = " $USER
echo "OMP_NUM_THREADS = " $OMP_NUM_THREADS
echo "wdir" $wdir
#!
##export wdir  $PBS_O_WORKDIR
echo "wdir" $wdir
echo "OMP_NUM_THREADS=" $OMP_NUM_THREADS
#!
cd $wdir
#!
echo $wdir
#!

if [ -e "$name.out" ]; then
   /bin/rm $name.out
fi


export JOBID=$SLURM_JOB_ID

echo -e "JobID: $JOBID\n======"
echo "Time: `date`"
echo "Running on master node: `hostname`"
echo "Current directory: `pwd`"
#!
$LAUNCH $pwd/$exec < $pwd/$name.inp > $pwd/$name.out


#!
echo "DONE"
