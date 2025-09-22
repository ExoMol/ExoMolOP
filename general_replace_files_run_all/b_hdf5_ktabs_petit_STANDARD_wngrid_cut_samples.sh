#!/bin/bash -l
export exec=create_ktables_hdf5_petit_STANDARD_wngrid_samples.py


export pwd=`pwd`
echo $pwd

export name=`echo $1 | sed -e 's/\.inp//'`
echo $name

export JOB=$3
echo $JOB

if [ -e "$name.o" ]; then
   /bin/rm $name.o
fi

if [ -e "$name.e" ]; then
   /bin/rm $name.e
fi

if [ -e "$name.out" ]; then
  if [ -e "$name.tmp" ]; then
    /bin/rm $name.tmp
  fi
  /bin/mv $name.out $name.tmp
fi

export nproc=36


export jobtype="skylake"
export MEM=175gb
export wclim=$2


echo "Nproc=" $nproc



echo "Nnodes=" 1, "Nproc=" $nproc, " Memory = "$MEM, "jobtype = " $jobtype, "wclimit = " $wclim
echo "Working dir is " $pwd



if [ ! -f $pwd/done/$name.done  ]; then 
qsub -A dp060 -S /bin/tcsh  -N $exec -o $exec.o -e $exec.e  -l "walltime=$wclim:00:00,nodes=1:ppn=$nproc"   \
     -v "name=$name,pwd=$pwd,nproc=$nprocs,exec=$exec"  \
     $pwd/run_hdf5_ktabs_petit_STANDARD_wngrid_cut_samples.csh
else
 echo "skip file, already done."
fi  

#    --workdir=$pwd --hint=compute_bound --no-requeue --mem=$MEM -p skylake \
     
     
     
