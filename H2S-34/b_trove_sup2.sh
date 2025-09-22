#!/bin/bash -l
export exec=j-xsec_2211_i17.x


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

export xname=`echo $name | sed -e 's/\_super2//'`
echo $xname


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



qsub -A dp060 -S /bin/tcsh -N $name -o $name.o -e $name.e  -l "walltime=$wclim:00:00,nodes=1:ppn=$nproc"   \
     -v "name=$name,xname=$xname,pwd=$pwd,nproc=$nprocs,exec=$exec"  \
     $pwd/run_trove_sup2.csh 

#    --workdir=$pwd --hint=compute_bound --no-requeue --mem=$MEM -p skylake \
     
     
     
