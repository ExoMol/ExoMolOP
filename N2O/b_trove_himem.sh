        #!/bin/bash -l
#
# Check if our working directory is on the central file server
#

export exec=j-trove_3009_NH3_Roman.x 


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
export MEM=1400gb
export wclim=$2


echo "Nproc=" $nproc



echo "Nnodes=" 1, "Nproc=" $nproc, " Memory = "$MEM, "jobtype = " $jobtype, "wclimit = " $wclim
echo "Working dir is " $pwd

qsub -A dp060  -N $name -o $name.o -e $name.e  -l "walltime=$wclim:00:00,nodes=1:ppn=$nproc"   \
     -q highmem \
     -v "name=$name,pwd=$pwd,nproc=$nproc,exec=$exec"  \
     $pwd/run_trove.sh 
     

