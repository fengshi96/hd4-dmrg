#!/bin/bash -x
zero=0.0; pi=3.14159265

#string="-m64 -std=c++11 -w -fPIC -I. -I/home/patel.3537/0.Codes/itensor -I/home/patel.3537/0.Codes/ -I/opt/intel/mkl/include  -O3 -DNDEBUG -Wall SpinCurr.cc -o SpinCurr -L/home/patel.3537/0.Codes/itensor/lib -litensor -L/opt/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_rt -lmkl_core -liomp5 -lpthread"
#echo "g++ $string"
#g++ $string

cd ..
for Hz in 5.05 5.10 5.15 5.20 5.25 5.30 5.35 5.40 5.45 5.50 5.55 5.60 5.65 5.70 5.75 5.80 5.85 5.90 5.95 6.00 #0.00 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95 1.00 1.05 1.10 1.15 1.20 1.25 1.30 1.35 1.40 1.45 1.50 1.55 1.60 1.65 1.70 1.75 1.80 1.85 1.90 1.95 2.00 2.05 2.10 2.15 2.20 2.25 2.30 2.35 2.40 2.45 2.50 2.55 2.60 2.65 2.70 2.75 2.80 2.85 2.90 2.95 3.00 3.05 3.10 3.15 3.20 3.25 3.30 3.35 3.40 3.45 3.50 3.55 3.60 3.65 3.70 3.75 3.80 3.85 3.90 3.95 4.00 4.05 4.10 4.15 4.20 4.25 4.30 4.35 4.40 4.45 4.50 4.55 4.60 4.65 4.70 4.75 4.80 4.85 4.90 4.95 5.00
do
	cd Hz_$Hz
	mkdir -p Observables
	cd Observables
	
rawjob="$(cat <<EOF
#!/bin/sh
#PBS -N LL-$Hz
#PBS -l nodes=1:ppn=8
#PBS -l walltime=24:00:00
#PBS -l vmem=64gb
#PBS -l mem=64gb
#PBS -m a
#PBS -M feng.934@osu.edu
#PBS -j oe
#PBS -q batch
hostname
#PBS -r n
module load intel
module load mkl
cd \$PBS_O_WORKDIR
time
../../Observables/SpinCurr ../input.inp &> Observables_Local2.dat
time
EOF
)"
echo "$rawjob" &> obs_Local.pbs

	 qsub obs_Local.pbs
	echo "done with $Hz" 
	cd ../../
done




