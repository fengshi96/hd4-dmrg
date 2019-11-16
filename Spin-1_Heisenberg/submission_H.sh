#!/bin/bash -x
zero=0.0; pi=3.14159265
JJ=1.0;


for Lambda in 0.00 #0.05 0.10 0.15 0.20 0.25 0.30 0.32 0.33 0.34 0.35 0.36 0.37 0.38 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95 1.00 1.05 1.10 1.12 1.14 1.15 1.16 1.18 1.19 1.20 1.21 1.22 1.23 1.24 1.25 1.30 1.35 1.40 1.45 1.50
 
do
  	
        mkdir -p Lambda_$Lambda
        cd Lambda_$Lambda

	cp ../Midinput.inp .

        sed -e 's/LLLLLL/'$Lambda'/g' \
            -e 's/JJJJJJ/'$JJ'/g' \
                <Midinput.inp  > input.inp
        rm Midinput.inp

mydir=$(pwd)
rawjob="$(cat <<EOF
#!/bin/sh
#PBS -N L-Spin-1-H-$Lambda
#PBS -l nodes=1:ppn=8 
#PBS -l walltime=64:00:00 
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
../dmrg input.inp &> runForinput.cout
time
EOF
)"

	echo "$rawjob" &> job.pbs
        #####qsub -W depend=afterok:181255 $sub_out
	qsub job.pbs


        echo -n "Lambda = $Lambda,    "

cd ..
done




