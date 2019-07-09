#!/bin/sh

echo "Running on Eddie..."

#load module R
#module load R

WORKINGDIR=/exports/home/cmclean5/STUDIES/HumanSynaptosome/localization/EDDIE

EXE=$WORKINGDIR/SCRIPTS/execute.sh

START=1
END=500

name="Human_Synaptosome_Dis_Loc"

chmod +x $EXE

for i in `seq $START $END`
    do
     qsub -N $name -l h_rt=48:00:00 -l h_vmem=8G $EXE
    done

echo "$0 done!"

exit 0
    
