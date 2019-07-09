#!/bin/sh
#$ -e /exports/eddie/scratch/cmclean5
#$ -o /exports/eddie/scratch/cmclean5

SUBDIR=$JOB_ID
echo "SUBDIR is $SUBDIR"

EXECDIR=/exports/home/cmclean5/STUDIES/HumanSynaptosome/localization
SCRIPTDIR=/exports/home/cmclean5/STUDIES/HumanSynaptosome/localization/EDDIE/SCRIPTS
DATADIR=/exports/home/cmclean5/STUDIES/HumanSynaptosome/localization/datasets

SEED=$SUBDIR
PERMS=10

cd $TMPDIR
echo "WORKING DIR " $TMPDIR

cp -r $SCRIPTDIR/calEDDIEDiseasePairs.R .
cp -r $DATADIR/PPI_Presynaptic_Published.gml .
#cp -r $DATADIR/PPI_PSP_clean_Published.gml .
#cp -r $DATADIR/PPI_PSP_reduced.gml .

# initiallise environment module
. /etc/profile.d/modules.sh

# load module R
module load R 

Rscript calEDDIEDiseasePairs.R $SEED $PERMS

OUTDIR=$EXECDIR/EDDIE/RESULTS/$SUBDIR

if [ ! -d $OUTDIR ]; then
    mkdir $OUTDIR
fi

cp -v *.csv $OUTDIR

echo "$0 done!"

exit 0
