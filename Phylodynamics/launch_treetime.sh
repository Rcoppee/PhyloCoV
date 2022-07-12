#! /bin/bash

#Script Romain Coppee
#Creation data: 02/10/2021
#Last modification: 04/10/2021

#FILES contains the BAM file for each sample
FILES_TREE=*.nwk

#initilization - DO NOT TOUCH
REPLICAT=1

#assuming that the alignments (.msa) and trees (.nwk) have the following name:
#2022_05_13-4Rre1HyiYx-1-selection
#2022_05_13-4Rre1HyiYx-2-selection
#...
for f in `ls $FILES_TREE`
do
	MY_TREE="2022_05_13-4Rre1HyiYx-${REPLICAT}-selection.nwk"
	MY_ALI="2022_05_13-4Rre1HyiYx-${REPLICAT}-selection.msa"
	MY_DATES="${REPLICAT}_date.txt"
	MY_DIRECT="replicat${REPLICAT}"
	#we fix the clock rate at 0.006 for dating the tree
	treetime --tree $MY_TREE \
             --aln $MY_ALI \
             --dates $MY_DATES \
             --clock-rate 0.006 \
             --reroot best \
             --no-tip-labels \
             --outdir $MY_DIRECT
	#we rename the trees for subsequent analyses
    mv ${MY_DIRECT}/timetree.nexus ${MY_DIRECT}/${REPLICAT}_tree_dated.nexus
    REPLICAT=$((REPLICAT+1))
done
