#!/bin/bash



folder=/pasteur/zeus/projets/p02/Hotpaleo/pierre/Projects/RAICES/

SinguExec="singularity exec --home $HOME:/home/$USER --bind /pasteur "

folderImg=/pasteur/zeus/projets/p02/Hotpaleo/common_data/VMs/singularity/
fsImage=$folderImg/evolbioinfo-finestructure-v4.1.1.img

folderScripts=/pasteur/zeus/projets/p02/Hotpaleo/pierre/Scripts/RAICES/fineStructure/

commandConvert="$SinguExec $fsImage impute2chromopainter.pl"
#commandConvert="perl $folderScripts/2b_impute2chromopainter.pl"
commandFs="$SinguExec $fsImage fs"
commandCombine="$SinguExec $fsImage chromocombine"

module load singularity
module load java


###Following: Ongaro et al. 2019. First cluster donor pops and then analyze admixed individuals as recipients.


###stage1
#estimated the mutation/emission and the switch rate parameters
#with ten steps of the Expectation–Maximization (E–M) algorithm 
#on a subset of chromosomes {4, 10, 15, 22}
#Using any individual except admixed individuals (including Native amerfican with <0.95 of native american ancestry) both as “donor” and “recipients.”
# For computation reason we randomly selected a subset of 600 donors for stage1



cd $folder/fineStructure/Outputs/stage1/
mu=$(awk '{if(NR==2)print $2}' stage1.Combined)
Ne=$(awk '{if(NR==2)print $1}' stage1.Combined)

mkdir stage3
cd stage3
if [ ! -e stage3.mcmc.xml ]
then
	fs fs \
		-m oMCMC \
		-x 1000000 \
		-y 2000000 \
		-z 10000 \
		../stage2/output.chunkcounts.out \
		stage3.mcmc.xml
else
	echo stage3.mcmc.xml already generated
fi
#continue the previous run for 100K additional steps, treating the original run as burnin.
if [ ! -e stage3.mcmc.longer.xml ]
then
	fs fs \
		-x 0 \
		-y 100000 \
		-z 10000 \
		../stage2/output.chunkcounts.out \
		stage3.mcmc.xml \
		stage3.mcmc.longer.xml
                       
else
	echo stage3.mcmc.longer.xml already generated
fi

## Infers a tree, using the best state seen in out.mcmc.xml as the initial state.
if [ ! -e stage3.tree.xml ]
then
	fs fs \
		-m T \
		-x 0 \
		../stage2/output.chunkcounts.out \
		stage3.mcmc.longer.xml \
		stage3.tree.xml
else
	echo stage3.tree.xml already generated
fi


#Infers a tree, using (-T 1) the maximum concordance state over out.mcmc.xml as the initial state. This is reported with full likelihood ordering (-k 2), useful for cutting at a given number of ppulations K (but may look bad in the GUI).
if [ ! -e stage3.tree.Other.xml ]
then
	fs fs \
		-m T \
		-k 2 \
		-T 1 \
		../stage2/output.chunkcounts.out \
		stage3.mcmc.longer.xml \
		stage3.tree.Other.xml
fi
