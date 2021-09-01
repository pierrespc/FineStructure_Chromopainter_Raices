#!/bin/bash



folder=/pasteur/zeus/projets/p02/Hotpaleo/pierre/Projects/RAICES/

SinguExec="singularity exec --home $HOME:/home/$USER --bind /pasteur "

folderImg=/pasteur/zeus/projets/p02/Hotpaleo/common_data/VMs/singularity/
fsImage=$folderImg/evolbioinfo-finestructure-v4.1.1.img

folderScripts=/pasteur/zeus/projets/p02/Hotpaleo/pierre/Scripts/RAICES/fineStructure/

#commandConvert="$SinguExec $fsImage perl impute2chromopainter.pl"
commandConvert="perl $folderScripts/2b_impute2chromopainter.pl"
commandFs="$SinguExec $fsImage fs"
module load singularity
module load java


###Following exactly Ongaro et al. 2019. First cluster donor pops and then analyze admixed individuals


###step1
#estimated the mutation/emission and the switch rate parameters
#with ten steps of the Expectation–Maximization (E–M) algorithm 
#on a subset of chromosomes {4, 10, 15, 22}
# using every individual except admixed individuals (including Native amerfican with <0.95 of native american ancestry) both as “donor” and “recipients.”


emStage1=10

mkdir $folder/fineStructure/Outputs/
cd $folder/fineStructure/Outputs/

#nR=$(cat $folder/fineStructure//Inputs/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP_REDUCED.MAF0.0000001.GENO0.02.MIND0.05.initialDonorList | wc -l)
#echo $nR
mkdir stage1
cd stage1 
if [ ! -e stage1.Combined ]
then
	for i in 22
	do
		if [ ! -s stage1.chr$i.EMprobs.out ]
		then
			$commandFs cp \
				-g $folder/fineStructure/Inputs/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP_REDUCED.MAF0.0000001.GENO0.02.MIND0.05.chr$i"_alignedRef_phased.phase" \
				-r $folder/fineStructure/Inputs/chr$i.recomb \
				-t $folder/fineStructure/Inputs/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP_REDUCED.MAF0.0000001.GENO0.02.MIND0.05_STEP1.ids \
				-i $emStage1 \
				-a 0 0 \
				-im \
				-in \
				-o stage1.chr$i 
				#--emfilesonly \
				#-s $NR 
		else
			echo stage1.chr$i already generated
		fi
	done
	exit
	###average Ne and m across chromosomes and indviduals
	Rscript 7b_averageStage1.R {22,15,10,4}
else
	echo stage1.Combined already generated
fi

mu=$(awk '{if(NR==2)print $2}' stage1.Combined)
Ne=$(awk '{if(NR==2)print $1}' stage1.Combined)

echo $mu $Ne

cd ..
mkdir stage2
cd stage2
if [ ! -e output.chunkcounts.out ]
then
	for i in {22..1}
	do
		if [ ! -s stage2.chr$i.regionsquaredchunkcounts.out ]
		then
			fs cp \
				-g ../Infiles/SouthAmerica_dataset_QCed.KinshipFilter.AdmFilter.MAF0.0000001.GENO0.02.MIND0.05.chr$i"_alignedRef_phased.phase" \
				-r ../Infiles/chr$i.recomb \
				-t ../Infiles/SouthAmerica_dataset_QCed.KinshipFilter.AdmFilter.MAF0.0000001.GENO0.02.MIND0.05.ids \
				-a 0 0 \
				-k 50 \
				-M $mu \
				-n $Ne \
				-s $nR \
				-o stage2.chr$i
		else
			echo stage2.chr$i already generated
		fi
	done
	/Users/pierrespc/Documents/PostDoc/scripts/Tools/chromopainter_mac/chromocombine  -l -m stage2.chr{1..22}
else
	echo output.chunkcounts.out already generated
fi

cd ..
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
