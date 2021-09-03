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


###Following: Ongaro et al. 2019. First cluster donor pops and then analyze admixed individuals as recipients.


###stage1
#estimated the mutation/emission and the switch rate parameters
#with ten steps of the Expectation–Maximization (E–M) algorithm 
#on a subset of chromosomes {4, 10, 15, 22}
#Using any individual except admixed individuals (including Native amerfican with <0.95 of native american ancestry) both as “donor” and “recipients.”
# For computation reason we randomly selected a subset of 600 donors for stage1


emStage1=10

mkdir $folder/fineStructure/Outputs/
cd $folder/fineStructure/Outputs/

#nR=$(cat $folder/fineStructure//Inputs/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP_REDUCED.MAF0.0000001.GENO0.02.MIND0.05.initialDonorList | wc -l)
#echo $nR
mkdir stage1
cd stage1 
if [ ! -e stage1.Combined ]
then
	for i in 22 18 
	do
		if [ ! -s stage1.chr$i.EMprobs.out ]
		then
			$commandFs cp \
				-g $folder/fineStructure/Inputs/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP_REDUCED.MAF0.0000001.GENO0.02.MIND0.05.chr$i"_alignedRef_phased.phase" \
				-r $folder/fineStructure/Inputs/chr$i.recomb \
				-t $folder/fineStructure/Inputs/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP_REDUCED.MAF0.0000001.GENO0.02.MIND0.05_STEP1.subSet.ids \
				-i $emStage1 \
				-a 0 0 \
				-iM \
				-in \
				-o stage1.chr$i > stage1.chr$i.LOG
				#-s $nR
				#--emfilesonly \
		else
			echo stage1.chr$i already generated
		fi
	done
	###average Ne and m across chromosomes and indviduals
	Rscript $folderScripts/3a_averageStage1.R $folder/fineStructure/Inputs/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP_REDUCED.MAF0.0000001.GENO0.02.MIND0.05.chr {22,18,10,7,3}
else
	echo stage1.Combined already generated
fi

mu=$(awk '{if(NR==2)print $2}' stage1.Combined)
Ne=$(awk '{if(NR==2)print $1}' stage1.Combined)

echo $mu $Ne

cd ..

####now stage 2: (describe by Ongaro et al. 2019): 
#We reconstruct each individual’s chromosomes as a series of genomic fragments inherited (copied) from a set of donor individuals,
# using the information on the allelic state of recipient and donors at each available position. 
#Briefly, we ‘painted’ the genomic profile of each donor as the combination of fragments received from other donor individuals.



mkdir stage2
cd stage2

totalInds=$(awk '{if($3==1)print $0}' $folder/fineStructure/Inputs/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP_REDUCED.MAF0.0000001.GENO0.02.MIND0.05_STEP1.ids | wc -l)
step=20
if [ ! -e output.chunkcounts.out ]
then
	for i in {22..1}
	do
		if [ ! -s stage2.chr$i.regionsquaredchunkcounts.out ]
		then
			mkdir stage2.chr$i.paral/
			mkdir stage2.chr$i.paral/logs/
			ind1=1
			ind2=$step
			while [ $ind1 -le $totalInds ]
			do
				echo $ind1 $ind2
				if [ ! -s stage2.chr$i.paral/$ind1.$ind2.regionsquaredchunkcounts.out ]
				then
					jobS2paral=$(sbatch -J S2.$i.$ind1.$ind2 -o stage2.chr$i.paral/logs/$ind1.$ind2.o -e stage2.chr$i.paral/logs/$ind1.$ind2.e --mem=4G \
						--wrap "$commandFs cp \
							-g $folder/fineStructure/Inputs/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP_REDUCED.MAF0.0000001.GENO0.02.MIND0.05.chr$i\"_alignedRef_phased.phase\" \
							-r $folder/fineStructure/Inputs/chr$i.recomb \
							-t $folder/fineStructure/Inputs/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP_REDUCED.MAF0.0000001.GENO0.02.MIND0.05_STEP1.ids \
							-a $ind1 $ind2 \
							-k 100 \
							-M $mu \
							-n $Ne \
							-o stage2.chr$i.paral/$ind1.$ind2")
				fi
				let ' ind1 = ind1 + step '
				let ' ind2 = ind2 + step '
				
			done
		else
			echo stage2.chr$i already generated
		fi
		exit
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
