#!/bin/bash


SinguExec="singularity exec --home $HOME:/home/$USER --bind /pasteur "

folderImg=/pasteur/zeus/projets/p02/Hotpaleo/common_data/VMs/singularity/
fsImage=$folderImg/evolbioinfo-finestructure-v4.1.1.img


commandConvert="$SinguExec $fsImage perl impute2chromopainter.pl"
commandFs="$SinguExec $fsImage fs"
module load singularity
module load java


rootRef=/pasteur/zeus/projets/p02/Hotpaleo/common_data/db/ref_datasets/Human/1000GP_Phase3/

folderScripts=/pasteur/zeus/projets/p02/Hotpaleo/pierre/Scripts/RAICES/

folder=/pasteur/zeus/projets/p02/Hotpaleo/pierre/Projects/RAICES/
cd $folder 

mkdir fineStructure
mkdir fineStructure/Inputs/


a=$(ls  shapeITphased/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP_REDUCED.MAF0.0000001.GENO0.02.MIND0.05.chr{1..22}_alignedRef.bim | wc -l)
echo $a
if [ $a -ne 22 ]
then
	a=$(ls  shapeITphased/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP_REDUCED.MAF0.0000001.GENO0.02.MIND0.05.chr{1..22}.bim | wc -l)
	if [ $a -eq 22 ]
	then
		suff=""
	else
		echo "not enough *bim files in $folder/shapeITphased/"
		answer="?"
		while [ $answer != "y" ] && [ $answer != "n" ]
		do 
			echo "copy $folder/shapeITphased/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP_REDUCED.MAF0.0000001.GENO0.02.MIND0.05.chr<CHR>.bim to $folder/shapeITphased/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP_REDUCED.MAF0.0000001.GENO0.02.MIND0.05.chr<CHR>_alignedRef.bim????"
			echo "type y/n"
			read answer
		done
		if [ $answer == "n" ]
		then
			echo "well fix this then"
			exit
		else
			suff="_alignedRef"
			for chr in {1..22}
			do
				if [ ! -e shapeITphased/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP_REDUCED.MAF0.0000001.GENO0.02.MIND0.05.chr$chr"_alignedRef.bim" ]
				then
					echo coyping for chr$chr
					cp shapeITphased/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP_REDUCED.MAF0.0000001.GENO0.02.MIND0.05.chr$chr.bim shapeITphased/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP_REDUCED.MAF0.0000001.GENO0.02.MIND0.05.chr$chr"_alignedRef.bim"
				else
					echo already there for chr$chr
				fi
			done
		fi
	fi
else
	suff="_alignedRef"
fi
echo "we have fixed bim files to :$folder/shapeITphased/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP_REDUCED.MAF0.0000001.GENO0.02.MIND0.05.chr<CHR>$suff.bim"
a=$(ls  shapeITphased/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP_REDUCED.MAF0.0000001.GENO0.02.MIND0.05.chr{1..22}$suff.bim | wc -l)
if [ $a != 22 ]
then
	echo "still not enough *bim files in $folder/shapeITphased/"
	exit
fi
a=$(ls  shapeITphased/withCM.Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP_REDUCED.MAF0.0000001.GENO0.02.MIND0.05.chr{1..22}$suff.bim | wc -l)
if [ $a != 22 ]
then
	Rscript $folderScripts/2a_putGenMap_inMapFile.R \
		T \
		$suff \
		shapeITphased/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP_REDUCED.MAF0.0000001.GENO0.02.MIND0.05.chr \
		T \
		T \
		$rootRef/genetic_map_chr \
		6 \
		shapeITphased/withCM.Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP_REDUCED.MAF0.0000001.GENO0.02.MIND0.05.chr
else
	echo cM ok!
fi
exit
for chr in {1..22}
do
	echo $chr
	echo "start.pos recom.rate.perbp" > fineStructure/chr$chr.recomb
	awk 'BEGIN{prevCM=0;prevBP=0} {rate=($3-prevCM)/($4-prevBP); print $4,rate; prevCM=$3; prevBP=$4}' shapeITphased/withCM.Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP_REDUCED.MAF0.0000001.GENO0.02.MIND0.05.chr$chr"_alignedRef.bim" >> fineStructure/Inputs/chr$chr.recomb
	if [  ! -s fineStructure/Inputs/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP_REDUCED.MAF0.0000001.GENO0.02.MIND0.05.chr$chr"_alignedRef_phased.phase" ]
        then
                $commandConvert \
			shapeITphased/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP_REDUCED.MAF0.0000001.GENO0.02.MIND0.05.chr$chr"_alignedRef_phased.haps" \
			fineStructure/Inputs/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP_REDUCED.MAF0.0000001.GENO0.02.MIND0.05.chr$chr"_alignedRef_phased"
        else
                echo fineStructure/Inputs/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP_REDUCED.MAF0.0000001.GENO0.02.MIND0.05.chr$chr"_alignedRef_phased.phase" already generated
        fi


done

awk '{if(NR>2)print $2,$1,1}' shapeITphased/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP_REDUCED.MAF0.0000001.GENO0.02.MIND0.05.chr22_alignedRef_phased.sample > fineStructure/Inputs/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP_REDUCED.MAF0.0000001.GENO0.02.MIND0.05.ids 

test=$(awk '{print $2}' fineStructure/Inputs/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP_REDUCED.MAF0.0000001.GENO0.02.MIND0.05.ids | uniq | wc -l)
test2=$(awk '{print $2}' fineStructure/Inputs/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP_REDUCED.MAF0.0000001.GENO0.02.MIND0.05.ids | sort | uniq | wc -l)
if [ $test != $test2 ]
then
	echo "your .phase file must be orderd by individuals... need to fix this"
	exit
fi
awk '{print $2}' fineStructure/Inputs/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP_REDUCED.MAF0.0000001.GENO0.02.MIND0.05.ids | uniq -c | awk '{print $2,$1}' > fineStructure/Inputs/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP_REDUCED.MAF0.0000001.GENO0.02.MIND0.05.initialDonorList


