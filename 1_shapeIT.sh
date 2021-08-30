#!/bin/bash



SinguExec="singularity exec --home $HOME:/home/$USER --bind /pasteur "

folderImg=/pasteur/zeus/projets/p02/Hotpaleo/common_data/VMs/singularity/
plinkImage=$folderImg/biocontainers-plink1.9-v1.90b6.6-181012-1-deb_cv1.img
shapeitImage=$folderImg/xiaoli-shapeit-latest.img 

commandPlink="$SinguExec $plinkImage plink1.9"
commandShapeit="$SinguExec $shapeitImage shapeit"
module load singularity
module load java


rootRef=/Volumes/MARIOLOLO/Data/1KG_Haplotypes/Phase3/

folderScript=/pasteur/zeus/projets/p02/Hotpaleo/pierre/Scripts/RAICES/
folder=/pasteur/zeus/projets/p02/Hotpaleo/pierre/Projects/RAICES/
folderIN=$folder/StartFilteredData/

mkdir $folder
cd $folder
mkdir shapeITphased
mkdir shapeITphased/Log
cd $folder/shapeITphased/Log
###parameters for QC filtering in plink
MIND=0.05
GENO=0.02
MAF=0.0000001
######

###Algorithm and Model parameters for shapeit
state=100
window=2
thread=1
burn=7
prune=8
main=20
effSize=15000
########


for chr in 22
#for chr in {22..1}
do
	if [ ! -e $folder/shapeITphased/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP_REDUCED.MAF$MAF.GENO$GENO.MIND$MIND.chr$chr"_alignedRef_phased.haps" ]
	then
		$commandPlink plink1.9 --bfile $folderIN/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP_REDUCED --chr $chr --maf $MAF --geno $GENO --mind $MIND --make-bed --out $folder/shapeITphased/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP_REDUCED.MAF$MAF.GENO$GENO.MIND$MIND.chr$chr
		perl $folderScript/1_runShapeIT.pl $folder/shapeITphased Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP_REDUCED.MAF$MAF.GENO$GENO.MIND$MIND.chr$chr $MIND $GENO $MAF $rootRef/genetic_map_chr$chr"_combined_b37.txt" $rootRef/1000GP_Phase3_chr$chr".hap.gz" $rootRef/1000GP_Phase3_chr$chr".legend.gz" $rootRef/1000GP_Phase3.sample $folder/shapeITphased $state $window $thread $burn $prune $main $effSize "$commandPlink" "$commandShapeit"
	else
		echo $folder/shapeITphased/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP_REDUCED.MAF$MAF.GENO$GENO.MIND$MIND.chr$chr"_alignedRef_phased.haps" already generated
	fi
done

