#!/bin/sh



rootRef=/Volumes/MARIOLOLO/Data/1KG_Haplotypes/Phase3/
folder=/pasteur/zeus/projets/p02/Hotpaleo/pierre/Projects/RAICES//fineStructure/
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
		~/src/plink1.9/plink --bfile $folder/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP_REDUCED --chr $chr --maf $MAF --geno $GENO --mind $MIND --make-bed --out $folder/shapeITphased/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP_REDUCED.MAF$MAF.GENO$GENO.MIND$MIND.chr$chr
		perl /Users/pierrespc/Documents/PostDoc/scripts/Tools/shapeIT/Prephase/1_runShapeIT.pl $folder/shapeITphased Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP_REDUCED.MAF$MAF.GENO$GENO.MIND$MIND.chr$chr $MIND $GENO $MAF $rootRef/genetic_map_chr$chr"_combined_b37.txt" $rootRef/1000GP_Phase3_chr$chr".hap.gz" $rootRef/1000GP_Phase3_chr$chr".legend.gz" $rootRef/1000GP_Phase3.sample $folder/shapeITphased $state $window $thread $burn $prune $main $effSize
	else
		echo $folder/shapeITphased/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP_REDUCED.MAF$MAF.GENO$GENO.MIND$MIND.chr$chr"_alignedRef_phased.haps" already generated
	fi
done

