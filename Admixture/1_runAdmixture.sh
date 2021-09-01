#!/bin/bash


outFolder=/pasteur/zeus/projets/p02/Hotpaleo/pierre/Projects/RAICES/
outPref=Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP_REDUCED.Filtered.pruned

SinguExec="singularity exec --home $HOME:/home/$USER --bind /pasteur "

folderImg=/pasteur/zeus/projets/p02/Hotpaleo/common_data/VMs/singularity/
admImage=$folderImg/evolbioinfo-admixture-v1.3.0.img
commandAdm="$SinguExec $admImage admixture"

mkdir $outFolder/Admixture

cd $outFolder/Admixture
for K in {2..13}
do
	cd $outFolder/Admixture
	mkdir K$K
	cd K$K	
	for rep in {1..10}
	do

		if [[ ! -e $outPref.rep$rep.${K}.Q ]]
        	then
		        $commandAdm $outFolder/StartFilteredData/$outPref.bed $K --seed $rep | tee $outPref.rep$rep.K${K}.out
			mv $outPref.${K}.Q $outPref.rep$rep.${K}.Q
			mv $outPref.${K}.P $outPref.rep$rep.${K}.P
	        else 
        	        echo $outPref.rep$rep.K${K}.out already exists!
	        fi
	done
done

