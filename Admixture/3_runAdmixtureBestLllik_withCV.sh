#!/bin/bash


outFolder=/pasteur/zeus/projets/p02/Hotpaleo/pierre/Projects/RAICES/
outPref=Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP_REDUCED.Filtered.pruned

mkdir $outFolder/Admixture/BestRUNperK/

cd $outFolder/Admixture/BestRUNperK/

ls 
SinguExec="singularity exec --home $HOME:/home/$USER --bind /pasteur "

folderImg=/pasteur/zeus/projets/p02/Hotpaleo/common_data/VMs/singularity/
admImage=$folderImg/evolbioinfo-admixture-v1.3.0.img
commandAdm="$SinguExec $admImage admixture"

for K in {2..13}
do
	
	rep=$(awk -v k=$K '{if($1==k)print $2}' ListBESTruns)	
	echo $rep
	if [[ ! -e $outPref.rep$rep.${K}.Q ]]
	then
		$commandAdm $outFolder/StartFilteredData/$outPref.bed $K --cv --seed $rep | tee $outPref.rep$rep.K${K}.out
		mv $outPref.${K}.Q $outPref.rep$rep.${K}.Q
		mv $outPref.${K}.P $outPref.rep$rep.${K}.P
        else 
       	        echo $outPref.rep$rep.K${K}.out already exists!
        fi
done

