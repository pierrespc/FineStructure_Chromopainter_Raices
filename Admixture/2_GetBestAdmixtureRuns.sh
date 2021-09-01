#!/bin/bash

folder=/pasteur/zeus/projets/p02/Hotpaleo/pierre/Projects/RAICES/
prefBED=Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP_REDUCED.Filtered.pruned
#prefED=FinalSetPierre.1KG.GENO.MIND.MAF
rootOut=$folder/Admixture/
mkdir $rootOut/BestRUNperK/
cd $rootOut

echo K RUN log CVscore > BestRUNperK/ListBESTruns
for k in {2..13}
do
	rm BestRUNperK/K$k.LogLik
	for run in {1..10}
	do 
		log=$(grep Loglik K$k/$prefBED.rep$run.K$k.out | grep -v delta | awk '{print $2}')
		echo $run $log >> BestRUNperK/K$k.LogLik
	done
	RUNbest=$(sort -k2 -n BestRUNperK/K$k.LogLik | tail -1 | awk '{print $1}')
	Logbest=$(sort -k2 -n BestRUNperK/K$k.LogLik | tail -1 | awk '{print $2}')
	#CVbest=$(grep "CV error" RUN$RUNbest/$prefBED.pruned.K$k.out | awk '{print $4}')

	#echo $k $RUNbest $Logbest $CVbest >> BestRUNperK/ListBESTruns
	echo $k $RUNbest $Logbest >> BestRUNperK/ListBESTruns

	#cp RUN$RUNbest/$prefBED.pruned.$k.P BestRUNperK/
	#cp RUN$RUNbest/$prefBED.pruned.$k.Q BestRUNperK/
	#cp RUN$RUNbest/$prefBED.pruned.K$k.out BestRUNperK/
done
	
	
		

