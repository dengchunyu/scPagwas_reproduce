#!/bin/sh
#PBS -N ad
#PBS -q workq
#PBS -j oe

cd /share2/pub/zhenggw/zhenggw/ldsc_AD/
source activate ldsc
for FG in ./AD_ldscores/*GeneSet
do
	for FC in ./1000G_EUR_Phase3_plink/*bim
		do
		index=`echo $FC|cut -d . -f 5`
		python ./ldsc/make_annot.py --gene-set-file $FG --gene-coord-file ENSG_coord.txt --windowsize 100000 --bimfile $FC --annot-file ./AD_ldscores/${FG##*\/}.$index.annot.gz
	done
done

#!/bin/sh
#PBS -N ad
#PBS -q workq
#PBS -j oe

cd /share2/pub/zhenggw/zhenggw/ldsc_AD/
source activate ldsc

for FG in ./AD_ldscores/*GeneSet
  do
    for FC in ./1000G_EUR_Phase3_plink/*bim
      do 
        index=`echo $FC|cut -d . -f 5`
        python ./ldsc/ldsc.py --l2 --bfile ./1000G_EUR_Phase3_plink/1000G.EUR.QC.$index --ld-wind-cm 1 --annot ./AD_ldscores/${FG##*\/}.$index.annot.gz --thin-annot --out ./AD_ldscores/${FG##*\/}.$index --print-snps ./hapmap3_snps/hm.$index.snp
       done 
done


#!/bin/sh
#PBS -N ct
#PBS -q workq
#PBS -j oe

cd /share2/pub/zhenggw/zhenggw/ldsc_AD/
source activate ldsc
for FG in ./CT_ldscores/*GeneSet
do
	for FC in ./1000G_EUR_Phase3_plink/*bim
		do
		index=`echo $FC|cut -d . -f 5`
		python ./ldsc/make_annot.py --gene-set-file $FG --gene-coord-file ENSG_coord.txt --windowsize 100000 --bimfile $FC --annot-file ./CT_ldscores/${FG##*\/}.$index.annot.gz
	done
done


#!/bin/sh
#PBS -N ct2
#PBS -q workq
#PBS -j oe

cd /share2/pub/zhenggw/zhenggw/ldsc_AD/
source activate ldsc

for FG in ./CT_ldscores/*GeneSet
  do
    for FC in ./1000G_EUR_Phase3_plink/*bim
      do 
        index=`echo $FC|cut -d . -f 5`
        python ./ldsc/ldsc.py --l2 --bfile ./1000G_EUR_Phase3_plink/1000G.EUR.QC.$index --ld-wind-cm 1 --annot ./CT_ldscores/${FG##*\/}.$index.annot.gz --thin-annot --out ./CT_ldscores/${FG##*\/}.$index --print-snps ./hapmap3_snps/hm.$index.snp
       done 
done


