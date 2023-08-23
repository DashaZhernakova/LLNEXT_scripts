#!/bin/bash
#SBATCH --job-name=SV
#SBATCH --output=clump.out
#SBATCH --error=clump.err
#SBATCH --time=10:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=10gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

ml PLINK/1.9-beta6-20190617

d=/groups/umcg-llnext/tmp01/umcg-dzhernakova/HMO_GWAS/

timepoints=("M1" "M2" "M3" "M6" "W2")
for t in ${timepoints[@]}
do
    cd ${d}/EMP_HMO_${t}
    echo "Starting with $d"
    hmos=( $( cut -f6 EMP_HMO_${t}.5e-08.genes.gwas_catalog.txt | tail -n+2 | sort | uniq ) )
    mkdir clumped
    for hmo in ${hmos[@]}
    do
     awk -v h=$hmo '{OFS=FS="\t"}; {if ($6 == h || NR == 1) print}' EMP_HMO_${t}.5e-08.genes.gwas_catalog.txt > clumped/tmp.${hmo}.eQTLs.1e-5.txt
     plink --bfile /groups/umcg-llnext/tmp01/umcg-dzhernakova/genotypes/1000G_snps_nodup \
     --clump clumped/tmp.${hmo}.eQTLs.1e-5.txt \
     --clump-snp-field rs_id \
     --clump-field PValue \
     --out clumped/${hmo}.1e-5 \
     --clump-r2 0.4
     
     awk -v h="${hmo}" '{if ($3 != "") print $3 "\t" h}' clumped/${hmo}.1e-5.clumped | tail -n+2 > clumped/${hmo}.1e-5.clumped.SNP-HMO.txt
     
     
    done
    cd ${d}
done