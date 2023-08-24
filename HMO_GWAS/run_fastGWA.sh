#!/bin/bash
#SBATCH --job-name=fastgwa
#SBATCH --output=logs/run_GWAS.out
#SBATCH --error=logs/run_GWAS.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=10gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

pheno=$1

timepoints=(W2 M1 M2 M3 M6)

d=/groups/umcg-llnext/tmp01/umcg-dzhernakova/HMO_GWAS/
gcta=/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/v2//tools/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1

geno_file=${d}/genotypes/all_chr.mothers.with_rel.flt
grm=${d}/genotypes/GRM/mothers_flt_pruned_sparse
gcta_mode="--fastGWA-mlm"
    
for tp in ${timepoints[@]}
do 
    mkdir -p ${d}/results/${tp}/  

    pheno_file=${d}/data/HMO_${tp}.txt
    res_file=${d}/results/${tp}/HMO_GWAS_res_${tp}_${pheno}
    
    #get the column number for the HMO in the HMO table
    col=`head -1 ${pheno_file} | sed "s:\t:\n:g" | tail -n+3 | grep -w -n ${pheno} | cut -d ":" -f1`

    $gcta \
      --bfile $geno_file \
      --maf 0.05 \
      --grm-sparse ${grm} \
      ${gcta_mode} \
      --pheno ${pheno_file%txt}noheader.txt \
      --mpheno $col \
      --out $res_file
    
    tmp=${pheno##*HMO_}
    pheno_short=${tmp%_invr}
    awk -v p=$pheno_short -v t=$tp 'BEGIN {FS=OFS="\t"}; {if (NR > 1 && $10 < 5e-8) print t, p, $0}' ${res_file}.fastGWA | sort -k12,12g > ${res_file}.fastGWA.5e-08.txt
    gzip -f ${res_file}.fastGWA 
    
done