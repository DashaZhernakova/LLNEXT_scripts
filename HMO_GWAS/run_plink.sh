#!/bin/bash
#SBATCH --job-name=plink
#SBATCH --output=logs/run_GWAS.out
#SBATCH --error=logs/run_GWAS.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=10gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

ml PLINK

timepoints=(W2 M1 M2 M3 M6)
d=/groups/umcg-llnext/tmp01/umcg-dzhernakova/HMO_GWAS/
geno_file=${d}/genotypes/all_chr.mothers.with_rel.flt
    
for tp in ${timepoints[@]}
do  
    # non-imputed 
    pheno_file=${d}/data/HMO_${tp}.txt
    res_file=${d}/results_plink/${tp}/HMO_GWAS_res_${tp}
    mkdir -p ${d}/results_plink/${tp}
    
    plink2 \
        --glm hide-covar \
        --bfile ${geno_file} \
        --pheno ${pheno_file} \
        --covar ${d}/genotypes/mothers.PC1-2.txt \
        --out ${res_file} \
        --pfilter 0.05
    
    for f in ${res_file}*linear
    do
        tmp=${f##*HMO_}
        pheno_short=${tmp%_invr.*}
        awk -v p=$pheno_short -v t=$tp 'BEGIN {FS=OFS="\t"}; {if (NR > 1 && $12 < 5e-8) print t, p, $0}' ${f} | sort -k14,14g > ${f}.5e-08.txt
        gzip -f $f 
    done
    
    # imputed
    pheno_file=${d}/data/HMO_${tp}_imp.txt
    res_file=${d}/results_plink_imp/${tp}/HMO_GWAS_res_${tp}
    mkdir -p ${d}/results_plink_imp/${tp} 
    
    plink2 \
        --glm hide-covar \
        --bfile ${geno_file} \
        --pheno ${pheno_file} \
        --covar ${d}/genotypes/mothers.PC1-2.txt \
        --out ${res_file} \
        --pfilter 0.05
    
    for f in ${res_file}*linear
    do
        tmp=${f##*HMO_}
        pheno_short=${tmp%_invr.*}
        awk -v p=$pheno_short -v t=$tp 'BEGIN {FS=OFS="\t"}; {if (NR > 1 && $12 < 5e-8) print t, p, $0}' ${f} | sort -k14,14g > ${f}.5e-08.txt
        gzip -f $f 
    done
    
done

