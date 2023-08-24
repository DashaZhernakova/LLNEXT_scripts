d=/groups/umcg-llnext/tmp01/umcg-dzhernakova/HMO_GWAS/
cd ${d}

#
# 1. Run GWAS using fastGWA
#
head -1 ${d}/data/HMO_M3.txt | sed "s:\t:\n:g" | tail -n+3 > ${d}/data/all_hmos.txt

cd ${d}/scripts/

# Submit jobs for CLR-transformed, non-zero CLR-transformed and binary tables
while read line
do 
    echo $line; 
    sbatch -o logs/imp_${line}.out -e logs/imp_${line}.err -J $line run_fastGWA.sh $line
done < ../data/all_hmos.txt 
#
# 2. Combine results
#

cd ${d}/results
sort -m -k12,12g -T $TMPDIR -S 10G --buffer-size=1000 */*fastGWA.5e-08.txt > HMO_GWAS_non_imputed.5e-08.txt

cd ${d}/results_imp
sort -m -k12,12g -T $TMPDIR -S 10G --buffer-size=1000 */*fastGWA.5e-08.txt > HMO_GWAS_imputed.5e-08.txt
