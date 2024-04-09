d=/groups/umcg-llnext/tmp01/umcg-dzhernakova/HMO_GWAS/
cd ${d}
script_dir=/groups/umcg-llnext/tmp01/umcg-dzhernakova/LLNEXT_scripts/HMO_GWAS/

#
# 1. Run GWAS using plink
#
head -1 ${d}/data/HMO_M3.txt | sed "s:\t:\n:g" | tail -n+3 > ${d}/data/all_hmos.txt

cd ${d}/scripts/

# Submit jobs 
while read line
do 
    echo $line; 
    sbatch -o logs/${line}.out -e logs/${line}.err -J $line run_plink.sh $line
done < ../data/all_hmos.txt 

#
# 2. Combine results
#

cd ${d}/results_plink

# Add rs ids and clump
conda activate py39
for f in */*.5e-08.txt
do
 if [ -s $f ]; then
     python ${script_dir}/utility_scripts/add_rs_by_position.py \
     $f \
     /groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/resources/All_20180423.vcf.gz \
     4 0 \
     > ${f}.rsids.txt
 
     Rscript ${script_dir}/utility_scripts/clump.R ${f}.rsids.txt
  fi 
done

# Combine per-HMO files
sort -m -k15,15g -T $TMPDIR -S 10G --buffer-size=1000 */*5e-08.txt.rsids.txt  \
> HMO_GWAS_plink_nonimputed.5e-08.no_clumping.txt

sort -m -k15,15g -T $TMPDIR -S 10G --buffer-size=1000 */*.rsids.txt.clumped_0.1.txt | \
grep -v -w pval > HMO_GWAS_plink_nonimputed.5e-08.clumped.txt


#
# 3. Annotate
#

col=3
genes=${d}/../resources/ensembl_b37_genes.bed
genes_prot=${d}/../resources/ensembl_b37_protcoding_genes.bed

module load BEDTools
genes=/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS//resources/ensembl_b37_genes.bed
genes_prot=/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS//resources/ensembl_b37_protcoding_genes.bed

f=HMO_GWAS_plink_nonimputed.5e-08.txt
col=6

# Add overlapping genes
awk -v c=${col} '{FS=OFS="\t"}; {split($c, snp, ":"); print snp[1], snp[2] - 1, snp[2], $c}' $f > ${f}.tmp.bed
echo -e "SNP\tOverlapping_genes" > ${f}.tmp.genes
intersectBed -a ${f}.tmp.bed -b $genes -wa -wb | \
  python2.7 ${script_dir}/utility_scripts/collapseIntersectBedRes.py - \
  >> ${f}.tmp.genes

 python3 ${script_dir}/add_columns_from_file_v2.py  \
  -i ${f} -f ${f}.tmp.genes -i_m 5 -f_m 0 -f_cols 1 \
  > ${f}.tmp.genes1.txt

# Add protein-coding genes within a 250kb window from the SNP
echo -e "SNP\tProtcoding_genes_within_250kb" > ${f}.tmp.genes2
windowBed -a ${f}.tmp.bed -b $genes_prot -w 250000 | \
python2.7 ${script_dir}/utility_scripts/collapseIntersectBedRes.py - \
  >> ${f}.tmp.genes2

python3 ${script_dir}/add_columns_from_file_v2.py  \
  -i ${f}.tmp.genes1.txt -f ${f}.tmp.genes2 -i_m 5 -f_m 0 -f_cols 1 \
  > ${f%txt}genes.txt

rm ${f}.tmp*
