module load BCFtools
module load PLINK

#
# 1. Convert vcf to plink
#
for chr in `seq 1 22`
do
 plink2 \
 --vcf /groups/umcg-llnext/tmp01/data/imputed/rawdata/LLnextv2_2920merged23122022.vcfs/${chr}.pbwt_reference_impute.vcf.gz \
 --extract-if-info "INFO > 0.5" \
 --make-bed --out ${chr}.info_filtered
 
 awk 'BEGIN {FS=OFS="\t"}; {if ($2 == ".") $2 = $1 ":" $4; print $0}' ${chr}.info_filtered.bim > ${chr}.info_filtered.bim.tmp
 mv ${chr}.info_filtered.bim.tmp ${chr}.info_filtered.bim
 
 awk 'BEGIN {FS=OFS="\t"}; {$2 = $1 ":" $4; print $0}' ${chr}.info_filtered.bim > ${chr}.info_filtered.bim.tmp
 mv ${chr}.info_filtered.bim.tmp ${chr}.info_filtered.bim
 
done


#
# 2. Merge per chromosome plink files
#

rm merge_list.txt
for chr in `seq 1 22`
do
 echo "${chr}.info_filtered" >> merge_list.txt
done

module load PLINK/1.9-beta6-20190617
plink --merge-list merge_list.txt --make-bed --out all_chr.info_filtered

cat all_chr.filtered-merge.missnp > multiallelic_SNPs.txt

for chr in `seq 1 22`
do
 plink --bfile ${chr}.info_filtered \
 --exclude multiallelic_SNPs.txt \
 --make-bed --out ${chr}.info_filtered.nodup
done

rm merge_list.txt
for chr in `seq 1 22`
do
 echo "${chr}.info_filtered.nodup" >> merge_list.txt
done

plink --merge-list merge_list.txt --make-bed --out all_chr.info_filtered


#
# 3. Extract mothers and filter SNPs
#

plink2 \
    --bfile all_chr.info_filtered \
    --keep qc/mothers.txt \
    --maf 0.05 --hwe 1e-6 --geno 0.05 \
    --make-bed \
    --out all_chr.mothers.with_rel.flt
    
