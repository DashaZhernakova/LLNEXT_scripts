ml PLINK

plink2 \
    --bfile all_chr.info_filtered \
    --keep qc/mothers.txt \
    --maf 0.05 --hwe 1e-6 --geno 0.05 \
    --make-bed \
    --out all_chr.mothers.with_rel.flt
    

#
# Make GRMs for fastGWA
#
cd /groups/umcg-llnext/tmp01/umcg-dzhernakova/HMO_GWAS/genotypes
ln -s /groups/umcg-llnext/tmp01/umcg-dzhernakova/genotypes/all_chr.mothers.with_rel.flt.bed
ln -s /groups/umcg-llnext/tmp01/umcg-dzhernakova/genotypes/all_chr.mothers.with_rel.flt.bim
ln -s /groups/umcg-llnext/tmp01/umcg-dzhernakova/genotypes/all_chr.mothers.with_rel.flt.fam



# make GRMs for GCTA

gcta=/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/v2//tools/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1
grm=GRM/mothers_flt_pruned
mkdir GRM

plink2 --bfile all_chr.mothers.with_rel.flt --indep-pairwise 250 50 0.2 --out tmp_pruning
$gcta --bfile all_chr.mothers.with_rel.flt  --extract tmp_pruning.prune.in  --make-grm --out $grm
$gcta --grm $grm --make-bK-sparse 0.05 --out ${grm}_sparse 
$gcta --bfile all_chr.mothers.with_rel.flt --extract tmp_pruning.prune.in --make-grm-gz --out ${grm}.text

rm tmp.pruning*
