d=/groups/umcg-llnext/tmp01/umcg-dzhernakova/HMO_GWAS/
cd ${d}
script_dir=/groups/umcg-llnext/tmp01/umcg-dzhernakova/LLNEXT_scripts/HMO_GWAS/scripts/

timepoints=(W2 M1 M2 M3 M6)

#
# 1. Manhattan plot
#
for tp in ${timepoints[@]}
do
 cp hmo_imp_plink_for_plots.bottom_p0.05.txt ${tp}.non_imputed.for_plotting.txt
 for f in ../results_plink/${tp}/HMO_GWAS_res_${tp}.*.glm.linear.gz
 do
  zcat $f | awk 'BEGIN {FS=OFS="\t"}; {if ($12 < 0.05) print}' >> ${tp}.non_imputed.for_plotting.txt
 done
done

ml RPlus/4.0.3-foss-2018c-v21.12.10
for tp in ${timepoints[@]}
do
  Rscript ${script_dir}/plotting_scripts/plot_manhattan.R ${tp}.non_imputed.for_plotting.txt
done

#
# 2. Corplot
#
python3 ${script_dir}/plotting_scripts/extract_pvals_for_corrplot.py
Rscript ${script_dir}/plotting_scripts/corplot_HMO.R 