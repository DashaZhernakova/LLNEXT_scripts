timepoints=("M1" "M2" "M3" "M6" "W2")
for t in ${timepoints[@]}
do
 awk -v tp=$t 'BEGIN {FS=OFS="\t"}; {if (($4 == tp && $3 == "first_pregnancy") || NR == 1) print $0}' 230425_HMO_data_milk_groups_ugml_cleaned_n1541_Imputation.txt | cut -f2,9- > HMO_${t}.txt 
 Rscript adjust_for_PCs.noINT.R  HMO_${t}.txt
done