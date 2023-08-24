
timepoints=(W2 M1 M2 M3 M6 M12)
data_types=(clr quant binary)


cd /groups/umcg-llnext/tmp01/umcg-dzhernakova/mGWAS/results
for dt in ${data_types[@]}
do
    res_f=${dt}.LCT_ABO_FUT2.5e-3.txt
    rm $res_f
    echo $dt
    while read line 
    do
        sp=$line
        echo $sp
        for tp in ${timepoints[@]}
        do
            f=${tp}/${dt}/babies_${tp}_${dt}_${sp}.fastGWA.gz
            zcat $f | awk -v tp=$tp -v sp=$sp 'BEGIN {FS=OFS="\t"}; {if ($1 == 2 && $3 > 136545410 && $3 < 136594750 && $10 < 0.0005) {print "LCT", tp, sp, $0} else if ($1 == 9 && $3 > 136125788 && $3 < 136150617 && $10 < 0.0005) {print "ABO", tp, sp, $0} else if ($1 == 19 && $3 > 49199228 && $3 < 49209207 && $10 < 0.0005) {print "FUT2", tp, sp, $0}}' >> $res_f
        done
    done < /groups/umcg-llnext/tmp01/umcg-dzhernakova/mGWAS/data/all_species.txt
done