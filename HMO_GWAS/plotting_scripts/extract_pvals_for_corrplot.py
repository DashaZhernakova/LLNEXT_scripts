import pandas as pd
import numpy as np

loci_coord = {"CADM1": [11,115039938, 115375675], "GRIN2B" : [12, 13693165, 14133053], "FUT2" : [19, 49199228, 49209207], "FUT3_FUT6" : [19, 5830621, 5851485], "ABO" : [9, 136125788, 136150617], "ST3GAL6" : [3, 98451080, 98540045]} 3:98608445

hmos = ["mother_milk_HMO_2FL_ugml_invr", "mother_milk_HMO_3FL_ugml_invr", "mother_milk_HMO_LDFT_ugml_invr", "mother_milk_HMO_3GL_ugml_invr", "mother_milk_HMO_6GL_ugml_invr", "mother_milk_HMO_A_tetra_ugml_invr", "mother_milk_HMO_3SL_ugml_invr", "mother_milk_HMO_LNT_ugml_invr", "mother_milk_HMO_LNnT_ugml_invr", "mother_milk_HMO_6SL_ugml_invr", "mother_milk_HMO_3F3SL_ugml_invr", "mother_milk_HMO_LNFP_V_ugml_invr", "mother_milk_HMO_LNnFP_V_ugml_invr", "mother_milk_HMO_LNFP_I_ugml_invr", "mother_milk_HMO_LNFP_III_ugml_invr", "mother_milk_HMO_LNFP_II_ugml_invr", "mother_milk_HMO_LNnDFH_ugml_invr", "mother_milk_HMO_LSTb_ugml_invr", "mother_milk_HMO_LNDH_I_ugml_invr", "mother_milk_HMO_LSTc_ugml_invr", "mother_milk_HMO_LNH_ugml_invr", "mother_milk_HMO_MFLNH_III_ugml_invr", "mother_milk_HMO_DSLNT_ugml_invr", "mother_milk_HMO_DFLNHa_ugml_invr", "mother_milk_HMO_Total_ugml_invr", "mother_milk_HMO_Fuc_ugml_invr", "mother_milk_HMO_Neut_ugml_invr", "mother_milk_HMO_Sia_ugml_invr"]

timepoints = ["W2", "M1", "M2", "M3", "M6"]

for tp in timepoints:
    for hmo in hmos:
        fname = "/groups/umcg-llnext/tmp01/umcg-dzhernakova/HMO_GWAS/results_plink_imp/" + tp + "/HMO_GWAS_res_" + tp + "." + hmo + ".glm.linear.gz"
        f = pd.read_csv(fname,  sep = "\t", header = 0)
        
        for (locus, coord) in loci_coord.items():
            res_per_locus = f.loc[(f['#CHROM'] == coord[0]) & (f['POS'] < (coord[2] + 250000)) & (f['POS'] > (coord[1] - 250000))]
            if len(res_per_locus.index) > 0:
                print(tp, hmo, locus, min(res_per_locus['P']))

