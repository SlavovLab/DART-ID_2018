DART-ID (2018)
Chen, Albert & Franks, Alexander & Slavov, Nikolai
--------------------------------------------------

Unless otherwise specified, data used for this paper is from the updated evidence.txt file (ev_updated.txt) in ftp://massive.ucsd.edu/MSV000083149/other/Alignments/SQC_20180815_2/. The MaxQuant output files are specified in the config file, SQC.yaml, and are available in ftp://massive.ucsd.edu/MSV000083149/other/MaxQuant/. Raw files are also available in .raw and mzXML form as part of the MassIVE submission.

---------
# Figures
---------

## Figure 1:

### Data

All data used for these diagrams are simulated. Values such as null distribution height and experiment shifts are exaggerated for visual clarity.

### Scripts: (GitHub/Rscripts)

* 1A: fig_scripts/fig_1a_canonical_demo.R
* 1B: fig_scripts/fig_1b_alignment_demo_v5.R
* 1C: fig_scripts/fig_1c_confidence_update_v2.R




## Figure 2:

### Data:

* MaxQuant search results: SQC_67_95_Varied -- ftp://massive.ucsd.edu/MSV000083149/other/MaxQuant/SQC_67_95_Varied/
* DART-ID Alignment, and results from other RT inference methods: SQC_varied_20180711_4/ -- ftp://massive.ucsd.edu/MSV000083149/other/Alignments/SQC_varied_20180711_4/
* iRT search and converted data: dat/irt.txt, dat/irt.xls (ftp://massive.ucsd.edu/MSV000083149/other/Rdata/)
* R objects: error_df.RData (ftp://massive.ucsd.edu/MSV000083149/other/Rdata/)

### Scripts: (GitHub/Rscripts)

* alignment_comparisons.R -- load and process alignment comparison data
* 2A-B: fig_scripts/fig_2_alignment_comparison.R




## Figure 3:

### Data:

Figure 3B. Contains only the first 50 experiments (sorted chronologically, from run date) of the entire dataset, for visual clarity and space constraints. Figures 3A, C, D, and E use data from the entire set.

The proportions and magnitudes of numbers from Figure 3B. do shift when including more experiments, but not significantly.

### Scripts: (GitHub/Rscripts)

* add_percolator.R -- add results from percolator output
* protein_quant.R -- Split data into PSM sets, get protein quantification from PSM quantification
* 3A: fig_scripts/fig_3a_pep_scatter_v3.R
* 3B: fig_scripts/fig_3b_protein_map_v3_4.R
* 3C: fig_scripts/fig_3c_fdr_fold_change
* 3D: fig_scripts/fig_3d_protein_ids_v2.R
* 3E: fig_scripts/fig_3e_pep_dists.R




## Figure 4:

### Data:

* MaxQuant search results: FP18 -- ftp://massive.ucsd.edu/MSV000083149/other/MaxQuant/FP18/
* DART-ID alignment results: All sets: ftp://massive.ucsd.edu/MSV000083149/other/Alignments/FP_validation_allsets_20180706_1/, Set 1: ftp://massive.ucsd.edu/MSV000083149/other/Alignments/FP_validation_set1_20180706_2/, Set 2: ftp://massive.ucsd.edu/MSV000083149/other/Alignments/FP_validation_set2_20180706_1/

### Scripts: (GitHub/Rscripts)

* 4B-C: fig_scripts/fig_4_validate_new_peptides_v2.R




## Figure 5:

### Scripts: (GitHub/Rscripts)

* add_percolator.R -- add results from percolator output
* 5B-C: fig_scripts/fig_5_cv_valiation.R




## Figure 6:

### Data:

Manually removed experiments that were missing a significant amount of TMT quantitation, in order to lessen the amount of missing data in the final expression matrix.

Experiments used for the PCA plot: 

180320S_QC_SQC67C8, 180324S_QC_SQC67C10, 180324S_QC_SQC67C11, 180324S_QC_SQC67C12, 180324S_QC_SQC67C13, 180324S_QC_SQC67C14, 180324S_QC_SQC67C15, 180324S_QC_SQC67C17, 180324S_QC_SQC67C18, 180324S_QC_SQC67C19, 180324S_QC_SQC67C9, 180324S_QC_SQC68B1, 180324S_QC_SQC68C1, 180324S_QC_SQC68C2, 180324S_QC_SQC68C3, 180324S_QC_SQC68D2, 180324S_QC_SQC68D3, 180324S_QC_SQC68D4, 180324S_QC_SQC68D4_20180330185431, 180324S_QC_SQC68E2, 180324S_QC_SQC68E3, 180324S_QC_SQC68E4, 180324S_QC_SQC68E5, 180324S_QC_SQC68E6, 180324S_QC_SQC68F1, 180324S_QC_SQC68F2, 180324S_QC_SQC68F3, 180324S_QC_SQC69A, 180324S_QC_SQC69B, 180324S_QC_SQC69C, 180324S_QC_SQC69D, 180324S_QC_SQC69E, 180324S_QC_SQC69F, 180324S_QC_SQC70D_20180330082946, 180324S_QC_SQC70E, 180324S_QC_SQC71D, 180324S_QC_SQC71E2, 180402S_QC_SQC72A1, 180402S_QC_SQC72B1, 180402S_QC_SQC72B10, 180402S_QC_SQC72B11, 180402S_QC_SQC72B12, 180402S_QC_SQC72B2, 180402S_QC_SQC72B3, 180402S_QC_SQC72B7, 180402S_QC_SQC72C1, 180402S_QC_SQC72C2, 180402S_QC_SQC72C3, 180402S_QC_SQC72C4, 180402S_QC_SQC72C5, 180402S_QC_SQC72C6, 180402S_QC_SQC72C7, 180402S_QC_SQC72C8, 180402S_QC_SQC72D1, 180402S_QC_SQC72D2, 180402S_QC_SQC72D3, 180402S_QC_SQC72D4, 180402S_QC_SQC73C1, 180402S_QC_SQC73D1, 180406S_QC_SQC74A, 180406S_QC_SQC74B, 180406S_QC_SQC74B2, 180406S_QC_SQC74C, 180406S_QC_SQC74E, 180406S_QC_SQC74F, 180406S_QC_SQC74G, 180406S_QC_SQC74H, 180406S_QC_SQC74J, 180409S_QC_SQC76D, 180409S_QC_SQC76E, 180420S_QC_SQC79A, 180420S_QC_SQC79B, 180420S_QC_SQC79C, 180420S_QC_SQC79D, 180420S_QC_SQC79F, 180420S_QC_SQC80B, 180420S_QC_SQC80C, 180420S_QC_SQC80D, 180420S_QC_SQC80E, 180420S_QC_SQC80F, 180420S_QC_SQC80G, 180424S_QC_SQC81A, 180424S_QC_SQC81C, 180502S_QC_SQC82B, 180502S_QC_SQC82C, 180502S_QC_SQC82D, 180502S_QC_SQC82E, 180502S_QC_SQC82F, 180502S_QC_SQC82G, 180502S_QC_SQC82H, 180502S_QC_SQC82I, 180503S_QC_SQC83A, 180503S_QC_SQC83B, 180503S_QC_SQC84A, 180503S_QC_SQC84A2, 180503S_QC_SQC84B, 180503S_QC_SQC85A1, 180503S_QC_SQC85A2, 180503S_QC_SQC85A3, 180503S_QC_SQC85B1, 180503S_QC_SQC85B2, 180503S_QC_SQC85B3, 180503S_QC_SQC86A, 180503S_QC_SQC86B, 180503S_QC_SQC86C, 180508S_QC_SQC87A1, 180508S_QC_SQC87A10, 180508S_QC_SQC87A2, 180508S_QC_SQC87A3, 180508S_QC_SQC87A4, 180508S_QC_SQC87A5, 180508S_QC_SQC87A6, 180508S_QC_SQC87A7, 180508S_QC_SQC87A8, 180508S_QC_SQC87A9, 180508S_QC_SQC87C1, 180508S_QC_SQC87C10, 180508S_QC_SQC87C2, 180508S_QC_SQC87C3, 180508S_QC_SQC87C4, 180508S_QC_SQC87C5, 180508S_QC_SQC87C6, 180508S_QC_SQC87C7, 180508S_QC_SQC87C8, 180508S_QC_SQC87C9

The code to filter the data is found in the figure generation script(s) below.

### Scripts: (GitHub/Rscripts)

* validation_cormats.R -- helper scripts for filtering/wrangling data
* 6A-B: fig_scripts/fig_6_cor_pca.R








---------------------
# Supporting Figures:
---------------------

## Supporting Figure 1:

### Data:

All data was simulated

### Scripts: (GitHub/Rscripts)

* fig_scripts/sfig_1_mixture_model_v2.R




## Supporting Figure 2:

### Data:

All data was simulated

### Scripts: (GitHub/Rscripts)

* fig_scripts/sfig_2_bootstrap.R




## Supporting Figure 3:

* fig_scripts/sfig_3_compare_linear.R




## Supporting Figure 4:

* fig_scripts/sfig_4_exp_rt_variance.R




## Supporting Figure 5:

* fig_scripts/sfig_5_distribution_choice.R




## Supporting Figure 6:

* fig_scripts/fig_ratio_scatter.R



