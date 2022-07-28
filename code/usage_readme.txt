##See data folder for details of input file
Rscript tum_normal_matched.R input_paired
Rscript tumour_single_sample.R input_single 
Edit paths in HRD_classification.R to point to appropriate inputs (generated in the previous steps) and output (prediction results)
