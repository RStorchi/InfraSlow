# InfraSlow

Infra Classification: 
Use batch_infra_class_preprocess.m to pre-process the data. Then use infra_classification_review.m to classify the data. 
The classification results will be saved in folder named “classification”. 


Fast Beta/Gamma Classification: 
Use batch_gamma_class_preprocess.m to pre-process the data. Then use gamma_classification_review.m to classify the data. 
The classification results will be saved in folder named “classification” . 


Phase Coupling Analysis: 
Use infra_phase_gamma_pow_review.m. Results are saved in folder named “Data” as RC_PAC_results.mat (for wild type animals) 
and as MELKO_PAC_results.mat for the melanopsin knockout genotype. 
To calculate cross-correlations use infra_crosscorr_review.m. Results are provided on the MATLAB shell. 


Permutation Test: 
Use permutation_test_review.m. 
To compare firing rates across neuron types (none, infra only, fast only and both infra and fast) use firing_test_review.m.


Fano Factor Analysis:
To determine Fano Factor as function of time windows use FF_review.m. Data are saved in folder "Data". 
To plot results use FF_figure_review.m.
