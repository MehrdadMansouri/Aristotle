# Aristotle
source code for Aristotle: Stratified Causal Discovery for Omics Data
Aristotle is an approach for Stratified Causal Discovery for omics data. Aristotle employs existing biological knowledge and a state-of-the-art patient stratification method to tackle the above challenges and applies a quasi-experimental design method to each stratum to find stratum-specific potential causes.

Files are described in the order that they should be run.

"Preprocessing.m" uses functions "Gene2SNP.m", "Pathway2SNP.m", "Patient2SNP.m", and "SNP2PathwayAll.m" to identify the pathway-based groupings via pathway-gene and gene-SNP relations, and outputs and exports files into directory "Dir" files that can be used by "Filtering.m". It was used for the Anthracycline application dataset and any alternative grouping can be used based on the application.

"Filtering.m" reads files from the directory "Dir" and groups features, and for each group, based on the weights, outputs the candidate features in table "Table", and exports it as 'Feature_Selections.txt'.

"QED.m" performs the quasi-experiment design on the candidate features and strata of "Filtering.m" and outputs their corresponding significance levels "P_Value" and q-q plot.

"FDR.m" performs Benjamini-Hochberg correction based on the input "P_values", discovers the causes that pass the multiple-hypotheses adjusted threshold and outputs the q-q plot.

"Num_Hyp_Estimation.m" estimates the number of hypotheses for a set of p-values "P" based on (Hwang 2014). It was only used in the evaluation to check which assumption about the number of hypotheses is correct.

"Synthetic_Simulation.m" generates the synthetic dataset with stratified population and causes. Parameters can be adusted at the beginning accordingly.

"Synthetic_Filtering.m" applies Aristotle to the series of datasets produced by "Synthetic_Simulation.m", and exports the output for each dataset respectively.

"Synthetic_Evaluation.m" evaluates the results of "Synthetic_Filtering.m" across different runs and compares them against the competing methods.
