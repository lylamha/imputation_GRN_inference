# imputation_GRN_inference
Project to study the effect of imputation methods on GRN reconstruction in single cell transcriptome data

### Imputation
We use dca, knn-smoothing, MAGIC and SAVER to impute data on experimental scRNAseq datasets


### GRN reconstruction
We use PIDC, GENIE3 and GRNBoost2 to reconstruct GRNs from imputed and unimputed data.
Here, we make use of the evaluation framework BEELINE by Pratapa et. al. in order to reconstruct the networks and to evaluate their performance using the STRING database - a functional protein protein interaction network.
