# Effect of imputation on network reconstruction in single cell RNAseq data
Project to study the effect of imputation methods on GRN reconstruction in single cell transcriptome data

### Data
We collected processed and normalised single-cell RNA-seq data from https://zenodo.org/record/3701939#.X6UYsHWYW-Y by Pratapa et. al., Nature Methods, 2020. 
Furtheremore, we downloaded the corresponding inferred pseudotime for each dataset/ cell type.

### Imputation
We use dca, knn-smoothing, MAGIC and SAVER to impute data on experimental scRNAseq datasets


### GRN reconstruction
We use PIDC, GENIE3 and GRNBoost2 to reconstruct GRNs from imputed and unimputed data.
Here, we make use of the evaluation framework BEELINE by Pratapa et. al. in order to reconstruct the networks and to evaluate their performance using the STRING database - a functional protein protein interaction network.
