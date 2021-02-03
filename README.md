# Effect of imputation on network reconstruction in single cell RNAseq data
Project to study the effect of imputation methods on GRN reconstruction in single cell transcriptome data

### Data

We collected processed and normalised single-cell RNA-seq data from [here](https://zenodo.org/record/3701939#.X6UYsHWYW-Y) (Pratapa et. al., Nature Methods, 2020). 
Here, each data folder consists of three files: ExpressionData.csv (normalised counts), GeneOrdering.csv (genes ordered by their variance and p-value) and Pseudotime.csv (inferred pseudotime for each cell).

### Imputation
We use dca (version 0.2.3), knn-smoothing (version 2.1), MAGIC ('Rmagic' R package version 2.0.3) and SAVER (R package version 1.1.2) to impute data on experimental scRNAseq datasets with the following example commands:

dca:
`dca </path/to/ExpressionData_raw.csv> </path/to/dca_result_folder>` (Please note, that dca only runs on raw count data.)

knn-smoothing:
```
python3 knn_smooth.py -k 15 -d 10 \
-f <path/to/ExpressionData.csv> \
-o <path/to/ExpressionData_knnsmooth_imputed.csv> --sep ,
```

For magic and saver see Rscripts.



### Gene filtering

We used `generateExpInputs.py` to filter the expression matrices according to the top500 most variable genes and significantly varying TFs (Bonferroni corrected p < 0.01).  

```
python generateExpInputs.py \
-e=ExpressionData.csv \
-g=GeneOrdering.csv \
-f=STRING-network.csv \
-i=mouse-tfs.csv \
-p=0.01 \
-c \
-n=500 \
-t \
-o=500HVG_TFs
```


### GRN reconstruction
We use PIDC, GENIE3 and GRNBoost2 to reconstruct GRNs from imputed and unimputed data.
Here, we make use of the evaluation framework BEELINE by Pratapa et. al. in order to reconstruct the networks and to evaluate their performance using the STRING database - a functional protein protein interaction network.
