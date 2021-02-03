# Effect of imputation on network reconstruction in single cell RNAseq data
Project to study the effect of imputation methods on gene regulatory network (GRN) reconstruction  in single cell transcriptome data

### Data

We collected processed and normalised single-cell RNA-seq data from [here](https://zenodo.org/record/3701939#.X6UYsHWYW-Y) (Pratapa et. al., Nature Methods, 2020). 
Here, each data folder consists of three files:
- ExpressionData.csv (normalised counts)
- GeneOrdering.csv (genes ordered by their variance and respective p-value)
- Pseudotime.csv (inferred pseudotime for each cell).

As a ground truth network we consider the STRING database obtained from the same link above as well as transcription factor information for both human and mouse.

### Imputation
We use dca (version 0.2.3), knn-smoothing (version 2.1), magic ('Rmagic' R package version 2.0.3) and saver (R package version 1.1.2) to impute data on experimental scRNAseq datasets with the following example commands:

**dca** (Please note, that dca only runs on raw count data.):
```
dca </path/to/ExpressionData_raw.csv> </path/to/dca_result_folder>
```

**knn-smoothing**:
```
python3 knn_smooth.py -k 15 -d 10 \
-f <path/to/ExpressionData.csv> \
-o <path/to/ExpressionData_knnsmooth_imputed.csv> --sep ,
```

For **magic** and **saver** see Rscripts in imputation folder.


### Gene filtering

We used `generateExpInputs.py` to filter the expression matrices according to the top500 most variable genes and significantly varying TFs (Bonferroni corrected p < 0.01). This is an example command that we applied for unimputed and imputed expression matrices. We obtain a subset of ExpressionData.csv called 500HVG_TFs-ExpressionData.csv with the mentioned filter criteria and the corresponding STRING network filtered on the same subset of genes (500HVG_TFs-network.csv).

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
Applying the command on all datasets (imputed as well as unimputed data) we can now use the GRN methods.


### GRN reconstruction
We use PIDC, GENIE3 and GRNBoost2 to reconstruct GRNs from imputed and unimputed data.
Here, we make use of the evaluation framework [BEELINE](https://github.com/Murali-group/Beeline) by Pratapa et. al. in order to reconstruct the networks:

We executed the following command as suggested in the BEELINE tutorial chainging the `config.yaml` files accordingly (see `config-files` folder).

```
python BLRunner.py --config config-files/config.yaml
```

We filtered out TF interactions only using `meta-scripts/get_TF_edges.R` by

```
Rscript --vanilla meta-scripts/get_TF_edges.R <folder/to/rankedEdges/directory>
```

Afterwards we evaluated the predicted networks with:
`python BLEvaluator.py --config config-files/config.yaml --epr`


### Evaluation

