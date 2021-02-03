
library(Rmagic)
library(Seurat)
library(reticulate)

# data files
data_dir <- "."
expr_files <- list.files(".", pattern = "ExpressionData.csv", recursive = TRUE, include.dirs = TRUE)

impList <- lapply(expr_files, function(file) {
  
  
  dat <- read.csv(paste0(data_dir, "/", file), row.names = 1)
  dat_name <- dirname(file)
  
  print(dat_name)
  
  # remove erccs
  ercc.idx <- grep("^ERCC-", rownames(dat))
  if( length(ercc.idx) ) { dat <- dat[-ercc.idx,] }
  
  nGenes <- nrow(dat)
  
  # seurat object
  so_dat <- CreateSeuratObject(counts = dat, min.cells = 1)
 
  # normalisation and scaling
  so_dat <- NormalizeData(so_dat, assay = "RNA")
  so_dat <- ScaleData(so_dat, assay = "RNA")
  
  # variance stabilisation
  so_dat <- FindVariableFeatures(so_dat, selection.method = "vst")

  # run magic
  so_dat <- magic(so_dat, assay="RNA", genes="all_genes")
  
  
  new_file <- gsub(".csv", "_magic_imputed.csv", file)
  new_file <- paste0(data_dir, "/", new_file)
  write.csv(as.data.frame(so_dat@assays$MAGIC_RNA@data),
            file = new_file, quote=F, row.names = T)
  

  return(new_file)
  
})

