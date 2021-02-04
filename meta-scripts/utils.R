


# ==========================================================
# parse config file
# ==========================================================


parse_config <- function(config) {
  
  in_dir  <- paste0(config$input_settings$input_dir, "/",
                    config$input_settings$dataset_dir)
  
  out_dir <- paste0(config$output_settings$output_dir, "/",
                    config$output_settings$output_prefix)
  
  datasets <- lapply(config$input_settings$datasets, function(ds_i) {
    
    return(cbind(name=ds_i$name,
                 exprData=ds_i$exprData,
                 trueEdges=ds_i$trueEdges))
    
  })
  
  df <- do.call(rbind, datasets)
  df <- as.data.frame(df)
  
  algos <- lapply(config$input_settings$algorithms, function(alg_i) {
    
    if (alg_i$params$should_run) {
      return(alg_i$name)
    }
    
  })
  
  algos <- unlist(algos)
  
  # crossing combinations of datasets and algorithms
  dsXalgos <- expand.grid(paste0(out_dir, "/", df$name), algos)
  dsXalgos <- dsXalgos[order(dsXalgos$Var1),]
  
  # create data.frame only with directories
  out_dirs <- paste0(dsXalgos$Var1, "/", dsXalgos$Var2)
  # in_dirs  <- paste0(in_dir, "/", df$name)
  in_dirs  <- gsub("outputs", "inputs", dsXalgos$Var1)
  
  df_dirs <- data.frame(name = gsub(paste0(out_dir, "/"), "", dirname(out_dirs)),
                        algorithm = basename(out_dirs), 
                        out_dirs  = out_dirs,
                        in_dirs   = in_dirs)
  
  # merge dfs (input files and corresponding dirs)
  df <- merge(df, df_dirs, by="name")
  
  return(df)
  
}



# ==========================================================
# Compute network density from config file considering
# TFs in a network
# ==========================================================


get_tfNet_density <- function(config_info) {
  
  # read in TF file
  species <- ifelse(startsWith(config_info$name[1], "h"), "human", "mouse")
  tfs <- fread(paste0("../Beeline-Networks/", species, "-tfs.csv"))
  
  net_dens <- apply(config_info, 1, function(run) {
    
    # print(paste0(run["name"], " - ", run["algorithm"]))
    
    expr_file <- paste0(run["in_dirs"], "/", run["exprData"])
    db_file   <- paste0(run["in_dirs"], "/", run["trueEdges"])
    
    expr <- read.csv(expr_file, row.names = 1)
    db <- fread(db_file, data.table = F)
    
    
    db_genes <- unique(c(db$Gene1, db$Gene2))
    db_genes <- db_genes[db_genes %in% rownames(expr)]
    
    db <- subset(db, Gene1 %in% db_genes & Gene2 %in% db_genes)
    
    numTFs <- sum(tfs$TF %in% db_genes)
    numGenes <- length(db_genes)
    
    # To compute network density, we computed the number of edges divided
    # by the total number of edges that can be outgoing from a TF
    net_dens <- nrow(db) / ((numTFs * numGenes) - numTFs)
    
    return(net_dens)
    
  })
  
  net_dens <- matrix(net_dens,
                     nrow = length(unique(config_info$name)),
                     ncol = length(unique(config_info$algorithm)),
                     byrow = T)
  net_dens <- as.data.frame(net_dens)
  rownames(net_dens) <- unique(config_info$name)
  colnames(net_dens) <- unique(config_info$algorithm)
  
  net_dens <- t(net_dens)
  
  return(net_dens)
  
}

get_tfNet_stats <- function(config_info) {
  
  # read in TF file
  species <- ifelse(startsWith(config_info$name[1], "h"), "human", "mouse")
  tfs <- fread(paste0("../Beeline-Networks/", species, "-tfs.csv"))
  
  stats_df <- apply(config_info, 1, function(run) {
    
    # print(paste0(run["name"], " - ", run["algorithm"]))
    
    expr_file <- paste0(run["in_dirs"], "/", run["exprData"])
    db_file   <- paste0(run["in_dirs"], "/", run["trueEdges"])
    
    expr <- read.csv(expr_file, row.names = 1)
    db <- fread(db_file, data.table = F)
    
    
    db_genes <- unique(c(db$Gene1, db$Gene2))
    db_genes <- db_genes[db_genes %in% rownames(expr)]
    
    db <- subset(db, Gene1 %in% db_genes & Gene2 %in% db_genes)
    
    numTFs   <- sum(tfs$TF %in% db_genes)
    numGenes <- length(db_genes)
    
    # To compute network density, we computed the number of edges divided
    # by the total number of edges that can be outgoing from a TF
    net_dens <- nrow(db) / ((numTFs * numGenes) - numTFs)
    
    stats_df <- data.frame(name     = run["name"],
                           numTFs   = numTFs,
                           numGenes = numGenes,
                           numEdges = nrow(db),
                           netDens  = net_dens)
    
    return(stats_df)
    
  })
  
  stats_df <- do.call(rbind, stats_df)
  rownames(stats_df) <- NULL
  
  return(stats_df)
  
}



# =====================================================================


# rank edges, assign the same rank if there is a tie
get_ranked_edges <- function(net) {
  
  net$EdgeRank <- rank(-net$EdgeWeight)
  net$EdgeRank <- as.numeric(factor(net$EdgeRank))
  
  return(net)
}


# identify the true Edges in predicted network
find_trueEdges <- function(net, db) {
  
  net_edges <- paste(net$Gene1, net$Gene2, sep = "|")
  db_edges  <- paste(db$Gene1, db$Gene2, sep = "|")
  
  # assign 1 if edge is found in db network
  net$TrueEdge <- 0
  net$TrueEdge[net_edges %in% db_edges] <- 1
  
  return(net)
}


# filter only for outgoing TF edges
get_TFEdges <- function(net, tfs) {
  
  return(net[net$Gene1 %in% tfs$TF,])
  
}
