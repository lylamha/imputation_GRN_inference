
library(reshape2)
library(ggplot2)
library(yaml)
library(igraph)
library(data.table)
library(RColorBrewer)


setwd(".")
source("meta-scripts/utils.R")

out_dir <- "my/output_dir/"

# ======================================================
# summary network statistics
# num_genes, num_tfs, network density
# ======================================================

net_stats <- lapply(config_files, function(file) {
  
  ds <- gsub(".yaml", "", basename(file))
  
  config      <- read_yaml(file)
  config_info <- parse_config(config)
  config_info$net <- "rankedEdges.csv"
  
  species <- ifelse(startsWith(ds, "h"), "human", "mouse")
  tfs <- fread(paste0("../Beeline-Networks/", species, "-tfs.csv"))
  
  config_info <- subset(config_info, algorithm != "PPCOR")
  
  net_stats <- apply(config_info, 1, function(run) {
    
    print(paste0(run["name"], " - ", run["algorithm"]))
    
    db_file   <- paste0(run["in_dirs"],  "/", run["trueEdges"])
    net_file  <- paste0(run["out_dirs"], "/", run["net"])
    
    
    db   <- fread(db_file,  data.table = F)
    net  <- fread(net_file, data.table = F)
    
    net <- get_TFEdges(net, tfs)
    net <- get_ranked_edges(net)
    
    
    net_topk <- subset(net, EdgeRank <= nrow(db))
    net_genes <- unique(c(net_topk$Gene1, net_topk$Gene2))
    
    num_genes <- length(net_genes)
    num_tfs   <- sum(net_genes %in% tfs$TF)
    
    net_dens <- nrow(net_topk) / ((num_tfs * num_genes) - num_tfs)
    
    stats_df <- data.frame(
      numGenes = num_genes,
      numTFs   = num_tfs,
      numEdges = nrow(net_topk),
      netDens  = net_dens
    )
    
    return(stats_df)
    
  })
  
  net_stats <- do.call(rbind, net_stats)
  
  config_info <- cbind(config_info, net_stats)
  config_info <- cbind(config_info,
                       colsplit(config_info$name, "_", c("data", "imputation")))
  
  net_stats <- subset(config_info, select=c("data",
                                            "imputation",
                                            "algorithm",
                                            "numGenes",
                                            "numTFs",
                                            "numEdges",
                                            "netDens"))
  
  return(net_stats)
  
})

net_stats_df <- do.call(rbind, net_stats)




# ==================================================
# plot distribution of edge weights
# ==================================================

plot_list <- lapply(net_files, function(file) {
  
  imp_name <- gsub(".*_(.*?)/.*", "\\1", file)
  alg_name <- basename(dirname(file))
  
  net <- fread(file)
  
  # get statistics
 
  p <- qplot(x = net$EdgeWeight, data=net, binwidth=.1) + 
    scale_x_log10() +
    theme_classic() +
    xlab("edge weights") +
    ggtitle(paste0(imp_name, " - ", alg_name)) +
    theme(plot.title = element_text(size = 4))
  
  return(p)
})

library(gridExtra)
png(paste0(out_dir, "edge_distributions.png"), width = 10, height = 6, res = 150, units = "in")
do.call("grid.arrange", c(plot_list, ncol=4))
dev.off()



# ==================================================
# plot BEELINE evaluation results
# ==================================================

config_files <- c("config-files/hESC.yaml",
                  "config-files/hHep.yaml",
                  "config-files/mESC.yaml",
                  "config-files/mDC.yaml",
                  "config-files/mHSC-E.yaml",
                  "config-files/mHSC-GM.yaml",
                  "config-files/mHSC-L.yaml")

res_files <- c("outputs/hESC/hESC-EPr.csv",
               "outputs/hHep/hHep-EPr.csv",
               "outputs/mESC/mESC-EPr.csv",
               "outputs/mDC/mDC-EPr.csv",
               "outputs/mHSC-E/mHSC-E-EPr.csv",
               "outputs/mHSC-GM/mHSC-GM-EPr.csv",
               "outputs/mHSC-L/mHSC-L-EPr.csv")

epr <- lapply(seq_along(config_files), function(i) {
  
  config_file <- config_files[i]
  res_file    <- res_files[i]
  
  print(i)
  
  # read config file
  config      <- read_yaml(config_file)
  config_info <- parse_config(config)
  
  # recreate df without factors
  config_info <- data.frame(lapply(config_info, as.character),
                            stringsAsFactors=FALSE)
  
  ds <- basename(dirname(res_file)) # dataset
  db <- gsub(".*_(.*)\\.csv", "\\1", basename(res_file))
  
  epr_res <- read.csv(res_file, row.names = 1, check.names = F)
  
  
  rnd_pred <- get_tfNet_density(config_info)
  rnd_pred <- rnd_pred[rownames(epr_res), colnames(epr_res)]
  
  # normalize area under epr by random predictor
  epr_res <- epr_res / rnd_pred
  
  df <- melt(as.matrix(epr_res))
  df$ds <- ds
  
  return(df)
  
})
  
epr <- do.call(rbind, epr)
epr$Var2 <- gsub(".*_(.*)", "\\1", epr$Var2)

## subsetting noimputation
noimp <- subset(epr, Var2 == "noimputation")

epr <- subset(epr, Var2 != "noimputation")

## prepare for plotting

tmp <- merge(noimp, epr, by = c("ds", "Var1"))
tmp$log2ratio <- log2(tmp$value.y / tmp$value.x)

ggplot(data = epr, aes(x=Var2, y=value, col=ds, shape=ds)) +
  geom_point(size=2) +
  geom_hline(data = noimp, aes(yintercept = value, col=ds),linetype = "dashed")  +
  # geom_hline(yintercept = 1, color="darkred", lwd=1.5)  +
  facet_wrap(~Var1, ncol=4) +
  theme_bw(base_size = 14) +
  labs(y = "EPR", 
       x = "imputation methods") +
  scale_color_manual(name="dataset", values=c("#46b4d3", "#34393c", "#db494b", "#eba432", "#86aa6d", "#ac3d1f", "#207288")) +
  scale_shape_manual(name="dataset", values=c(0, 1, 2, 3, 4, 8, 5)) +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  theme(legend.title = element_text(size = 10), 
        legend.text = element_text(size = 10))
ggsave(paste0(out_dir, "eval_epr_ratio_plot_allds.png"), width = 6, height = 4)


ggplot(data = tmp, aes(x=Var2.y, y=log2ratio, col=ds, shape=ds)) +
  geom_point(size=2) +
  geom_hline(yintercept = 0, color="darkgrey", linetype = "dashed")  +
  facet_wrap(~Var1, ncol=4) +
  theme_bw(base_size = 14) +
  labs(y = "log2ratio (epr_imputed / epr_noimp)", 
       x = "imputation methods") +
  scale_color_manual(name="dataset", values=c("#46b4d3", "#34393c", "#db494b", "#eba432", "#86aa6d", "#ac3d1f", "#207288")) +
  scale_shape_manual(name="dataset", values=c(0, 1, 2, 3, 4, 8, 5)) +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  theme(legend.title = element_text(size = 10), 
        legend.text = element_text(size = 10))
ggsave(paste0(out_dir, "eval_log2ratio_epr_ratio_plot_allds.png"), width = 6, height = 4)


#### get variances of auepr

vars_df <- lapply(unique(tmp$ds), function(ds) {
 
  # get variance if imputation (Var2.y) is fixed across algorithm
  res1 <- aggregate(log2ratio ~ Var2.y, tmp[tmp$ds == ds,], var)
  colnames(res1) <- c("variable", "variance")
  res1$class <- "imputation"
  
  # get variance if algorithm (Var1) is fixed across imputation
  res2 <- aggregate(log2ratio ~ Var1, tmp[tmp$ds == ds,], var)
  colnames(res2) <- c("variable", "variance")
  res2$class <- "algorithm"
  
  res <- rbind(res1, res2)
  res$ds <- ds
  
  return(res)
})

vars_df <- do.call(rbind, vars_df)

pVal <- wilcox.test(variance ~ class, vars_df, alternative="greater")$p.value
pVal <- signif(pVal,3)
ggplot(vars_df, (aes(x=class, y=variance))) +
  geom_violin(col="grey", fill="lightgrey", alpha=.4) +
  geom_jitter(aes(shape=ds, col=ds),width = 0.1) +
  scale_color_manual(name="dataset", values=c("#46b4d3", "#34393c", "#db494b", "#eba432", "#86aa6d", "#ac3d1f", "#207288")) +
  scale_shape_manual(name="dataset", values=c(0, 1, 2, 3, 4, 8, 5)) +
  theme_classic(base_size = 14) +
  xlab("fixed method") +
  theme(legend.position = c(0.7, 0.7)) +
  ggtitle(label="AUEPR log-fold-ratio variances", subtitle = paste0("wilcoxon p-val = ", pVal))
ggsave(paste0(out_dir, "eval_log2ratio_auepr_ratio_variances.png"), width = 4, height = 5)




stats_epr_res <- lapply(seq_along(config_files), function(i) {
  
  config_file <- config_files[i]
  res_file    <- res_files[i]
  
  print(i)
  
  # read config file
  config      <- read_yaml(config_file)
  config_info <- parse_config(config)
  config_info$trueEdges <- "STRING-network.csv"
  
  
  ds <- basename(dirname(res_file)) # dataset
  db <- gsub(".*_(.*)\\.csv", "\\1", basename(res_file))
  
  epr_res <- read.csv(res_file, row.names = 1, check.names = F)
  
  net_stats <- get_tfNet_stats(config_info)
  net_stats <- unique(net_stats)
  
  
  ### collect all results in a data.frame
  res <- melt(as.matrix(epr_res))
  colnames(res) <- c("algorithm", "name", "epr")
  res <- merge(res, net_stats, by = "name")
  
  # normalize area under epr by random predictor
  res$epr_ratio <- res$epr / res$netDens
  
  return(res)
  
})


stats_epr_res <- do.call(rbind, stats_epr_res)

