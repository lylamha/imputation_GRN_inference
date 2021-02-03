#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("The directoriy folder must be supplied.n", call.=FALSE)
}


library(data.table)


out_dir <- args[1]

net_files <- list.files(path = out_dir, pattern = "rankedEdges.csv", full.names = T, recursive = T)
net_files <- grep("tmp/rankedEdges.csv", net_files, invert = T, value = T)

print(paste0(length(net_files), " found in ", out_dir, "."))


for (file in net_files) {
  
  alg_name <- basename(dirname(file))
  dat_name <- basename(out_dir)
  # print(paste0(dat_name, "-", alg_name))
  
  fcp <- paste0("./", dirname(file), "/tmp/rankedEdges.csv") 
  
  ### copy existing rankedEdges.csv file into tmp
  if( file.exists(fcp) ) {
    print(paste0(basename(file), " exists in /tmp already. Skip copying."))
  } else {
    tmp_dir <- paste0("./", dirname(file), "/tmp/")
    if (!dir.exists(tmp_dir)) { dir.create(tmp_dir) }
    fcp <- file.copy(file, dirname(fcp))
  }
  
  
  ### read output network file created by beeline
  net <- fread(file, data.table = F)
  
  
  ### get TF info
  species <- ifelse(startsWith(dat_name, "h"), "human", "mouse")
  
  tf_file <- paste0("../Beeline-Networks/", species, "-tfs.csv")
  
  if(!file.exists(tf_file)) stop(tf_file, " does not exist.")
  
  tfs <- fread(tf_file)
  
  ### consider only outgoing edges from TF
  net$gene1_tf <- ifelse(net$Gene1 %in% tfs$TF, T, F)
  
  net <- subset(net, gene1_tf, select = c("Gene1", "Gene2", "EdgeWeight"))
  # print(nrow(net))

  # overwrite
  print(paste("overwriting", file))
  write.table(net, file=file, quote=F, row.names = F, sep = "\t")
  
}
