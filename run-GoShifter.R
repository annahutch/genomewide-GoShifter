# script to generate GoShifter statistics for each LD block 
#  (to be ran as a slurm array job for 1 to the number of LD blocks)

#####################

### replace these with the correct file names

p_file <- "replace_this_with_path_to_p.RDS" # vector of GWAS p-values
q_file <- "replace_this_with_path_to_q.RDS" # vector of binary annotations (for the SNPs in the same order as p_file)
blocks_file <- "replace_this_with_path_to_blocks.RDS" # vector of LD block each SNP resides (for the SNPs in the same order as p_file)
gene_dist_file <- "replace_this_with_path_to_genedist.RDS" # vector of distance to nearest gene for each SNP (for the SNPs in the same order as p_file) (should be positive values)
MAF_file <-  "replace_this_with_path_to_MAF.RDS" # vector of MAF values for each SNP (for the SNPs in the same order as p_file)

#####################

library(data.table)
source("GoShifter-funcs.R")

p_all <- readRDS(p_file)
q_all <- readRDS(q_file)
blocks <- readRDS(blocks_file)
gene_dist_all <- readRDS(gene_dist_file)
MAF_all <- readRDS(MAF_file)

#####################

# set up to run as a slurm array job for each LD block
task_id_string <- Sys.getenv("SLURM_ARRAY_TASK_ID")
i <- as.numeric(task_id_string) # 1 to the number of LD blocks

# first stratify by LD block
my_block <- unique(blocks)[i]
ind <- which(blocks == my_block)

p <- p_all[ind]
q <- q_all[ind]

# then by gene distance and MAF
gene_dist <- gene_dist_all[ind]
MAF <- MAF_all[ind]

genedist_breaks <- as.vector(quantile(gene_dist[-which(gene_dist == 0)], probs = seq(0, 1, 1/2)))
MAF_breaks <- as.vector(quantile(na.omit(MAF), probs = seq(0, 1, 1/2)))

if(all(is.na(genedist_breaks)) | any(genedist_breaks == 0) | any(duplicated(genedist_breaks))){
  
  gene_dist_bin <- cut(gene_dist, breaks = c(-Inf, Inf))
  
} else{
  
  gene_dist_bin <- cut(gene_dist, breaks = c(-Inf, genedist_breaks))
  
}

MAF_bin <- cut(MAF, breaks = c(-Inf, MAF_breaks))

# for NAs - randomly allocate a bin
MAF_bin[which(is.na(MAF_bin))] <- levels(MAF_bin)[sample(1:length(levels(MAF_bin)), length(which(is.na(MAF_bin))), replace = T)]

p_binned <- split(p, list(gene_dist_bin, MAF_bin), drop = TRUE)
q_binned <- split(q, list(gene_dist_bin, MAF_bin), drop = TRUE)

my_res <- mapply(my_GoShifter, my_p = p_binned, my_q = q_binned, SIMPLIFY = FALSE)

saveRDS(my_res, paste0("res/", my_block,".RDS"))
