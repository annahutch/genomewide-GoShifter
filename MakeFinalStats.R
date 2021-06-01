library(dplyr)
library(data.table)

# list all result files
files <- list.files(path = "res", pattern="*.RDS", full=TRUE)

# make one big list
res_list_of_lists <- lapply(files, readRDS)
res_list <- do.call(c, res_list_of_lists)

# first row of each list is for the actual SNP-annotation configuration
res_obs <- lapply(res_list, function(x) x[1,]) %>% rbindlist()

# test statistic for actual SNP-annotation configuration
X_obs <- (sum(res_obs$n1 * res_obs$p1) / sum(res_obs$n1)) / (sum(res_obs$n0 * res_obs$p0) / sum(res_obs$n0))

# remaining rows are statistics for random permutations of the annotation vector to generate a null distribution
res_null <- lapply(res_list, function(x) x[-1,])

# randomly select a row to sample many times and generate a null distribution
random_rows <- t(replicate(10000, lapply(res_null, function(x) x[sample.int(nrow(x), 1), ]) %>% rbindlist()))

# test statistics computed under the null
X_null <- apply(random_rows, 1, function(x) (sum(x$n1 * x$p1) / sum(x$n1)) / (sum(x$n0 * x$p0) / sum(x$n0)))

saveRDS(list(X_obs = X_obs, X_null = X_null), "final_test_statistics.RDS")