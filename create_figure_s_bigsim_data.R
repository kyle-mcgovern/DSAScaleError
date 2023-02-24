source("gsea.R")

Ds <- c(15, 50,  150)
D <- c(500, 10000, 20000)
epsilons <- seq(-2,2,0.2)
dist_types <- c("norm", "uniform", "right_skew", "left_skew")

set.seed(500)

m <- c()
for(ds in Ds) {
	for(d in D) {
		for (dist_type in dist_types) {
			if (dist_type=="norm") {
				lfc <- rnorm(d)
				lfc <- lfc - mean(lfc)
			} else if (dist_type=="uniform") {
				lfc <- runif(d, -1, 1)
				lfc <- lfc - mean(lfc)
			} else if (dist_type=="right_skew") {
				lfc <- runif(d*0.1, 0, 1)
				lfc <- c(lfc, rep(-0.1, d*0.9))
				lfc <- lfc - mean(lfc)
			} else if (dist_type=="left_skew") {
				lfc <- runif(d*0.1, -1, 0)
				lfc <- c(lfc, rep(0.1, d*0.9))
				lfc <- lfc - mean(lfc)
			}
			rds_file <- paste0("./data/", ds, "_", d, "_", dist_type, ".RDS")
			saveRDS(lfc, rds_file)

			null_inds <- replicate(5000, sample(1:length(lfc), ds, replace=F))
			test_inds <- replicate(10000, sample(1:length(lfc), ds, replace=F))
			cl <- makeCluster(7)
			registerDoParallel(cl)
			clusterCall(cl, function() { source("gsea_C.R") })
			e_m <- foreach(epsilon=epsilons, .combine=rbind, .noexport=c("gsea_score_C")) %dopar% {
				curr_m <- c()
				null_dist <- c()
				for(i in 1:ncol(null_inds)) {
					null_dist <- c(null_dist, gsea_base(lfc+epsilon, null_inds[,i]))
				}
				for(i in 1:ncol(test_inds)) {
					res <- gsea(lfc+epsilon, test_inds[,i], null_scores=null_dist)
					curr_m <- rbind(curr_m, c(dist_type, ds, d, epsilon, i, res$p_value, res$obs_score))
				}
				curr_m
			}
			stopCluster(cl)
			m <- rbind(m, e_m)
			print(dim(m))
		}
	}
}

colnames(m) <- c("dist_type", "ds", "d", "epsilon", "index", "p_value", "obs_score")
write.table(m, "./data/figure_s1_data.txt")
