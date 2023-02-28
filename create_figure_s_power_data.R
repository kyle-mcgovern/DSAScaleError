#!/usr/bin/env Rscript

library(fgsea)
source("gsea.R")

set.seed(2002)

cores <- 39
niter <- 5000

pathways <- gmtPathways("c2.all.v7.4.symbols.gmt")
lfc_table <- read.table("./throid_songbird/differentials.tsv", header=T)
dim(lfc_table)
lfc <- lfc_table[,3]
names(lfc) <- lfc_table[,1]
head(lfc)

path_inds <- list()
path_names <- list()
for (p in names(pathways)) {
	inds <- match(pathways[[p]], names(lfc))
	inds <- inds[!is.na(inds)]
	if((length(inds)<15)|(length(inds)>500)) {
		next
	}
	path_inds[p] <- list(inds)
	path_names[[p]] <- pathways[[p]][pathways[[p]]%in%names(lfc)]
}

# Grid Search
epsilons <- c(seq(-200, -11, 1), seq(-10, 10, 0.1), seq(11,200,1))
signif_lst <- lfc_sensitivity_testing(lfc, path_inds, niter, epsilons, cores)
signif_lst

# Significant at values of tilde theta
m <- c()
for(e in seq(-2,2,0.2)) {
	print(e)
	signif <- gsea_parallel(lfc+e, path_inds, iters=niter, cores=cores)$p_values <= 0.05
	m <- cbind(m, signif)
}
m

m <- cbind(signif_lst, m)
row.names(m) <- names(path_inds)
write.table(m, "./data/figure_s3_table.txt")

