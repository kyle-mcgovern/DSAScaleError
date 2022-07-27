source("gsea.R")
library(gtools)

# GSEA-LFC-S

set.seed(3872)
m <- log(t(rdirichlet(6, c(1,1,1,1,1,1))))
set.seed(NULL)
lfc_1 <- apply(m, 1, function(row) mean(row[c(4,5,6)])-mean(row[c(1,2,3)]))
lfc_2 <- apply(m, 1, function(row) mean(row[c(2,5,6)])-mean(row[c(1,3,4)]))
names(lfc_1) <- c("A", "B", "C", "D", "E", "F")
names(lfc_2) <-  c("A", "B", "C", "D", "E", "F")
sorted_lfc_1 <- sort(lfc_1, decreasing=T)
sorted_lfc_2 <- sort(lfc_2, decreasing=T)
print(m)
print(lfc_1)
print(lfc_2)
print(sorted_lfc_1)
print(sorted_lfc_2)

get_lfc <- function(m) {
	apply(m, 1, function(row) mean(row[c(4,5,6)])-mean(row[c(1,2,3)]));
}

all_gsea_s <- function(m, inds, epsilon) {
	lfc <- get_lfc(m)
	combs <- combinations(length(lfc), length(inds), 1:length(lfc))
	scores <- c()
	obs_score <- gsea_base(lfc+epsilon, inds)
	print(obs_score)
	for(i in 1:nrow(combs)) {
		scores <- c(scores, gsea_base(lfc+epsilon, combs[i,]))
	}
	scores <- scores[sign(scores)==sign(obs_score)]
	sum((abs(scores) >= obs_score)) / length(scores)
}

m <- cbind(c(1,1,1,1,1,10),c(1,1,1,1,38.3,12),c(1,1,1,1,1,1),c(25,16,12,8,4,0.02),c(20,16,32,8,4,0.01),c(50,35,15,20,1,0.01))
m <- round(apply(m, 2, function(col) log(col/sum(col))),2)
m
lfc_1 <- apply(m, 1, function(row) mean(row[c(4,5,6)])-mean(row[c(1,2,3)]))
lfc_2 <- apply(m, 1, function(row) mean(row[c(2,5,6)])-mean(row[c(1,3,4)]))
names(lfc_1) <- 1:6
names(lfc_2) <- 1:6
lfc_1
lfc_2
sort(lfc_1, decreasing=T)
sort(lfc_2, decreasing=T)

for (e in seq(-5,5,0.1)) {
	print(c(e, all_gsea_s(m, c(1,2,3),e)))
}


# GSEA-LFC

all_gsea <- function(lfc, inds, epsilon) {
	combs <- combinations(length(lfc), length(inds), 1:length(lfc))
	scores <- c()
	obs_score <- gsea_base(lfc+epsilon, inds)
	for(i in 1:nrow(combs)) {
		scores <- c(scores, gsea_base(lfc+epsilon, combs[i,]))
	}
	scores <- scores[sign(scores)==sign(obs_score)]
	sum((abs(scores) >= obs_score)) / length(scores)
}

lfc <- c(5,4,3,2,1,0,-1)
for (e in seq(-5,5,0.1)) {
	print(c(e, all_gsea(lfc, c(1,2,3),e)))
}
