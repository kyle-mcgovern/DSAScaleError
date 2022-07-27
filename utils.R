rand_strs <- function(n=5000) {
	a <- do.call(paste0, replicate(7, sample(LETTERS, n, TRUE), FALSE))
	paste0(a, sprintf("%04d", sample(9999, n, TRUE)), sample(LETTERS, n, TRUE))
}
