## check reproducibility differences

tmp <- readRDS("inst/estimation/loso-predictions/kde-Region1-2009-2010-loso-predictions.rds")
tmp2 <- readRDS("inst/estimation/loso-predictions-check/kde-Region1-2009-2010-loso-predictions.rds")

diff_idx <- which(tmp!=tmp2, arr.ind=TRUE)
tmp[diff_idx[1,1], diff_idx[1,2]]-tmp2[diff_idx[1,1], diff_idx[1,2]]
sum(exp(as.numeric(as.matrix(tmp[,c(-1, -2)])))-exp(as.numeric(as.matrix(tmp2[,c(-1, -2)]))))
