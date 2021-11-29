library(pcalg)
library(graph)
library(stabs)

set.seed(1)

args <- commandArgs()
cutoff <- as.numeric(args[4])
print(cutoff)
source("newclass.R")

#load the drosophila data
load("../data/ovarian.rda")
p <- ncol(data[[1]])

#learn causal networks
stabs.ges <- function(x, y, q, ...){
	#Y is the label of the classes, X is the input matrix
	dt <- lapply(data, function(sing.dt) {totcol <- nrow(sing.dt); sing.dt[sample(1:totcol, as.integer(0.9 * totcol), replace=FALSE), ]})
	#dt <- lapply(1:max(y), function(i) x[y==i, 1:p])
	print(ncol(x))
	lambdas <- c(2,3,4,5)
	model.lambda <- function(lambda){
		l0score <- new("MultiGaussL0pen", data = dt, lambda = lambda * log(ncol(dt[[1]])), intercept = TRUE, use.cpp = FALSE)
		ges.fit <- ges(l0score)
		dag <- as(ges.fit$essgraph, "matrix")
		as.vector(dag != 0)
	}
	path <- sapply(lambdas, model.lambda)
	selected <- rowSums(path) != 0
	return(list(selected=selected, path=path))
}

#construct x
x <- do.call(rbind, data)
x <- cbind(x, matrix(0, nrow=nrow(x), ncol=p * (p-1)))

#construct y
y <- c()
for(i in 1:length(data)){
	y <- c(y, rep(i, nrow(data[[i]])))
}

stab.result <- stabsel(x = x, y = y, fitfun = stabs.ges, cutoff = cutoff, PFER = 1)
save(stab.result, file=paste("sep/ovarian_result.rda", sep=""))
