library(pcalg)
library(graph)

set.seed(1)

args <- commandArgs()
load(args[4])
lambda <- as.numeric(args[5])
#lambda <- lambda^2

#do joint estimation given single data
ges.alg <- function(data){
	l0score <- lapply(data, function(x) new("GaussL0penObsScore", data=x, lambda = lambda * log(ncol(x)), intercept = FALSE, use.cpp = TRUE))
	ges.fit <- lapply(l0score, function(x) ges(x))
	return(lapply(ges.fit, function(x) as(x$essgraph, "graphNEL")))
}

#run jobs on cluster
gesdag.list <- lapply(data.list, ges.alg)

#save to file
save(gesdag.list, dag.list, file=paste("sep/", strsplit(basename(args[4]), ".rda", fixed=TRUE)[[1]][1], "_lambda_", toString(lambda), ".rda", sep=""))
