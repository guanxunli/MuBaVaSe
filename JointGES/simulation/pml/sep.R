library(pcalg)
library(graph)

set.seed(1)

args <- commandArgs()
load(args[4])
lambda <- as.numeric(args[5])
#lambda <- lambda^2

#do joint estimation given single data
ges.alg <- function(data){
	data <- do.call(rbind, data)
	l0score <- new("GaussL0penObsScore", data=data, lambda = lambda * log(ncol(data)), intercept = FALSE, use.cpp = TRUE)
	ges.fit <- ges(l0score)
	return(as(ges.fit$essgraph, "graphNEL"))
}

#run jobs on cluster
gesdag.list <- lapply(data.list, ges.alg)

#save to file
save(gesdag.list, dag.list, file=paste("sep/", strsplit(basename(args[4]), ".", fixed=TRUE)[[1]][1], "_lambda_", toString(lambda), ".rda", sep=""))
