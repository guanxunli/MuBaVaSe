library(pcalg)
library(graph)
library(leaps)

set.seed(1)

args <- commandArgs()
load(args[4])
lambda <- as.numeric(args[5])
load(paste("../pml/joint/", strsplit(basename(args[4]), ".", fixed=TRUE)[[1]][1], "_lambda_", toString(lambda), ".rda", sep=""))

#do joint estimation given single data
ges.alg <- function(data, dag){
	k <- length(data)
	return(lapply(1:k, function(t) dag))
}

gesdag.list <- lapply(1:dagnum, function(i) ges.alg(data.list[[i]], gesdag.list[[i]]))
save(gesdag.list, dag.list, file=paste("joint/", strsplit(basename(args[4]), ".", fixed=TRUE)[[1]][1], "_lambda_", toString(lambda), ".rda", sep=""))
