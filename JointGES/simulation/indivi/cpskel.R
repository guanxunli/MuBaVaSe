library(pcalg)
library(graph)

set.seed(1)

source("newclass.R")

args <- commandArgs()
load(args[4])
lambda <- as.numeric(args[5])
load(paste("joint/", strsplit(basename(args[4]), ".", fixed=TRUE)[[1]][1], "_lambda_", toString(lambda), ".rda", sep=""))
#lambda <- lambda^2

#do joint estimation given single data
ges.alg <- function(data, dag){
	gaps <- as(dag, "matrix") != 0
	gaps <- gaps & t(gaps)
	joint.l0score <- lapply(data, function(x) new("GaussL0penObsScore", data=x, lambda = lambda * log(ncol(x)), intercept = FALSE, use.cpp = TRUE))
	joint.fit <- lapply(joint.l0score, function(x) ges(x, fixedGaps=gaps))
	return(lapply(joint.fit, function(fit) as(fit$essgraph, "graphNEL")))
}

gesdag.list <- lapply(1:dagnum, function(i) ges.alg(data.list[[i]], gesdag.list[[i]]))
save(gesdag.list, dag.list, file=paste("joint/", strsplit(basename(args[4]), ".", fixed=TRUE)[[1]][1], "_lambda_", toString(lambda), "_final", ".rda", sep=""))
