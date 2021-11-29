library(pcalg)
library(graph)
library(glmnet)

set.seed(3)

args <- commandArgs()
load(args[4])
lambda <- as.numeric(args[5])
load(paste("../pml/joint/", strsplit(basename(args[4]), ".rda", fixed=TRUE)[[1]][1], "_lambda_", toString(lambda), ".rda", sep=""))

subset <- function(y, x, data){
	t <- rep(0, ncol(data))
	if(length(x) <= 1)
		t[x] <- 1
	else{
		model <- cv.glmnet(as.matrix(data[,x]), data[,y], family="gaussian", intercept=FALSE)
		nonz <- which(as.vector(coef(model)) != 0) - 1
		t[x[nonz]] <- 1
	}
	return(t)
}

#do joint estimation given single data
ges.alg <- function(data, dag){
	in.mat <- as(pdag2dag(dag)$graph, "matrix")
	joint.mat <- lapply(data, function(dt) sapply(1:ncol(dt), function(i) subset(i, which(in.mat[,i] != 0), dt)))
	return(lapply(joint.mat, function(sing.mat) dag2cpdag(as(sing.mat, "graphNEL"))))
}

gesdag.list <- lapply(1:dagnum, function(i) try(ges.alg(data.list[[i]], gesdag.list[[i]])))
idxs <- sapply(1:dagnum, function(i) !is.character(gesdag.list))
gesdag.list <- gesdag.list[idxs]
dag.list <- dag.list[idxs]
save(gesdag.list, dag.list, file=paste("joint/", strsplit(basename(args[4]), ".rda", fixed=TRUE)[[1]][1], "_lambda_", toString(lambda), ".rda", sep=""))
