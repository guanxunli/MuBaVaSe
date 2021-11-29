library(pcalg)
library(graph)
library(glmnet)

set.seed(1)

args <- commandArgs()
load(args[4])
lambda1 <- as.numeric(args[5])
lambda2 <- as.numeric(args[6])
load(paste("../pml/joint/", strsplit(basename(args[4]), ".rda", fixed=TRUE)[[1]][1], "_lambda_", toString(lambda1), ".rda", sep=""))

subset <- function(y, x, data){
	t <- rep(0, ncol(data))
	if(length(x) <= 1)
		t[x] <- 1
	else{
		model <- glmnet(as.matrix(data[,x]), data[,y], family="gaussian", intercept=FALSE, lambda=lambda2)
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

gesdag.list <- lapply(1:dagnum, function(i) ges.alg(data.list[[i]], gesdag.list[[i]]))
save(gesdag.list, dag.list, file=paste("joint_lambda_2/", strsplit(basename(args[4]), ".rda", fixed=TRUE)[[1]][1], "_lambda_", toString(lambda2), ".rda", sep=""))
