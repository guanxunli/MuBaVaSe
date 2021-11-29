library(pcalg)
library(graph)
library(leaps)

set.seed(1)

args <- commandArgs()
load(args[4])
lambda <- as.numeric(args[5])
load(paste("joint/", strsplit(basename(args[4]), ".", fixed=TRUE)[[1]][1], "_lambda_", toString(lambda), ".rda", sep=""))
penal <- lambda * log(ncol(data.list[[1]][[1]])) / sum(sapply(data.list[[1]], nrow))

subset <- function(y, x, data){
	t <- rep(0, ncol(data))
	if(length(x) == 0)
		return(t)
	else{
		model <- leaps(as.matrix(data[,x]), data[,y], names=x, method="r2", nbest=1, int=FALSE)
		model.optidx <- which.min(c(0, log(1 - model$r2) + penal * model$size))
		if(model.optidx == 1)
			return(t)
		else{
			t[x[model$which[model.optidx-1,]]] <- 1
			return(t)
		}
	}
}

#do joint estimation given single data
ges.alg <- function(data, dag){
	in.mat <- as(pdag2dag(dag)$graph, "matrix")
	joint.mat <- lapply(data, function(dt) sapply(1:ncol(dt), function(i) subset(i, which(in.mat[,i] != 0), dt)))
	return(lapply(joint.mat, function(sing.mat) dag2cpdag(as(sing.mat, "graphNEL"))))
}

gesdag.list <- lapply(1:dagnum, function(i) ges.alg(data.list[[i]], gesdag.list[[i]]))
save(gesdag.list, dag.list, file=paste("joint/", strsplit(basename(args[4]), ".", fixed=TRUE)[[1]][1], "_lambda_", toString(lambda), "_final", ".rda", sep=""))
