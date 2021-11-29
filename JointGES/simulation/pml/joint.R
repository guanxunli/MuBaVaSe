library(pcalg)
library(graph)
library(parallel)

set.seed(1)

#source("newclass.R")

args <- commandArgs()
load(args[4])
lambda <- as.numeric(args[5])

#do joint estimation given single data
ges.alg <- function(data){
	source("newclass.R")
	l0score <- new("MultiGaussL0pen", data = data, lambda = lambda * log(ncol(data[[1]])), intercept = FALSE, use.cpp = FALSE)
	ges.fit <- ges(l0score)
	return(as(ges.fit$essgraph, "graphNEL"))
}

#prepare for running parallel
#cor.num <- detectCores()
cl <- makeCluster(getOption("cl.cores", 20))
clusterExport(cl, c("lambda"))
clusterEvalQ(cl, c(library(pcalg), library(graph)))

#run jobs on cluster
gesdag.list <- parLapply(cl, data.list, ges.alg)

#stop the cluster after finishing jobs
stopCluster(cl)

#gesdag.list <- lapply(1:dagnum, function(i) ges.alg(data.list[[i]]))
save(gesdag.list, dag.list, file=paste("joint/", strsplit(basename(args[4]), ".rda", fixed=TRUE)[[1]][1], "_lambda_", toString(lambda), ".rda", sep=""))
