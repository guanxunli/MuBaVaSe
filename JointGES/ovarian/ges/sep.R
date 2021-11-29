library(pcalg)
library(graph)
library(stabs)

set.seed(1)

load("ovarian/data/ovarian.rda")
cutoff <- 0.6

#learn causal networks
stabs.ges <- function(x, y, q, ...){
	#Y is the label of the classes, X is the input matrix
	idx <- y[1]
	totcol <- nrow(data[[idx]])
	dt <- data[[idx]][sample(1:totcol, as.integer(0.9 * totcol), replace=FALSE), ]

	#train the model
	lambdas <- c(1,2,3,4,5)
	model.lambda <- function(lambda){
		l0score <- new("GaussL0penObsScore", data=dt, lambda = lambda * log(ncol(dt)), intercept = FALSE, use.cpp = TRUE)
		ges.fit <- ges(l0score)
		dag <- as(ges.fit$essgraph, "matrix")
		as.vector(dag != 0)
	}

	#get the path and selected variables
	path <- sapply(lambdas, model.lambda)
	selected <- rowSums(path) != 0
	return(list(selected=selected, path=path))
}

#run jobs on cluster
stab.inputlist <- lapply(1:length(data), function(i) {p <- ncol(data[[i]]); 
	dt <- cbind(as.matrix(data[[i]]), matrix(0, nrow=nrow(data[[i]]), ncol=p * (p-1))); 
	list(x=dt, y=rep(i, nrow(data[[i]])))
	})
gesdag.list <- lapply(stab.inputlist, function(stab.input) stabsel(x = stab.input$x, y = stab.input$y, fitfun = stabs.ges, cutoff = cutoff, PFER = 1))

#save to file
save(gesdag.list, file=paste("sep/ovarian_result.rda", sep=""))
