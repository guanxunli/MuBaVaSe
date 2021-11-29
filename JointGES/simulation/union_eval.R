library(pcalg)
library(graph)

ges.shd <- function(prefix){
	load(paste(prefix, ".rda", sep=""))
	get.shd <- function(g1, g2.list){
		p <- length(g2.list[[1]]@edgeL)
		elist <- lapply(1:p, function(i) {
			t <- lapply(g2.list, function(g2) g2@edgeL[[i]]$edges)
			t <- unique(unlist(t))
			if(length(t) == 0) {return(integer(0))} else {return(t)}})
		names(elist) <- 1:p
		g2 <- new("graphNEL", node=sapply(1:p, toString), edgeL=elist, edgemode="directed")
		shd(g1, dag2cpdag(g2))
	}
	mean(sapply(1:length(dag.list), function(i) get.shd(gesdag.list[[i]], dag.list[[i]])))
}

#true positive rate for directed edges
ges.tprate <- function(prefix){
	load(paste(prefix, ".rda", sep=""))
	get.tprate <- function(g1, g2.list){
		p <- length(g2.list[[1]]@edgeL)
		elist <- lapply(1:p, function(i) {
			t <- lapply(g2.list, function(g2) g2@edgeL[[i]]$edges)
			t <- unique(unlist(t))
			if(length(t) == 0) {return(integer(0))} else {return(t)}})
		names(elist) <- 1:p
		g2 <- new("graphNEL", node=sapply(1:p, toString), edgeL=elist, edgemode="directed")
		g1 <- as(g1, "matrix") != 0
		g1 <- xor(g1 , g1 & t(g1))
		g2 <- as(g2, "matrix") != 0
		sum(g1 & g2) / sum(g2)
		#sum((g1 | t(g1)) & (g2 | t(g2))) / sum(g2 | t(g2))
	}
	mean(sapply(1:length(dag.list), function(i) get.tprate(gesdag.list[[i]], dag.list[[i]])))
}

ps <- c(100)
neighs <- c(2)
ns <- c(600)
#ks <- c(3, 4, 5, 6)
ks <- c(3, 5)
#lambdas <- c(0.3, 1, 2, 3, 5, 10)
lambdas <- c(1, 2, 3, 4, 5)

for(p in ps) for(n in ns) for(neigh in neighs){

	#get the prefix of data
	prefix <- paste("p_", toString(p), "_neigh_", toString(neigh), "_n_", toString(n), sep="")
	print(prefix)

	#SHD
	png(paste("figure/", prefix, "_union.shd.png", sep=""))
	par(mar=c(5,5.5,3,1))
	a <- sapply(lambdas, function(lambda) ges.shd(paste("pml/joint/", prefix, "_k_", toString(ks[1]), "_lambda_", toString(lambda), sep="")))
	print(a)
	plot(x=lambdas, y=a, ylim=c(0, 250), axes=FALSE, xlab="scaling constant", ylab="average SHD for CP-DAG", col="blue",
		 type="o", cex.lab=2, lwd=3, pch=NA)
	axis(1, las=1, cex.axis=1.4)
	axis(2, las=1, cex.axis=1.4)
	box()
	a <- sapply(lambdas, function(lambda) ges.shd(paste("pml/sep/", prefix, "_k_", toString(ks[1]), "_lambda_", toString(lambda), sep="")))
	lines(x=lambdas, y=a, col="green", type="o", lwd=3, pch=NA)
	print(a)
	a <- sapply(lambdas[2:length(lambdas)], function(lambda) ges.shd(paste("pml/joint/", prefix, "_k_", toString(ks[2]), "_lambda_", toString(lambda), sep="")))
	lines(x=lambdas[2:length(lambdas)], y=a, col="blue", type="o", lwd=3, lty=2, pch=NA)
	print(a)
	a <- sapply(lambdas, function(lambda) ges.shd(paste("pml/sep/", prefix, "_k_", toString(ks[2]), "_lambda_", toString(lambda), sep="")))
	lines(x=lambdas, y=a, col="green", type="o", lwd=3, lty=2, pch=NA)
	print(a)
	legend("topright", legend=c("joint GES k = 3", "GES k = 3", "joint GES k = 5", "GES k = 5"), cex=1.4, 
		col=c("blue", "green", "blue", "green"), lty=c(1,1,2,2), lwd=rep(3,4))
	dev.off()

	#true positve rate for directed edges
	png(paste("figure/", prefix, "_union.tprate.png", sep=""))
	par(mar=c(5,5.5,3,1))
	a <- sapply(lambdas, function(lambda) ges.tprate(paste("pml/joint/", prefix, "_k_", toString(ks[1]), "_lambda_", toString(lambda), sep="")))
	print(a)
	plot(x=lambdas, y=a, ylim=c(0, 1), axes=FALSE, xlab="scaling constant", ylab="true positve rate for directed edges", col="blue",
		 type="o", cex.lab=2, lwd=3, pch=NA)
	axis(1, las=1, cex.axis=1.4)
	axis(2, las=1, cex.axis=1.4)
	box()
	a <- sapply(lambdas, function(lambda) ges.tprate(paste("pml/sep/", prefix, "_k_", toString(ks[1]), "_lambda_", toString(lambda), sep="")))
	lines(x=lambdas, y=a, col="green", type="o", lwd=3, pch=NA)
	print(a)
	a <- sapply(lambdas[2:length(lambdas)], function(lambda) ges.tprate(paste("pml/joint/", prefix, "_k_", toString(ks[2]), "_lambda_", toString(lambda), sep="")))
	lines(x=lambdas[2:length(lambdas)], y=a, col="blue", type="o", lwd=3, lty=2, pch=NA)
	print(a)
	a <- sapply(lambdas, function(lambda) ges.tprate(paste("pml/sep/", prefix, "_k_", toString(ks[2]), "_lambda_", toString(lambda), sep="")))
	lines(x=lambdas, y=a, col="green", type="o", lwd=3, lty=2, pch=NA)
	print(a)
	legend("topright", legend=c("joint GES k = 3", "GES k = 3", "joint GES k = 5", "GES k = 5"), cex=1.4, 
		col=c("blue", "green", "blue", "green"), lty=c(1,1,2,2), lwd=rep(3,4))
	dev.off()
}
