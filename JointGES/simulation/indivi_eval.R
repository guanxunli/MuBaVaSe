library(pcalg)
library(graph)

#SHD for GES
ges.shd <- function(prefix){
	load(paste(prefix, ".rda", sep=""))
	mean(sapply(1:length(dag.list), function(i) mean(sapply(1:length(dag.list[[i]]), function(t) shd(gesdag.list[[i]][[t]], dag2cpdag(dag.list[[i]][[t]]))))))
}

#true positive rate for directed edges
ges.tprate <- function(prefix){
	load(paste(prefix, ".rda", sep=""))
	tprate <- function(g1, g2){
		g1 <- as(g1, "matrix") != 0
		g2 <- as(g2, "matrix") != 0
		sum(g1 & g2) / sum(g2)
	}
	mean(sapply(1:length(dag.list), function(i) mean(sapply(1:length(dag.list[[i]]), function(j) tprate(gesdag.list[[i]][[j]], dag.list[[i]][[j]])))))
}

#false positive rate for directed edges
ges.fprate <- function(prefix){
	load(paste(prefix, ".rda", sep=""))
	fprate <- function(g1, g2){
		g1 <- as(g1, "matrix") != 0
		g2 <- as(g2, "matrix") != 0
		fp <- sum(g1 | t(g1)) / 2 - sum(g1 & g2)
		p <- dim(g1)
		fp / (p * (p - 1) / 2 - sum(g2))
	}
	mean(sapply(1:length(dag.list), function(i) mean(sapply(1:length(dag.list[[i]]), function(j) fprate(gesdag.list[[i]][[j]], dag.list[[i]][[j]])))))
}

ps <- c(100)
neighs <- c(2)
ns <- c(600)
#ks <- c(3, 4, 5, 6)
ks <- c(3, 5)
lambdas <- c(1, 2, 3, 4, 5)
#lambdas <- c(2, 3, 4, 5)

for(p in ps) for(n in ns) for(neigh in neighs){

	#get the prefix of data
	#prefix <- paste("p_", toString(p), "_neigh_", toString(neigh), "_n_", toString(n), "_k_", toString(k), sep="")
	prefix <- paste("p_", toString(p), "_neigh_", toString(neigh), "_n_", toString(n), sep="")
	print(prefix)

	shd.plot <- function(){
		#SHD
		png(paste("figure/", prefix, ".shd.png", sep=""))
		par(mar=c(5,5.5,3,1))
		a <- sapply(lambdas, function(lambda) ges.shd(paste("indivi/joint/", prefix, "_k_", toString(ks[1]), "_lambda_", toString(lambda), sep="")))
		print(a)
		plot(x=lambdas, y=a, ylim=c(0, 180), axes=FALSE, xlab="scaling constant (c)", ylab="structural Hamming distance", col="blue",
			 type="o", cex.lab=2, lwd=3, pch=NA)
		axis(1, las=1, cex.axis=1.4)
		axis(2, las=1, cex.axis=1.4)
		box()
		a <- sapply(lambdas, function(lambda) ges.shd(paste("indivi/sep/", prefix, "_k_", toString(ks[1]), "_lambda_", toString(lambda), sep="")))
		lines(x=lambdas, y=a, col="green", type="o", lwd=3, pch=NA)
		print(a)
		a <- sapply(lambdas, function(lambda) ges.shd(paste("indivi/joint/", prefix, "_k_", toString(ks[2]), "_lambda_", toString(lambda), sep="")))
		lines(x=lambdas, y=a, col="blue", type="o", lwd=3, lty=2, pch=NA)
		print(a)
		a <- sapply(lambdas, function(lambda) ges.shd(paste("indivi/sep/", prefix, "_k_", toString(ks[2]), "_lambda_", toString(lambda), sep="")))
		lines(x=lambdas, y=a, col="green", type="o", lwd=3, lty=2, pch=NA)
		print(a)
		legend("topright", legend=c("jointGES K = 3", "GES K = 3", "jointGES K = 5", "GES K = 5"), cex=1.4, 
			col=c("blue", "green", "blue", "green"), lty=c(1,1,2,2), lwd=rep(3,4))
		dev.off()
	}

	tprate.plot <- function(){
		#true positve rate for directed edges
		png(paste("figure/", prefix, ".tprate.png", sep=""))
		par(mar=c(5,5.5,3,1))
		a <- sapply(lambdas, function(lambda) ges.tprate(paste("indivi/joint/", prefix, "_k_", toString(ks[1]), "_lambda_", toString(lambda), sep="")))
		print(a)
		plot(x=lambdas, y=a, ylim=c(0, 1), axes=FALSE, xlab="scaling constant (c)", ylab="true positive rate", col="blue",
			 type="o", cex.lab=2, lwd=3, pch=NA)
		axis(1, las=1, cex.axis=1.4)
		axis(2, las=1, cex.axis=1.4)
		box()
		a <- sapply(lambdas, function(lambda) ges.tprate(paste("indivi/sep/", prefix, "_k_", toString(ks[1]), "_lambda_", toString(lambda), sep="")))
		lines(x=lambdas, y=a, col="green", type="o", lwd=3, pch=NA)
		print(a)
		a <- sapply(lambdas, function(lambda) ges.tprate(paste("indivi/joint/", prefix, "_k_", toString(ks[2]), "_lambda_", toString(lambda), sep="")))
		lines(x=lambdas, y=a, col="blue", type="o", lwd=3, lty=2, pch=NA)
		print(a)
		a <- sapply(lambdas, function(lambda) ges.tprate(paste("indivi/sep/", prefix, "_k_", toString(ks[2]), "_lambda_", toString(lambda), sep="")))
		lines(x=lambdas, y=a, col="green", type="o", lwd=3, lty=2, pch=NA)
		print(a)
		legend("bottomleft", legend=c("jointGES K = 3", "GES K = 3", "jointGES K = 5", "GES K = 5"), cex=1.4, 
			col=c("blue", "green", "blue", "green"), lty=c(1,1,2,2), lwd=rep(3,4))
		dev.off()
	}

	fprate.plot <- function(){
		#false positve rate for directed edges
		png(paste("figure/", prefix, ".fprate.png", sep=""))
		par(mar=c(5,5.5,3,1))
		a <- sapply(lambdas, function(lambda) ges.fprate(paste("indivi/joint/", prefix, "_k_", toString(ks[1]), "_lambda_", toString(lambda), sep="")))
		print(a)
		plot(x=lambdas, y=a, axes=FALSE, xlab="scaling constant", ylab="true positve rate for directed edges", col="blue",
			 type="o", cex.lab=2, lwd=3, pch=NA)
		axis(1, las=1, cex.axis=1.4)
		axis(2, las=1, cex.axis=1.4)
		box()
		a <- sapply(lambdas, function(lambda) ges.fprate(paste("indivi/sep/", prefix, "_k_", toString(ks[1]), "_lambda_", toString(lambda), sep="")))
		lines(x=lambdas, y=a, col="green", type="o", lwd=3, pch=NA)
		print(a)
		a <- sapply(lambdas, function(lambda) ges.fprate(paste("indivi/joint/", prefix, "_k_", toString(ks[2]), "_lambda_", toString(lambda), sep="")))
		lines(x=lambdas, y=a, col="blue", type="o", lwd=3, lty=2, pch=NA)
		print(a)
		a <- sapply(lambdas, function(lambda) ges.fprate(paste("indivi/sep/", prefix, "_k_", toString(ks[2]), "_lambda_", toString(lambda), sep="")))
		lines(x=lambdas, y=a, col="green", type="o", lwd=3, lty=2, pch=NA)
		print(a)
		legend("bottomleft", legend=c("joint GES k = 3", "GES k = 3", "joint GES k = 5", "GES k = 5"), cex=1.4, 
			col=c("blue", "green", "blue", "green"), lty=c(1,1,2,2), lwd=rep(3,4))
		dev.off()
	}

	#shd.plot()
	tprate.plot()
}

