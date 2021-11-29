library(pcalg)

#result of ges in terms of different lambda
ges.shd <- function(prefix){
	load(paste(prefix, ".rda", sep=""))
	mean(sapply(1:length(dag.list), function(i) mean(sapply(1:length(dag.list[[i]]), function(j) shd(gesdag.list[[i]][[j]], dag2cpdag(dag.list[[i]][[j]]))))))
}

ges.joint.shd <- function(prefix){
	load(paste(prefix, ".rda", sep=""))
	mean(sapply(1:length(dag.list), function(i) mean(sapply(1:length(dag.list[[i]]), function(j) shd(gesdag.list[[i]], dag2cpdag(dag.list[[i]][[j]]))))))
}

ges.tprate <- function(prefix){
	load(paste(prefix, ".rda", sep=""))
	tprate <- function(g1, g2){
		g1 <- as(g1, "matrix") != 0
		g1 <- xor(g1 , g1 & t(g1))
		g2 <- as(g2, "matrix") != 0
		sum(g1 & g2) / sum(g2)
		#sum((g1 | t(g1)) & (g2 | t(g2))) / sum(g2 | t(g2))
	}
	mean(sapply(1:length(dag.list), function(i) mean(sapply(1:length(dag.list[[i]]), function(j) tprate(gesdag.list[[i]][[j]], dag.list[[i]][[j]])))))
}

ges.joint.tprate <- function(prefix){
	load(paste(prefix, ".rda", sep=""))
	tprate <- function(g1, g2){
		g1 <- as(g1, "matrix") != 0
		g2 <- as(g2, "matrix") != 0
		sum(g1 & g2) / sum(g2)
		#sum((g1 | t(g1)) & (g2 | t(g2))) / sum(g2 | t(g2))
	}
	mean(sapply(1:length(dag.list), function(i) mean(sapply(1:length(dag.list[[i]]), function(j) tprate(gesdag.list[[i]], dag.list[[i]][[j]])))))
}

ges.skel <- function(prefix){
	load(paste(prefix, ".rda", sep=""))
	tprate <- function(g1, g2){
		g1 <- as(g1, "matrix") != 0
		g2 <- as(g2, "matrix") != 0
		sum((g1 | t(g1)) & (g2 | t(g2))) / sum(g2 | t(g2))
	}
	mean(sapply(1:length(dag.list), function(i) mean(sapply(1:length(dag.list[[i]]), function(j) tprate(gesdag.list[[i]][[j]], dag.list[[i]][[j]])))))
}

ges.joint.skel <- function(prefix){
	load(paste(prefix, ".rda", sep=""))
	tprate <- function(g1, g2){
		g1 <- as(g1, "matrix") != 0
		g2 <- as(g2, "matrix") != 0
		sum((g1 | t(g1)) & (g2 | t(g2))) / sum(g2 | t(g2))
	}
	mean(sapply(1:length(dag.list), function(i) mean(sapply(1:length(dag.list[[i]]), function(j) tprate(gesdag.list[[i]], dag.list[[i]][[j]])))))
}

ps <- c(100)
neighs <- c(2)
ns <- c(600)
ks <- c(3, 4, 5, 6)
lambdas <- c(0.3, 1, 2, 3, 5, 10)

for(p in ps) for(n in ns) for(neigh in neighs) for(k in ks){

	#get the prefix of data
	prefix <- paste("p_", toString(p), "_neigh_", toString(neigh), "_n_", toString(n), "_k_", toString(k), sep="")
	print(prefix)

	#shd
	png(paste("figure/", prefix, ".shd.png", sep=""))
	par(mar=c(5,5.5,3,1))
	#shd for normal ges
	a <- sapply(lambdas, function(lambda) ges.shd(paste("pml/sep/", prefix, "_lambda_", toString(lambda), sep="")))
	print(a)
	plot(x=lambdas, y=a, ylim=c(0, 200), axes=FALSE, xlab="constant", ylab="average SHD for CP-DAG", col="blue",
		 type="o", cex.lab=2, lwd=3, pch=NA)
	axis(1, las=1, cex.axis=1.4)
	axis(2, las=1, cex.axis=1.4)
	box()
	#shd for joint GES
	a <- sapply(lambdas, function(lambda) ges.joint.shd(paste("pml/joint/", prefix, "_lambda_", toString(lambda), sep="")))
	print(a)
	lines(x=lambdas, y=a, col="green", type="o", lwd=3, pch=NA)
	legend("topright", legend=c("separate GES", "joint GES"), col=c("blue", "green"), lty=c(1,1))
	dev.off()

	#true positive rate
	png(paste("figure/", prefix, ".tprate.png", sep=""))
	par(mar=c(5,5.5,3,1))
	#shd for normal ges
	a <- sapply(lambdas, function(lambda) ges.tprate(paste("pml/sep/", prefix, "_lambda_", toString(lambda), sep="")))
	plot(x=lambdas, y=a, ylim=c(0, 1), axes=FALSE, xlab="constant", ylab="true positive rate for CP-DAG", col="blue",
		 type="o", cex.lab=2, lwd=3, pch=NA)
	axis(1, las=1, cex.axis=1.4)
	axis(2, las=1, cex.axis=1.4)
	box()
	#shd for joint GES
	a <- sapply(lambdas, function(lambda) ges.joint.tprate(paste("pml/joint/", prefix, "_lambda_", toString(lambda), sep="")))
	lines(x=lambdas, y=a, col="green", type="o", lwd=3, pch=NA)
	legend("topright", legend=c("separate GES", "joint GES"), col=c("blue", "green"), lty=c(1,1))
	dev.off()

	#true positive rate skeleton
	png(paste("figure/", prefix, ".skel.png", sep=""))
	par(mar=c(5,5.5,3,1))
	#shd for normal GES
	a <- sapply(lambdas, function(lambda) ges.skel(paste("pml/sep/", prefix, "_lambda_", toString(lambda), sep="")))
	plot(x=lambdas, y=a, ylim=c(0, 1), axes=FALSE, xlab="constant", ylab="true positive rate for skeleton", col="blue",
		 type="o", cex.lab=2, lwd=3, pch=NA)
	axis(1, las=1, cex.axis=1.4)
	axis(2, las=1, cex.axis=1.4)
	box()
	#shd for joint GES
	a <- sapply(lambdas, function(lambda) ges.joint.skel(paste("pml/joint/", prefix, "_lambda_", toString(lambda), sep="")))
	lines(x=lambdas, y=a, col="green", type="o", lwd=3, pch=NA)
	legend("topright", legend=c("separate GES", "joint GES"), col=c("blue", "green"), lty=c(1,1))
	dev.off()
}
