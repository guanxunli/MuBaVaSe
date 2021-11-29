library(graph)
library(pcalg)

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

#plot new figures
get.tpplot <- function(tp.plot, prefix, title){
	maxy = max(sapply(tp.plot, function(t) max(t$tp)))
	maxx = max(sapply(tp.plot, function(t) max(t$fp)))
	png(paste("figure/", prefix, "_", title, ".png", sep=""))
	par(mar=c(5,5.5,3,1))
	plot(0,0, xlim=c(0, maxx), ylim=c(0, maxy), axes=FALSE,  xlab="Number of false positives", ylab="Number of true positives", cex.lab=2, type="n", cex=2)
	axis(1, las=1, cex.axis=1.4)
	axis(2, las=1, cex.axis=1.4)
	box()
	for(i in 1:length(tp.plot))
		lines(x=tp.plot[[i]]$fp, y=tp.plot[[i]]$tp, col=cols[i], type="p", pch=pchs[i], cex=2)
	legend("bottomright", legend=c("joint GES K = 3", "GES K = 3", "joint GES K = 5", "GES K = 5"), col=cols, pch=pchs, cex=1.9)
	dev.off()
}

ps <- c(100)
ns <- c(600)
neighs <- c(2)
ks <- c(3, 5)
cols <- c("blue", "green", "blue", "green")
pchs <- c(1, 1, 20, 20)

for(p in ps) for(n in ns) for(neigh in neighs){

	#get the prefix of data
	prefix <- paste("p_", toString(p), "_neigh_", toString(neigh), "_n_", toString(n), sep="")
	print(prefix)

	#set up the significance levels
	#get results for joint GES
	lambdas <- c(0.00001, 0.00005, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5)
	jointk3.list <- lapply(lambdas, function(lambda){
		tp <- sapply(lambdas, function(lambda) ges.tprate(paste("indivi/joint/", prefix, "_k_", toString(ks[1]), "_lambda_", toString(lambda), sep="")))
		fp <- sapply(lambdas, function(lambda) ges.fprate(paste("indivi/joint/", prefix, "_k_", toString(ks[1]), "_lambda_", toString(lambda), sep="")))
		list(tp=tp, fp=fp)
	})

	jointk5.list <- lapply(lambdas, function(lambda){
		tp <- sapply(lambdas, function(lambda) ges.tprate(paste("indivi/joint/", prefix, "_k_", toString(ks[2]), "_lambda_", toString(lambda), sep="")))
		fp <- sapply(lambdas, function(lambda) ges.fprate(paste("indivi/joint/", prefix, "_k_", toString(ks[2]), "_lambda_", toString(lambda), sep="")))
		list(tp=tp, fp=fp)
	})
	#get results for separate GES
	lambdas <- c(0.25, 0.5, 1, 3, 5, 7, 9, 11, 15)
	sepk3.list <- lapply(lambdas, function(lambda){
		tp <- sapply(lambdas, function(lambda) ges.tprate(paste("indivi/sep/", prefix, "_k_", toString(ks[1]), "_lambda_", toString(lambda), sep="")))
		fp <- sapply(lambdas, function(lambda) ges.fprate(paste("indivi/sep/", prefix, "_k_", toString(ks[1]), "_lambda_", toString(lambda), sep="")))
		list(tp=tp, fp=fp)
	})
	sepk5.list <- lapply(lambdas, function(lambda){
		tp <- sapply(lambdas, function(lambda) ges.tprate(paste("indivi/sep/", prefix, "_k_", toString(ks[2]), "_lambda_", toString(lambda), sep="")))
		fp <- sapply(lambdas, function(lambda) ges.fprate(paste("indivi/sep/", prefix, "_k_", toString(ks[2]), "_lambda_", toString(lambda), sep="")))
		list(tp=tp, fp=fp)
	})

	#plot the results, second comes the rates for skeleton
	di.tpplot <- lapply(list(jointk3.list, sepk3.list, jointk5.list, sepk3.list,), function(tp.list) 
		list(tp=sapply(tp.list, function(t) t$tp), fp=sapply(tp.list, function(t) t$fp)))
	get.tpplot(di.tpplot, prefix, "directed")
}
