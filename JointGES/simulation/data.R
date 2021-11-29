library(pcalg)
library(graph)
library(MASS)

#parameter used for generating data
ps <- c(100)
neighs <- c(2)
ks <- c(3, 5)
ns <- c(600, 900, 1200)
dagnum <- 100
uniq_edg <- 30 # set it to 30 if you want Fig 2. If you want Fig 3, set it to 60.

#set random seed
set.seed(1)

#generate random graph for all K class
randgraph <- function(p, neigh) {
	edl <- randomDAG(p, prob = neigh / p)
	nedl <- lapply(edl@edgeL, function(x) x$edges)
	return(nedl)
}

#generate random graph for each specific K
randngraph <- function(p, edl, uniq_edg, k){
	#get the list of unselected edges
	fedl <- combn(1:p, 2, simplify=FALSE)
	idx <- sapply(fedl, function(tt) !(tt[2] %in% edl[[tt[1]]]))
	fedl <- fedl[idx]

	#random select a list of edges
	fedl <- sample(fedl, uniq_edg)

	#random assign selected into K DAGs
	didx <- sample(1:k, length(fedl), replace=TRUE)
	kedl <- lapply(1:k, function(t) edl)
	for(idx in 1:length(fedl)){
		kedl[[didx[idx]]][[fedl[[idx]][1]]] <- c(kedl[[didx[idx]]][[fedl[[idx]][1]]], fedl[[idx]][2])
	}
	kedl <- lapply(kedl, function(kkedl) lapply(1:p, function(x) list(edges=kkedl[[x]], weights=sample(c(1,-1), length(kkedl[[x]]), replace=TRUE) * runif(length(kkedl[[x]]), min=0.1, max=1))))
	for(t in 1:k) 
		names(kedl[[t]]) <- 1:p
	return(kedl)
}

#generate random graphs
for(p in ps){
for(neigh in neighs){

	#get list of intersected dags
	ori.dag.list <- lapply(1:dagnum, function(i) randgraph(p, neigh))
	permut <- lapply(1:dagnum, function(i) sample(p, p))
	for(n in ns){
	for(k in ks){
		#get all K dags
		dag.list <- lapply(ori.dag.list, function(dag) {kedgeL <- randngraph(p, dag, uniq_edg, k); lapply(1:k, function(t) { new("graphNEL", node=sapply(1:p, toString), 
			edgeL=kedgeL[[t]], edgemode="directed")})})

		print("done with DAGs")
		#get data for all K dags
		data.list <- lapply(dag.list, function(dags) lapply(dags, function(dag) mvrnorm(as.integer(n / k), mu=rep(0,p), Sigma=trueCov(dag))))

		print("done with data")
		#random permutation of nodes as well as data for each node
		dag.list <- lapply(1:dagnum, function(i) lapply(dag.list[[i]], function(dag) {
			dagmat <- as(dag, "matrix")
			dagmat <- dagmat[permut[[i]], permut[[i]]]
			rownames(dagmat) <- 1:p
			colnames(dagmat) <- 1:p
			return(as(abs(dagmat), "graphNEL"))}))
		data.list <- lapply(1:dagnum, function(i) lapply(data.list[[i]], function(data) {
			t <- data[, permut[[i]]]
			colnames(t) <- 1:p
			return(t)}))

		#get title
		title <- paste("data/p_", toString(p), "_neigh_", toString(neigh), "_n_", toString(n), "_k_", toString(k), ".rda", sep="")
		#save dag and data information
		save(data.list, dag.list, dagnum, title, file=title)
	}}
}}
