library(pcalg)
library(graph)

load("data.rda")
load("main.rda")

pardag2dag <- function(graph){
	#construct dag for PAG of joint graph
	p <- length(graph)
	edl <- lapply(1:p, function(x) list())
	names(edl) <- 1:p
	#convert graph connection from "in-edge" to "out-edge"
	for(i in 1:p) {
		in.edge.num <- length(graph[[i]]$edges)
		if(in.edge.num > 0)
			for(j in 1:in.edge.num) {
				node <- graph[[i]]$edges[j]
				if(length(edl[[node]]) == 0) 
					edl[[node]] <- list(edges=i, weights=graph[[i]]$weights[j])
				else
					edl[[node]] <- list(edges=c(edl[[node]]$edges, i), weights=c(edl[[node]]$weights, graph[[i]]$weights[j]))
			}
	}

	#get the new dag estimate from the joint estimation method
	return(new("graphNEL", node=sapply(1:p, toString), edgeL=edl, edgemode="directed"))
}

#convert graph learned in "main.rda" into graphnel object
joint_dag <- lapply(joint_graph, pardag2dag)
dagl <- lapply(indivi_graph, pardag2dag)

a <- sapply(1:length(joint_dag), function(x) c(shd(joint_dag[[x]], grs[[x]]), shd(dagl[[x]], grs[[x]])))
print(a)
print(rowMeans(a))
