library(pcalg)
library(graph)

#load the ovarian cancer data
load("../data/ovarian.rda")
genenames <- colnames(data[[1]])
p <- length(genenames)

#load the prediction result
args <- commandArgs()
load(args[4])
cutoff <- 0.75
dag1 <- matrix(as.vector(gesdag.list[[1]]$max > cutoff), nrow=p, ncol=p)
dag2 <- matrix(as.vector(gesdag.list[[2]]$max > cutoff), nrow=p, ncol=p)
dag <- dag1 | dag2
mat <- dag | t(dag)

#get the statistics
conn <- sapply(1:nrow(mat), function(i) sum(mat[i, ]))
print(conn)
genenames[which(conn > 6)]
