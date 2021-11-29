library(pcalg)
library(graph)

#load the ovarian cancer data
load("../data/ovarian.rda")
genenames <- colnames(data[[1]])
p <- length(genenames)

#load the prediction result
args <- commandArgs()
load(args[4])
#cutoff <- as.numeric(args[5])
cutoff <- 0.6
dag <- matrix(as.vector(stab.result$max > cutoff), nrow=p, ncol=p)

mat <- as(dag, "matrix")
mat <- mat | t(mat)
conn <- sapply(1:nrow(mat), function(i) sum(mat[i, ]))
print(conn)
genenames[which(conn > 5)]
