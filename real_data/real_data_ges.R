## load data
dta <- readRDS("dta_use.rds")
dta_group <- read.csv("raw_data/dta.csv")
dta_group <- dta_group[, -1]
rownames(dta_group) <- dta_group$AOCSID
sample_group <- dta_group[rownames(dta), 8]
sample_group <- ifelse(sample_group == "1", 1, 0)
dta_1 <- dta[which(sample_group == 1), ]
dta_2 <- dta[which(sample_group == 0), ]

## GES method
library(pcalg)
# data set 1
score_ges <- new("GaussL0penObsScore", data = dta_1, intercept = FALSE, 
                 lambda = sqrt(2 * log(ncol(dta_1)) / nrow(dta_1)))
ges_fit1 <- ges(score_ges)
ges_adj1 <- as(ges_fit1$repr, "matrix")
ges_adj1 <- ifelse(ges_adj1 == TRUE, 1, 0)
ges_adj1_u <- ceiling((ges_adj1 + t(ges_adj1)) / 2)
sum(ges_adj1_u) / 2 # 99
# data set 2
score_ges <- new("GaussL0penObsScore", data = dta_2, intercept = FALSE,
                 lambda = sqrt(2 * log(ncol(dta_2)) / nrow(dta_2)))
ges_fit2 <- ges(score_ges)
ges_adj2 <- as(ges_fit2$repr, "matrix")
ges_adj2 <- ifelse(ges_adj2 == TRUE, 1, 0)
ges_adj2_u <- ceiling((ges_adj2 + t(ges_adj2)) / 2)
sum(ges_adj2_u) / 2 # 197
# Intersection
length(intersect(which(ges_adj1_u == 1), which(ges_adj2_u == 1))) / 2 # 99

## get center
ges_adj <- ceiling((ges_adj1 + ges_adj2) / 2)
c_row <- rowSums(ges_adj)
c_col <- colSums(ges_adj)
c_sum <- c_row + c_col
names(c_sum) <- colnames(dta_1)
head(sort(c_sum, decreasing = TRUE))
