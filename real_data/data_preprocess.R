setwd("~/Programming/R/MuBaVaSe/real_data")
## load group data
dta_group1 <- read.csv("raw_data/GSE9899_clinical_anns.csv")
dta_group2 <- read.csv("raw_data/dta.csv")
dta_group2 <- dta_group2[, -1]
dta_group <- merge(dta_group1, dta_group2, by = c("AOCSID", "Primary.Site", "Type", "Subtype", "StageCode", "Consolidated.Grade"))
head(dta_group)
dta_group <- dta_group[-which(dta_group$group == "NC"), ]
table(dta_group$group)

## load raw data 
dta_raw <- readRDS("raw_data/dta_raw.rds")
dta_raw <- t(dta_raw)
colnames(dta_raw) <- dta_raw[1, ]
dta_raw <- dta_raw[-1, ]
colnames(dta_raw)[1:6] <- c("AOCSID", "Primary.Site", "Type", "Subtype", "StageCode", "Consolidated.Grade")
# remove ":"
for (iter in 2:6) {
  tmp <- dta_raw[, iter]
  tmp_tmp <- rep(NA, length(tmp))
  tmp <- strsplit(tmp, ":")
  for (i in seq_len(length(tmp))) {
    tmp_tmp[i] <- tmp[[i]][[2]]
  }
  tmp_tmp <- gsub(" ", "", tmp_tmp, fixed = TRUE)
  dta_raw[, iter] <- tmp_tmp
}

## merge data with group
dta <- merge(dta_group, dta_raw, by = c("AOCSID", "Primary.Site", "Type", "Subtype"))

## select genes
KEGG <- fgsea::gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2019_Human')
gene_list <- KEGG$Apoptosis
gene_list <- as.matrix(gene_list)
colnames(gene_list) <- "SYMBOL"

## transform gene names
gene_tmp <- colnames(dta)[-seq_len(11)]
gene_tmp <- AnnotationDbi::select(hgu133a2.db::hgu133a2.db, gene_tmp, c("SYMBOL"))

## merge genes
gene_merge <- merge(gene_tmp, gene_list, by = "SYMBOL")
gene_num <- as.matrix(table(gene_merge$SYMBOL))
gene_use <- names(gene_num[which(gene_num < 3), ])
gene_use <- as.matrix(gene_use)
colnames(gene_use) <- "SYMBOL"
gene_use <- merge(gene_merge, gene_use, by = "SYMBOL")

## dta used
dta_use <- dta[, gene_use$PROBEID]
colnames(dta_use) <- gene_use$SYMBOL
saveRDS(dta_use, "dta_use.rds")


