setwd("~/Programming/R/MuBaVaSe/real_data")
## load group data
dta_group1 <- read.csv("raw_data/GSE9899_clinical_anns.csv")
dta_group2 <- read.csv("raw_data/dta.csv")
dta_group2 <- dta_group2[, -1]
dta_group <- merge(dta_group1, dta_group2, by = c("AOCSID", "Primary.Site", "Type", "Subtype", "StageCode", "Consolidated.Grade"))
head(dta_group)
dta_group <- dta_group[-which(dta_group$group == "NC"), ]
table(dta_group$group)
# generate data frame
df_group <- as.data.frame(matrix(NA, nrow = nrow(dta_group), ncol = ncol(dta_group)))
colnames(df_group) <- colnames(dta_group)
df_group[, 1:8] <- apply(dta_group, 2, as.character)

## load raw data 
dta_raw <- readRDS("raw_data/dta_raw.rds")
dta_raw <- t(dta_raw)
colnames(dta_raw) <- dta_raw[1, ]
dta_raw <- dta_raw[-1, ]
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

# generate data frame
df_raw <- as.data.frame(matrix(NA, nrow = nrow(dta_raw), ncol = ncol(dta_raw)))
df_raw[, 8:ncol(dta_raw)] <- apply(dta_raw[, 8:ncol(dta_raw)], 2, as.numeric)
df_raw[, 1:7] <- apply(dta_raw[, 1:7], 2, as.character)
colnames(df_raw)[1:6] <- c("AOCSID", "Primary.Site", "Type", "Subtype", "StageCode", "Consolidated.Grade")
colnames(df_raw)[7:ncol(df_raw)] <- colnames(dta_raw)[7:ncol(dta_raw)]
rm(dta_raw)
## merge data with group
dta <- merge(df_group, df_raw, by = c("AOCSID", "Primary.Site", "Type", "Subtype"))
rownames(dta) <- dta$AOCSID

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
gene_name <- unique(gene_use$SYMBOL)
n_gene <- length(gene_name)
dta_use <- matrix(NA, nrow = nrow(dta), ncol = n_gene)
rownames(dta_use) <- rownames(dta)
for (i in seq_len(n_gene)) {
  gene_tmp <- gene_name[1]
  index <- which(gene_use$SYMBOL == gene_tmp)
  probeid_tmp <- gene_use$PROBEID[index]
  dta_tmp <- matrix(as.numeric(dta[, probeid_tmp]), nrow = nrow(dta))
  dta_use[, i] <- rowSums(dta_tmp)
}
colnames(dta_use) <- gene_name
saveRDS(dta_use, "dta_use.rds")