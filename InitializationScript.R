library(limma)
library(DESeq2)
library(ggplot2)
library(tidyverse)
library(readr)
library(ggrepel)
library(clusterProfiler)
library(org.Mm.eg.db)   # change to org.Hs.eg.db if human
library(enrichplot)
library(DOSE)
library(readxl)
library(AnnotationDbi)
library(stringr)

FastP_featurecounts <- read.csv("FastP_featurecounts.csv")
metadata <- read.csv("metadata.csv")

#head(FastP_featurecounts)

fc <- FastP_featurecounts[,-(2:6)]
colnames(fc) <- as.character(unlist(fc[1,]))
fc <- fc[-1,]
rownames(fc) <- as.character(unlist(fc[,1]))
fc <- fc[,-1]
#head(fc)

count_matrix <- as.matrix(fc)
head(count_matrix)
storage.mode(count_matrix) <- "numeric"
print(is.numeric(count_matrix))

ids <- data.frame(ensembl_ids = rownames(fc), stringsAsFactors = FALSE)

# 2. Clean Ensembl IDs (remove .version)
ids$clean_ensembl <- str_remove(ids$ensembl_ids, "\\.\\d+$")

# 3. Query annotations
annotations <- select(org.Mm.eg.db,
                      keys = unique(ids$clean_ensembl),
                      columns = c("SYMBOL", "GENENAME", "ENTREZID"),
                      keytype = "ENSEMBL")

# 4. Enforce one-to-one mapping
annot_unique <- annotations[!duplicated(annotations$ENSEMBL), ]
rownames(annot_unique) <- annot_unique$ENSEMBL

# 5. Build master table by matching CLEAN IDs
ids$SYMBOL   <- annot_unique[ids$clean_ensembl, "SYMBOL"]
ids$GENENAME <- annot_unique[ids$clean_ensembl, "GENENAME"]
ids$ENTREZID <- annot_unique[ids$clean_ensembl, "ENTREZID"]

# 6. Inspect
head(ids)

# Check head of read in files, and then ccreate copies ot modify
head(metadata)
meta <- metadata

# Correct unique grouping so that we can quickly split on it
meta <- meta %>%
  mutate(
    GroupID = group_indices(., Age, Exposure, Diet, Treatment)
  )

meta %>%
  distinct(Age, Exposure, Diet, Treatment, GroupID) %>%
  arrange(GroupID)

meta <- as.data.frame(meta)
meta$SampleID <- sprintf("Yim-%03d", as.integer(rownames(meta)))
rownames(meta) <- meta$SampleID
meta$SampleID <- NULL
