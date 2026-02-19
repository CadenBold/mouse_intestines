head(meta)
control_group <- "WSPDsa 1 Day"
test_group <- "FAWDsa 1 Day"

# Helper Function that looks at the control group and test group input
# and is able to extract the group and age parts of the string
parse_group_label <- function(label) {
  parts <- strsplit(label, " ")[[1]]
  list(
    group = parts[1],
    age   = tolower(paste(parts[2], parts[3]))
  )
}


ctrl <- parse_group_label(control_group)
test <- parse_group_label(test_group)

mask_control <-
  meta$Group == ctrl$group &
  meta$Age   == ctrl$age

mask_test <-
  meta$Group == test$group &
  meta$Age   == test$age

keep <- mask_control | mask_test

# Removing Sample 9
keep[9] <- FALSE

meta_test <- meta[keep,]

# Examining The sample groups
table(meta_test$Age)
meta_test

meta_test$Group <- relevel(
  factor(meta_test$Group),
  ref = ctrl$group
)


count_test <- count_matrix[,rownames(meta_test)]

dds <- DESeqDataSetFromMatrix(
  countData = count_test,
  colData   = meta_test,
  design    = ~ Group
)

#Potential Filtering Option
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

interaction_term <- setdiff(
  resultsNames(dds),
  "Intercept"
)


interaction_term
res <- results(dds)
head(res)

res <- results(dds, name = interaction_term)
res <- as.data.frame(res)
res$ENSEMBL <- rownames(res)

# keep only significant genes
sig <- res[!is.na(res$padj) & res$padj < 0.05, ]

# add absolute log2 fold change
sig$absLog2FC <- abs(sig$log2FoldChange)

# sort by absolute effect size
sig_sorted <- sig[order(-sig$absLog2FC), ]
sig_sorted$clean_ensembl <- str_remove(sig_sorted$ENSEMBL, "\\.\\d+$")
sig_sorted$SYMBOL <- ids$SYMBOL[match(sig_sorted$clean_ensembl, ids$clean_ensembl)]

# Only get mapped symbols
sig_mapped <- sig_sorted[!is.na(sig_sorted$SYMBOL) & sig_sorted$SYMBOL != "", ]

# look at top hits
head(sig_sorted, 10)

top20 <- head(sig_mapped, 20)

# clean Ensembl again just in case
top20$clean_ensembl <- str_remove(top20$ENSEMBL, "\\.\\d+$")

# map gene symbols
top20$SYMBOL <- ids$SYMBOL[match(top20$clean_ensembl, ids$clean_ensembl)]

# if some symbols are missing, fall back to Ensembl
top20$LABEL <- ifelse(is.na(top20$SYMBOL), top20$clean_ensembl, top20$SYMBOL)

# classify direction
top20$Direction <- ifelse(top20$log2FoldChange > 0, paste0("Upregulated in ",test_group), 
                          paste0("Downregulated in ", test_group))


ggplot(top20, aes(x = reorder(LABEL, absLog2FC),
                  y = absLog2FC,
                  fill = Direction)) +
  geom_col() +
  coord_flip() +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6))+
  labs(
    x = "Gene Name",
    y = "| log2(FoldChange) |",
    title = "Top 20 Log Fold Genes (Mapped Gene Symbols Only)",
    subtitle = paste0(control_group," vs ",test_group)
  ) +
  theme_classic(base_size = 13)

sig$clean_ensembl <- str_remove(sig$ENSEMBL, "\\.\\d+$")
sig$ENTREZID <- ids$ENTREZID[match(sig$clean_ensembl, ids$clean_ensembl)]

# remove genes without Entrez IDs
sig <- sig[!is.na(sig$ENTREZID), ]

sig_up   <- sig[sig$log2FoldChange > 0, ]
sig_down <- sig[sig$log2FoldChange < 0, ]

up_entrez   <- unique(sig_up$ENTREZID)
down_entrez <- unique(sig_down$ENTREZID)

# background universe
bg_entrez <- unique(ids$ENTREZID)
bg_entrez <- bg_entrez[!is.na(bg_entrez)]

go_up <- enrichGO(
  gene          = up_entrez,
  universe      = bg_entrez,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

go_down <- enrichGO(
  gene          = down_entrez,
  universe      = bg_entrez,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

go_up_df <- as.data.frame(go_up)

# sort by Count descending
go_up_df <- go_up_df[order(-go_up_df$Count), ]

# pick top N categories
topN <- 15
go_up_top <- go_up_df[1:topN, ]


go_up_top$Description <- str_wrap(go_up_top$Description, width = 40)

ggplot(go_up_top, aes(x = reorder(Description, Count), y = Count)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    x = "GO Term",
    y = "Gene Count",
    title = paste0("GO Comparison ",control_group," - Upregulated in ", test_group),
    subtitle=paste0(control_group, " vs ", test_group)
  ) +
  theme_minimal(base_size = 13)

# Wrap text for top 15 downregulated GO terms
go_down_df <- as.data.frame(go_down)
go_down_df <- go_down_df[order(-go_down_df$Count), ]
go_down_top <- go_down_df[1:15, ]

go_down_top$Description <- str_wrap(go_down_top$Description, width = 30)

ggplot(go_down_top, aes(x = reorder(Description, Count), y = Count)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    x = "GO Term",
    y = "Gene Count",
    title = paste0("GO Comparison ",control_group," - Downregulated in ", test_group),
    subtitle=paste0(control_group, " vs ", test_group)
  ) +
  theme_minimal(base_size = 13)