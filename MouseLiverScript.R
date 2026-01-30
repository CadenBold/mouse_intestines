head(meta)

meta_test <- meta[
  meta$Age == "1 day" &
    meta$Treatment == "sa" &
    meta$Diet %in% c("PD", "WD") &
    meta$Exposure == "FA",
]

meta_test <- meta_test[-c(5),]

table(meta_test$Age)
meta_test

count_test <- count_matrix[,rownames(meta_test)]

dds <- DESeqDataSetFromMatrix(
  countData = count_test,
  colData   = meta_test,
  design    = ~ Group
)

#Potential Filtering Option
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

resultsNames(dds)

res <- results(dds)
head(res)

res <- results(dds, name = "Group_FAWDsa_vs_FAPDsa")
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
top20$Direction <- ifelse(top20$log2FoldChange > 0, "Upregulated in FAWDsa 1 Day", 
                          "Downregulated in FAWDsa 1 Day")


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
    subtitle = "FAWDsa 1 Day vs True Control (FAPDsa 1 Day) with Sample 9 Removed"
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
    title = "Top 15 GO Terms – Upregulated in FAWDsa 1 Day",
    subtitle="Sample 9 Removed from FAWDsa 1 Day"
  ) +
  theme_minimal(base_size = 13)

# Wrap text for top 15 downregulated GO terms
go_down_df <- as.data.frame(go_down)
go_down_df <- go_down_df[order(-go_down_df$Count), ]
go_down_top <- go_down_df[1:15, ]

go_down_top$Description <- str_wrap(go_down_top$Description, width = 40)

ggplot(go_down_top, aes(x = reorder(Description, Count), y = Count)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    x = "GO Term",
    y = "Gene Count",
    title = "Top 15 GO Terms – Downregulated in FAWDsa 1 Day",
    subtitle="Sample 9 Removed from FAWDsa 1 Day"
  ) +
  theme_minimal(base_size = 13)