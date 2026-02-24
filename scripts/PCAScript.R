# Inputs to change entire script
graph_title <- "PCA of Samples without Sample 9 and 48"
pc <- "PC1"
keep<- meta$Liver.. != 9 & meta$Liver.. != 48

print(keep)

meta_sub <- meta[keep, ]
vst_sub <- vst[, rownames(meta_sub)]

gene_vars <- apply(vst_sub, 1, var)

top_n <- 1000
top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:top_n]

vst_top <- vst_sub[top_genes, ]
pca_input <- t(vst_top)

pca <- prcomp(pca_input, center = TRUE, scale. = TRUE)

pca_df <- data.frame(
  Sample = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  meta_sub[rownames(pca$x), ]
)

pve <- (pca$sdev^2) / sum(pca$sdev^2)

pca_df$GroupFactor <- interaction(
  pca_df$Age,
  pca_df$Exposure,
  pca_df$Diet,
  pca_df$Treatment,
  drop = TRUE
)

pca_df$Outlier <- with(
  pca_df,
  abs(scale(PC1)) > 2 | abs(scale(PC2)) > 2
)



ggplot(pca_df, aes(PC1, PC2)) +
  geom_point(
    aes(fill = GroupFactor),
    shape = 21,
    size = 7,
    color = "black",
    alpha = 0.9
  ) +
  geom_text(
    aes(label = str_remove(Sample,"^Yim-")),
    color = "Black",
    size = 3,
    fontface = "bold"
  ) +
  labs(
    x = paste0("PC1 (", round(pve[1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(pve[2] * 100, 1), "%)"),
    fill = "Group",
    title = graph_title
  ) +
  theme_classic()


pc1_loadings <- data.frame(
  Gene = colnames(pca_input),
  Loading = pca$rotation[,pc]
)

pc1_loadings$clean_ensembl <- str_remove(pc1_loadings$Gene, "\\.\\d+$")

# Map SYMBOLs
pc1_loadings$SYMBOL <- ids$SYMBOL[match(pc1_loadings$clean_ensembl, ids$clean_ensembl)]
pc1_loadings$AbsLoading <- abs(pc1_loadings$Loading)

pc1_negative <- subset(pc1_loadings, Loading > 0)
pc1_negative <- pc1_negative[order(pc1_negative$AbsLoading, decreasing=TRUE), ]

#pc1_negative <- head(pc1_negative, 500)


pc1_negative$clean_ensembl <- str_remove(rownames(pc1_negative),"\\.\\d+$")
pc1_negative$SYMBOL <- ids$SYMBOL[match(pc1_negative$clean_ensembl,ids$clean_ensembl)]
head(pc1_negative$SYMBOL)

pc1_negative$ENTREZID <- ids$ENTREZID[match(pc1_negative$clean_ensembl,ids$clean_ensembl)]
pc1_negative <- subset(pc1_negative, !is.na(ENTREZID))
length(unique(pc1_negative$ENTREZID))

background_genes <- data.frame(clean_ensembl = str_remove(rownames(vst_sub),"\\.\\d+$"))
background_genes$ENTREZID <- ids$ENTREZID[
  match(background_genes$clean_ensembl, ids$clean_ensembl)
]

background_entrez <- unique(na.omit(background_genes$ENTREZID))


go_pc <- enrichGO(
  gene = unique(pc1_negative$ENTREZID),
  universe = background_entrez,
  OrgDb         = org.Mm.eg.db,   # org.Hs.eg.db if human
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

go_pc
barplot(go_pc, showCategory = 10) + ggtitle(paste0(pc, " Loadings"))
head(pc1_negative$SYMBOL,20)



