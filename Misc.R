terms <- setdiff(resultsNames(dds), "Intercept")

# Baseline term (to exclude shared genes)
res_base <- results(dds, name = terms[1]) %>%
  as.data.frame() %>%
  rownames_to_column("ENSEMBL")

exclude <- res_base %>%
  filter(!is.na(padj), padj < 0.05)

# Exposure term of interest
res_exp <- results(dds, name = terms[2]) %>%
  as.data.frame() %>%
  rownames_to_column("ENSEMBL")

sig <- res_exp %>%
  filter(!is.na(padj), padj < 0.05)

# Exposure-specific genes
exposure_genes <- sig %>%
  filter(!ENSEMBL %in% exclude$ENSEMBL) %>%
  arrange(padj)

# -----------------------------
# 2. Map ENSEMBL → SYMBOL
# -----------------------------
exposure_genes <- exposure_genes %>%
  mutate(
    clean_ensembl = str_remove(ENSEMBL, "\\.\\d+$"),
    SYMBOL = ids$SYMBOL[match(clean_ensembl, ids$clean_ensembl)]
  ) %>%
  filter(!is.na(SYMBOL), SYMBOL != "")

# -----------------------------
# 3. Extract VST expression for top genes
# -----------------------------
expr <- assay(vsd)
rownames(expr) <- str_remove(rownames(expr), "\\.\\d+$")

expr <- expr[
  ids$clean_ensembl[ids$SYMBOL %in% exposure_genes$SYMBOL],
  ,
  drop = FALSE
]

rownames(expr) <- ids$SYMBOL[
  match(rownames(expr), ids$clean_ensembl)
]



ggplot(expr_long_df, aes(x = plot_group, y = expression, color = plot_group)) +
  stat_summary(
    fun = mean,
    geom = "crossbar",
    width = 0.6,
    color = "black",
    linewidth = 0.5
  ) +
  geom_jitter(
    width = 0.15,
    size = 3,
    alpha = 0.8
  ) +
  facet_wrap(~ SYMBOL, scales = "free_y") +
  theme_classic(base_size = 13) +
  theme(
    panel.grid.major.y = element_line(color = "grey85", linewidth = 0.5),
    panel.grid.minor.y = element_line(color = "grey90", linewidth = 0.25),
    panel.grid.major.x = element_line(color = "grey88", linewidth = 0.4),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 1.09, hjust = 1),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = c(1, 0),      # bottom right
    legend.justification = c(1, 0),
    legend.background = element_rect(fill = "white", color = "black"),  # optional border
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 11)
  )  +
  labs(
    x = "Exposure / Age Group",
    y = "Expression (VST Normalized)",
    color = "Group",
    title = "Top 5 Western Diet genes with FA 1D ≈ FA 28D (p-value)"
  )