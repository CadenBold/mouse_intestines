top_genes_all <- lapply(
  setdiff(resultsNames(dds), "Intercept"),
  function(coef_name) {
    
    res <- results(dds, name = coef_name) %>%
      as.data.frame() %>%
      rownames_to_column("ENSEMBL") %>%
      filter(!is.na(padj), padj < 0.05) %>%
      arrange(desc(abs(log2FoldChange))) %>%
      slice_head(n = 10) %>%
      mutate(term = coef_name)
    
    res
  }
)

top_genes_all <- bind_rows(top_genes_all)

print(top_genes_all)
top_genes_all$clean_ensembl <- str_remove(top_genes_all$ENSEMBL, "\\.\\d+$")
top_genes_all$SYMBOL <- ids$SYMBOL[match(top_genes_all$clean_ensembl, ids$clean_ensembl)]

top_genes_all <- top_genes_all[!is.na(top_genes_all$SYMBOL) & top_genes_all$SYMBOL != "", ]
print(top_genes_all$SYMBOL)

vst_sub <- vst_df[rownames(vst_df) %in% top_genes_all$ENSEMBL, , drop = FALSE]

rownames(vst_sub) <- top_genes_all$SYMBOL[
  match(rownames(vst_sub), top_genes_all$ENSEMBL)
]

vst_scaled <- t(scale(t(vst_sub)))

vst_long <- vst_sub %>%
  as.data.frame() %>%
  mutate(SYMBOL = rownames(.)) %>%
  pivot_longer(
    cols = -SYMBOL,
    names_to = "SampleID",
    values_to = "expression"
  )

head(meta)
meta$plot_group
# Join metadata (plot_group and groupID)

vst_long <- vst_long %>%
  mutate(SampleID_numeric = as.character(as.integer(str_remove(SampleID, "^Yim-"))))


vst_long <- vst_long %>%
  left_join(dplyr::select(meta, SampleID, plot_group, GroupID), 
            by = c("SampleID_numeric" = "SampleID"))


print(head(vst_long))

avg_df <- vst_long %>%
  group_by(SYMBOL, plot_group) %>%
  summarise(avg_expression = mean(expression), .groups = "drop")

avg_mat <- avg_df %>%
  pivot_wider(names_from = plot_group, values_from = avg_expression) %>%
  column_to_rownames("SYMBOL") %>%
  as.matrix()

avg_scaled <- t(scale(t(avg_mat)))


pheatmap(avg_scaled)