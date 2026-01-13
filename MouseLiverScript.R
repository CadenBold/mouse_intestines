install.packages("BiocManager")

BiocManager::install(c(
  "clusterProfiler",
  "org.Mm.eg.db",   # mouse
  "org.Hs.eg.db",   # human (if needed)
  "enrichplot",
  "DOSE"
))
