# R/deg_functions.R

run_deseq2 <- function(counts, group, comparison) {
  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = data.frame(condition = group),
    design = ~ condition
  )
  
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("condition", comparison, levels(group)[1]))
  
  norm_counts <- counts(dds, normalized = TRUE)
  
  return(list(
    deg_results = as.data.frame(res),
    norm_counts = norm_counts,
    dds = dds
  ))
}

get_top_genes <- function(deg_results, n = 50, pval_cutoff = 0.05, lfc_cutoff = 1) {
  sig_genes <- deg_results[deg_results$padj < pval_cutoff & 
                           abs(deg_results$log2FoldChange) > lfc_cutoff, ]
  
  # Order by p-value and absolute log2FC
  sig_genes <- sig_genes[order(sig_genes$padj, -abs(sig_genes$log2FoldChange)), ]
  
  if (nrow(sig_genes) == 0) {
    return(character(0))
  }
  
  n <- min(n, nrow(sig_genes))
  rownames(sig_genes)[1:n]
}
