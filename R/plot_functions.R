# R/plot_functions.R

plot_volcano <- function(deg_results, pval_cutoff = 0.05, lfc_cutoff = 1, title = "Volcano Plot") {
  deg_results$significant <- ifelse(
    deg_results$padj < pval_cutoff & abs(deg_results$log2FoldChange) > lfc_cutoff,
    ifelse(deg_results$log2FoldChange > lfc_cutoff, "Up", "Down"),
    "Not significant"
  )
  
  deg_results$log10_padj <- -log10(deg_results$padj)
  
  p <- ggplot(deg_results, aes(x = log2FoldChange, y = log10_padj, 
                              color = significant, 
                              text = paste("Gene:", rownames(deg_results),
                                          "<br>log2FC:", round(log2FoldChange, 2),
                                          "<br>p.adj:", format.pval(padj)))) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c("Up" = "red", "Down" = "blue", "Not significant" = "grey")) +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed") +
    geom_hline(yintercept = -log10(pval_cutoff), linetype = "dashed") +
    labs(title = title, x = "log2 Fold Change", y = "-log10(Adjusted p-value)") +
    theme_minimal() +
    theme(legend.title = element_blank())
  
  ggplotly(p, tooltip = "text")
}

plot_pca <- function(norm_counts, metadata, condition_col, title = "PCA Plot") {
  pca <- prcomp(t(norm_counts), scale. = TRUE)
  percent_var <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
  
  pca_data <- data.frame(
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    Condition = metadata[[condition_col]]
  )
  
  p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Condition, 
                           text = paste("Sample:", rownames(pca_data)))) +
    geom_point(size = 3) +
    xlab(paste0("PC1 (", percent_var[1], "%)")) +
    ylab(paste0("PC2 (", percent_var[2], "%)")) +
    ggtitle(title) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  ggplotly(p, tooltip = "text")
}

plot_ma <- function(deg_results, pval_cutoff, lfc_cutoff, title = "MA Plot") {
  # Ensure deg_results is properly formatted
  req(deg_results)
  
  # Safely calculate mean expression
  deg_results$mean_expr <- tryCatch({
    if ("baseMean" %in% colnames(deg_results)) {
      deg_results$baseMean  # DESeq2 output
    } else if ("AveExpr" %in% colnames(deg_results)) {
      deg_results$AveExpr   # limma output
    } else {
      # Fallback calculation if neither column exists
      if (is.null(dim(deg_results$counts))) {
        stop("Count matrix has lost its dimensions")
      }
      rowMeans(deg_results$counts, na.rm = TRUE)
    }
  }, error = function(e) {
    stop("Failed to calculate mean expression: ", e$message)
  })
  
  # Create significance factor
  deg_results$significant <- ifelse(
    deg_results$padj < pval_cutoff & abs(deg_results$log2FoldChange) > lfc_cutoff,
    ifelse(deg_results$log2FoldChange > lfc_cutoff, "Up", "Down"),
    "Not significant"
  )
  
  # Create the plot
  p <- ggplot(deg_results, aes(x = mean_expr, y = log2FoldChange, 
                              color = significant,
                              text = paste("Gene:", rownames(deg_results),
                                          "<br>Mean expr:", round(mean_expr, 2),
                                          "<br>log2FC:", round(log2FoldChange, 2),
                                          "<br>p.adj:", format.pval(padj)))) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c("Up" = "red", "Down" = "blue", "Not significant" = "grey")) +
    geom_hline(yintercept = 0, color = "black") +
    geom_hline(yintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed") +
    labs(title = title, x = "Mean expression", y = "log2 Fold Change") +
    theme_minimal() +
    theme(legend.title = element_blank())
  
  ggplotly(p, tooltip = "text")
}
