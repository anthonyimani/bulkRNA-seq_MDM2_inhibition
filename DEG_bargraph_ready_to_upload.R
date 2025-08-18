library(dplyr)
library(ggplot2)
library(scales)

res_df <- as.data.frame(res)
res_df$gene_name <- rownames(res_df)

genes_positive <- c(
  'CCR7','TOX2','NR4A1','TP53I3','CDK1','CDK2',
  'CDKN1A','MKI67','CCNB1','CCNA2','SOX4','IFNG','GZMA','PRF1','EOMES','CCL1','CSPG4','CD160','CD244','AREG','CYP4F3','SULF2','BUB1B','DLGAP5',
  'KIF15','FASLG'
)

genes_negative <- c(
  'CCR7','TOX2','NR4A1','TP53I3','CDK1','CDK2',
  'CDKN1A','MKI67','CCNB1','CCNA2','SOX4','IFNG','GZMA','PRF1','EOMES','CCL1','CSPG4','CD160','CD244','AREG','CYP4F3','SULF2','BUB1B','DLGAP5',
  'KIF15','FASLG'
)

combined_df <- res_df %>%
  filter(!is.na(log2FoldChange)) %>%
  filter(
    (gene_name %in% genes_positive & log2FoldChange > 0) |
      (gene_name %in% genes_negative & log2FoldChange < 0)
  ) %>%
  dplyr::select(gene_name, log2FoldChange) %>%
  arrange(desc(log2FoldChange)) %>%
  mutate(
    gene_name = factor(gene_name, levels = gene_name),
    group = ifelse(log2FoldChange > 0, "Up", "Down")
  )

ggplot(combined_df, aes(x = gene_name, y = log2FoldChange, fill = log2FoldChange)) +
  geom_col() +
  scale_fill_gradientn(
    colors = c("dodgerblue", "skyblue", "white", "pink", "red"),
    values = rescale(c(min(combined_df$log2FoldChange),
                       -2, 0, 2,
                       max(combined_df$log2FoldChange))),
    name = "log2FC"
  ) +
  labs(
    x = "",
    y = "log2(Fold Change)",
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 16),
    axis.text.y = element_text(size = 13),
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14)
  )

