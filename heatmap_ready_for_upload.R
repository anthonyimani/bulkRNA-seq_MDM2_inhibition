library(readr)
library(dplyr)
library(pheatmap)

clean_data <- data %>%
  filter(gene_name != "---") %>%                  
  distinct(gene_name, .keep_all = TRUE)           

expression_data <- clean_data %>%
  dplyr::select(`2-2`, `1-1`, `2-3`, `1-2`, `1-3`, `2-1`) %>%
  as.data.frame()
rownames(expression_data) <- clean_data$gene_name

top_var_genes <- apply(expression_data, 1, var)
top_genes <- order(top_var_genes, decreasing = TRUE)[1:2000]  
expression_top <- expression_data[top_genes, ]

pheatmap(expression_top,
         scale = "row",                               
         clustering_distance_rows = "euclidean",      
         clustering_distance_cols = "euclidean",      
         clustering_method = "ward.D2",               # Ward's method for clustering
         color = colorRampPalette(c("dodgerblue1", "white", "red"))(100),
         fontsize_row = 6,
         show_rownames = FALSE,                       
         main = "Unsupervised Hierarchical Clustering of Gene Expression")
