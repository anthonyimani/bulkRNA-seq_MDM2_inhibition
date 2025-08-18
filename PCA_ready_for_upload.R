library(readr)
library(dplyr)
library(ggplot2)
library(plotly)

data <- read_csv("7.9.2025.DESeq2_normalized_counts.csv")

expression_data <- data %>%
  select(`2-2`, `1-1`, `2-3`, `1-2`, `1-3`, `2-1`)

expression_matrix <- t(as.matrix(expression_data))
sample_names <- rownames(expression_matrix)

pca_result <- prcomp(expression_matrix, scale. = TRUE)
pca_df <- as.data.frame(pca_result$x)
pca_df$Sample <- sample_names

group_assignments <- data.frame(
  Sample = c("1-1", "1-2", "1-3", "2-1", "2-2", "2-3"),
  Group = c("0 uM NAV", "0 uM NAV", "0 uM NAV", "2 uM NAV", "2 uM NAV", "2 uM NAV")
)

pca_df <- merge(pca_df, group_assignments, by = "Sample")

fig <- plot_ly(
  pca_df, 
  x = ~PC1, y = ~PC2, z = ~PC3, 
  type = "scatter3d", 
  mode = "markers+text",
  color = ~Group,
  colors = c("0 uM NAV" = "purple3", "2 uM NAV" = "darkorange1"),
  text = ~Sample,
  textposition = "top center",
  marker = list(
    size = 10,                       # Marker size
    opacity = 0.8,                  # Transparency
    line = list(color = "black", width = 1)  # Black border
  )
) %>%
  layout(title = "")

fig