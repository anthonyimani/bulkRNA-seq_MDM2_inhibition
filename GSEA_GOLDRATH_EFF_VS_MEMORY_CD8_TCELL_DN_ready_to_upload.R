library(readr)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(enrichplot)
library(ggplot2)

res_df <- read_csv('7.9.2025.DESeq2_results_full.csv')

colnames(res_df)[1] <- "gene_name"

res_df <- res_df %>%
  filter(!is.na(pvalue), !is.na(log2FoldChange)) %>%
  distinct(gene_name, .keep_all = TRUE) %>%
  mutate(
    neg_log10p = -log10(pvalue),
    neg_log10p = ifelse(is.infinite(neg_log10p), 320, neg_log10p),
    signed_stat = sign(log2FoldChange) * neg_log10p
  )

gene_map <- bitr(
  res_df$gene_name,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

res_df <- res_df %>%
  left_join(gene_map, by = c("gene_name" = "SYMBOL")) %>%
  filter(!is.na(ENTREZID))

gene_list <- res_df$signed_stat
names(gene_list) <- res_df$ENTREZID
gene_list <- gene_list[!is.na(names(gene_list))]              
gene_list <- gene_list[!duplicated(names(gene_list))]         
gene_list <- sort(gene_list, decreasing = TRUE)

c7_sets <- msigdbr(species = "Homo sapiens", category = "C7") %>%
  dplyr::select(gs_name, entrez_gene) %>%
  filter(!is.na(entrez_gene))

hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, entrez_gene) %>%
  filter(!is.na(entrez_gene))

em_c7 <- GSEA(
  geneList     = gene_list,
  TERM2GENE    = c7_sets,
  pvalueCutoff = 1
)

hallmark <- GSEA(
  geneList     = gene_list,
  TERM2GENE    = hallmark_sets,
  pvalueCutoff = 1
)

# GOLDRATH_NAIVE_VS_MEMORY_CD8_TCELL_DN - positive
# GOLDRATH_NAIVE_VS_EFF_CD8_TCELL_DN - negative
# GOLDRATH_EFF_VS_MEMORY_CD8_TCELL_UP - negative
# GOLDRATH_EFF_VS_MEMORY_CD8_TCELL_DN - positive, not sig

gseaplot2(
  em_c7,
  geneSetID = "GOLDRATH_EFF_VS_MEMORY_CD8_TCELL_DN",
  title = "",
  base_size = 12,
  ES_geom = "line",
  color = "green",
  pvalue_table = TRUE
) 
