library(DESeq2)
library(readr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

raw_counts <- read_tsv("rounded_raw_counts.tsv")
raw_counts <- as.data.frame(raw_counts)  

count_data <- raw_counts[, c("gene_name", "2-2", "1-1", "2-3", "1-2", "1-3", "2-1")]

count_data <- count_data[count_data$gene_name != "---", ]

count_data <- count_data[!duplicated(count_data$gene_name), ]

rownames(count_data) <- count_data$gene_name
count_data$gene_name <- NULL


sample_names <- colnames(count_data)
condition <- ifelse(grepl("^1-", sample_names), "control",
                    ifelse(grepl("^2-", sample_names), "treated", NA))
coldata <- data.frame(row.names = sample_names,
                      condition = factor(condition, levels = c("control", "treated")))

print(coldata)

dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = coldata,
                              design = ~ condition)

print(dds)

dds <- dds[rowSums(counts(dds)) > 10, ]


dds <- DESeq(dds)

summary(results(dds))

normalized_counts <- counts(dds, normalized = TRUE)
write.csv(as.data.frame(normalized_counts), "7.9.2025.DESeq2_normalized_counts.csv")

res <- results(dds)
res <- res[order(res$padj), ]
write.csv(as.data.frame(res), "7.9.2025.DESeq2_results_full.csv")

library(EnhancedVolcano)

library(EnhancedVolcano)

res$pvalue[res$pvalue == 0] <- 1e-320
res$padj[res$padj == 0] <- 1e-320

labeled_genes <- c('CCR7','TOX2','NR4A1','TP53I3','CDK1','CDK2',
                   'CDKN1A','MKI67','CCNB1','CCNA2','SOX4','IFNG','GZMA','PRF1','EOMES','CCL1','CSPG4','CD160','CD244','AREG','CYP4F3','SULF2','BUB1B','DLGAP5',
                   'KIF15','FASLG')

keyvals <- ifelse(
  rownames(res) %in% labeled_genes, 'purple',
  ifelse(res$log2FoldChange < -1 & res$padj < 0.01, 'dodgerblue',
         ifelse(res$log2FoldChange > 1 & res$padj < 0.01, 'red',
                'black')))
keyvals[is.na(keyvals)] <- 'black'

names(keyvals)[keyvals == 'purple'] <- 'Genes of Interest'
names(keyvals)[keyvals == 'red'] <- 'Upregulated'
names(keyvals)[keyvals == 'dodgerblue'] <- 'Downregulated'
names(keyvals)[keyvals == 'black'] <- 'Not significant'

# EnhancedVolcano plot
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                xlim = c(-12, 12),
                selectLab = labeled_genes,
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 0.01,
                FCcutoff = 1.0,
                pointSize = 3.5,          
                labSize = 6.0,            
                labCol = 'purple',
                labFace = 'bold',
                boxedLabels = TRUE,       
                colCustom = keyvals,      
                colAlpha = 0.4,           
                legendPosition = 'right',
                drawConnectors = TRUE,    
                widthConnectors = 0.25,    
                colConnectors = 'black'
)



