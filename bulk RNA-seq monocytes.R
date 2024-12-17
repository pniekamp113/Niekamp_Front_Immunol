#bulk RNA sequencing analysis of monocytes treated with MSC-CM for 24hrs

#load required libraries 


library(httr)
library(clusterProfiler)
library(EnhancedVolcano)
library(DESeq2)
library(tidyverse)
library(dplyr)
library(org.Hs.eg.db)
library(GSVA)
library(VennDiagram)
library(ggVennDiagram)
library(RColorBrewer)
library(gplots)
library(ggplot2)
library(reshape2)
library(pheatmap)

setwd("") #set wd to read in raw FPKM data file

# Read your CSV file
rna_seq_data <- read.csv("gene_fpkm_edit for R analysis.csv", header = TRUE)

# Define a threshold for filtering
threshold <- 2

# Filter out rows where more than 13 columns have FPKM values below the threshold (at least 3 columns have to have >2FPKM)
filtered_data <- rna_seq_data[rowSums(rna_seq_data[, -1] < threshold) <= 13, ]
filtered_data <- filtered_data[rowSums(filtered_data[, -1] == 0) < 3, ] 

# Save filtered_data to a CSV file
write.csv(filtered_data, "filtered_data.csv", row.names = FALSE)


# Multiply FPKM values by a scaling factor (e.g., 10,000) to convert to integers. Necessary as dds needs integers
scaling_factor <- 1e4
filtered_data[, -1] <- round(filtered_data[, -1] * scaling_factor)

write.csv(filtered_data, file = "filtered_data")
csv_file_path <- "C:\\Users\\Patrick Niekamp\\Desktop\\Experiments\\Exp. 18_RNA seq\\R files\\filtered_data.csv"

# Create a sample information dataframe
sample_info <- data.frame(
  sample_id = colnames(filtered_data)[-1],  # assuming the first column is gene names
  group = c(rep("Ctrl", 5), rep("MSC-CM", 5), rep("M1", 3), rep("M2", 3))
)

# Convert 'group' column to a factor with numeric levels
sample_info$group <- factor(sample_info$group, levels = c("Ctrl", "MSC-CM", "M1", "M2"))

# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = as.matrix(filtered_data[, -1]), colData = sample_info, design = ~ group)

# Run DESeq2 analysis
dds <- DESeq(dds)

# Get differential expression results for each comparison
results_ctrl_msc <- results(dds, contrast = c("group", "Ctrl", "MSC-CM"))
results_ctrl_m1 <- results(dds, contrast = c("group", "Ctrl", "M1"))
results_ctrl_m2 <- results(dds, contrast = c("group", "Ctrl", "M2"))
results_msc_m1 <- results(dds, contrast = c("group", "MSC-CM", "M1"))
results_msc_m2 <- results(dds, contrast = c("group", "MSC-CM", "M2"))
results_m1_m2 <- results(dds, contrast = c("group", "M1", "M2"))
results_msc_ctrl <- results(dds, contrast = c("group", "MSC-CM", "Ctrl"))
results_m1_ctrl <- results(dds, contrast = c("group", "M1", "Ctrl"))
results_M2_Ctrl <- results(dds, contrast = c("group", "M2", "Ctrl"))


# Extract gene names from the first column
gene_names <- filtered_data$gene_name

# Assign gene names to results
results_ctrl_msc$gene_name <- gene_names
results_ctrl_m1$gene_name <- gene_names
results_ctrl_m2$gene_name <- gene_names
results_msc_m1$gene_name <- gene_names
results_msc_m2$gene_name <- gene_names
results_m1_m2$gene_name <- gene_names
results_msc_ctrl$gene_name <- gene_names
results_m1_ctrl$gene_name <- gene_names
results_M2_Ctrl$gene_name <- gene_names

# Save DESeq2 results to files
saveRDS(results_ctrl_msc, "C:\\Users\\Patrick Niekamp\\Desktop\\Experiments\\Exp. 18_RNA seq\\R files\\results_ctrl_msc.rds")
saveRDS(results_ctrl_m1, "C:\\Users\\Patrick Niekamp\\Desktop\\Experiments\\Exp. 18_RNA seq\\R files\\results_ctrl_m1.rds")
saveRDS(results_ctrl_m2, "C:\\Users\\Patrick Niekamp\\Desktop\\Experiments\\Exp. 18_RNA seq\\R files\\results_ctrl_m2.rds")
saveRDS(results_msc_m1, "C:\\Users\\Patrick Niekamp\\Desktop\\Experiments\\Exp. 18_RNA seq\\R files\\results_msc_m1.rds")
saveRDS(results_msc_m2, "C:\\Users\\Patrick Niekamp\\Desktop\\Experiments\\Exp. 18_RNA seq\\R files\\results_msc_m2.rds")
saveRDS(results_m1_m2, "C:\\Users\\Patrick Niekamp\\Desktop\\Experiments\\Exp. 18_RNA seq\\R files\\results_m1_m2.rds")
saveRDS(results_msc_ctrl, "C:\\Users\\Patrick Niekamp\\Desktop\\Experiments\\Exp. 18_RNA seq\\R files\\results_msc_ctrl.rds")
saveRDS(results_m1_ctrl, "C:\\Users\\Patrick Niekamp\\Desktop\\Experiments\\Exp. 18_RNA seq\\R files\\results_m1_ctrl.rds")
saveRDS(results_M2_Ctrl, "C:\\Users\\Patrick Niekamp\\Desktop\\Experiments\\Exp. 18_RNA seq\\R files\\results_m2_ctrl.rds")


# Gene Ontology analysis

results_msc_ctrl <- readRDS("C:\\Users\\Patrick Niekamp\\Desktop\\Experiments\\Exp. 18_RNA seq\\R files\\results_msc_ctrl.rds")

# Extract differentially expressed genes with log2 fold change above the threshold
de_genes_up <- results_msc_ctrl$gene_name[results_msc_ctrl$log2FoldChange > 1]

# Perform gene ontology enrichment analysis for up-regulated genes
enrich_result_up <- enrichGO(
  gene = de_genes_up,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pvalueCutoff = 0.01,
  qvalueCutoff = 0.2
)

write.csv(enrich_result_up, file = "pathways_enriched_in_msc_vs_ctrl.csv")

# Fig 4A
fit <- plot(barplot(enrich_result_up, showCategory = 12))


# extract genes for specific pathways (e.g. leukocyte migration)
head(as.data.frame(enrich_result_up), n = 20)

#extract genes for leukocyte migration
df <- as.data.frame(enrich_result_up)
gene_list <- df[df$ID == "GO:0050900", "geneID"]
genes <- unlist(strsplit(gene_list, "/"))

# load filtered data
filtered_data <- read.csv("C:\\Users\\Patrick Niekamp\\Desktop\\Experiments\\Exp. 18_RNA seq\\R files\\filtered_data.csv")
head(filtered_data)

# Extract only rows corresponding to the genes of interest
filtered_genes_data <- filtered_data %>%
  filter(gene_name %in% genes)

gene_names <- filtered_genes_data$gene_name
numeric_data <- dplyr::select(filtered_genes_data, -gene_name)

average_OH1_OH5 <- rowMeans(numeric_data[, c("OH01", "OH02", "OH03", "OH04", "OH05")])
log2FC_data <- log2(numeric_data / average_OH1_OH5)

# Add back gene names for reference
log2FC_data <- cbind(gene_name = gene_names, log2FC_data)

row_names <- log2FC_data$gene_name
heatmap_data <- log2FC_data %>% 
  dplyr::select(-gene_name) %>% 
  as.matrix()

rownames(heatmap_data) <- row_names

heatmap_data[is.na(heatmap_data)] <- 0
heatmap_data[is.nan(heatmap_data)] <- 0
heatmap_data[is.infinite(heatmap_data)] <- 0


# Convert matrix to data frame for ggplot
df_heatmap <- as.data.frame(heatmap_data)
df_heatmap$gene_name <- rownames(heatmap_data) # Add gene names as a new column

write.csv(df_heatmap, file = "leukocyte migration_heatmap.csv")

# Reshape data from wide to long format using reshape2::melt()
df_long <- melt(df_heatmap, id.vars = "gene_name", variable.name = "Sample", value.name = "Log2FC")

# Fig 4B
ggplot(df_long, aes(x = Sample, y = gene_name, fill = Log2FC)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = "Sample", y = "Gene", title = "Heatmap of Log2FC Values")


# Pathway "phagocytosis
#extract genes 
df <- as.data.frame(enrich_result_up)
gene_list <- df[df$ID == "GO:0006909", "geneID"]
genes <- unlist(strsplit(gene_list, "/"))

# load filtered data
filtered_data <- read.csv("C:\\Users\\Patrick Niekamp\\Desktop\\Experiments\\Exp. 18_RNA seq\\R files\\filtered_data.csv")
head(filtered_data)

# Extract only rows corresponding to the genes of interest
filtered_genes_data <- filtered_data %>%
  filter(gene_name %in% genes)

gene_names <- filtered_genes_data$gene_name
numeric_data <- dplyr::select(filtered_genes_data, -gene_name)

average_OH1_OH5 <- rowMeans(numeric_data[, c("OH01", "OH02", "OH03", "OH04", "OH05")])
log2FC_data <- log2(numeric_data / average_OH1_OH5)

# Add back gene names for reference
log2FC_data <- cbind(gene_name = gene_names, log2FC_data)

row_names <- log2FC_data$gene_name
heatmap_data <- log2FC_data %>% 
  dplyr::select(-gene_name) %>% 
  as.matrix()

rownames(heatmap_data) <- row_names

heatmap_data[is.na(heatmap_data)] <- 0
heatmap_data[is.nan(heatmap_data)] <- 0
heatmap_data[is.infinite(heatmap_data)] <- 0


# Convert matrix to data frame for ggplot
df_heatmap <- as.data.frame(heatmap_data)
df_heatmap$gene_name <- rownames(heatmap_data) # Add gene names as a new column

write.csv(df_heatmap, file = "phagocytosis_heatmap.csv")

# Reshape data from wide to long format using reshape2::melt()
df_long <- melt(df_heatmap, id.vars = "gene_name", variable.name = "Sample", value.name = "Log2FC")

# Fig 5A
ggplot(df_long, aes(x = Sample, y = gene_name, fill = Log2FC)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = "Sample", y = "Gene", title = "Heatmap of Log2FC Values")

