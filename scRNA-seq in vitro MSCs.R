#scRNA seq analysis of cultured MSCs

#load all required libraries
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(sctransform)
library(rgl)
library(cowplot)
library(ggplot2)
library(ape)
library(devtools)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)


packageVersion("Seurat")

#Load data and pre-process
setwd("C:/Users/Patrick Niekamp/Desktop/temp") #set working directory for metadata and H5 files
Sys.setenv('R_MAX_VSIZE'=67108864000)
64000*1024^2
options(future.globals.maxSize= 67108864000)
setwd("../")
sample_summary<-read.csv("ossium_batch2_metadata _12_samples.csv",header = TRUE)

for (i in 1:length(sample_summary$file)){
  cmd_1<-paste(sample_summary$label[i]," <- Read10X_h5('",sample_summary$file[i],"')",sep="")
  eval(parse(text = cmd_1)) 
  cmd_2<-paste(sample_summary$label[i]," <- CreateSeuratObject(",sample_summary$label[i],"$`Gene Expression`,project = '",sample_summary$label[i],"')",sep="")
  eval(parse(text = cmd_2)) 
}

# add metadata and change orig.ident
for (i in 1:dim(sample_summary)[1]){
  cmd_1<-paste(sample_summary$label[i],"@meta.data$donor <- '",sample_summary$donor[i],"'",sep="")
  eval(parse(text = cmd_1)) 
  cmd_2<-paste(sample_summary$label[i],"@meta.data$time <- '",sample_summary$time[i],"'",sep="")
  eval(parse(text = cmd_2)) 
  cmd_3<-paste(sample_summary$label[i],"@meta.data$batch <- 'second'",sep="")
  eval(parse(text = cmd_3)) 
  cmd_4<-paste(sample_summary$label[i],"@meta.data$annotation <- '",sample_summary$pheno[i],"'",sep="")
  eval(parse(text = cmd_4)) 
}

#create mito content
for (i in 1:length(sample_summary$label)){
  cmd_1<-paste(sample_summary$label[i],'[["percent.mt"]] <- PercentageFeatureSet(',sample_summary$label[i],', pattern = "^MT-")',sep="")
  eval(parse(text = cmd_1)) 
}


#filter seurat objects on mito<20% and nGene>400 nUMI < 50000 and nUMI>2000  and nGene < 7500
for (i in 1:length(sample_summary$label)){
  cmd_1<-paste(sample_summary$label[i]," <- subset(",sample_summary$label[i],", subset = nFeature_RNA > 400 &  nFeature_RNA < 8000 & percent.mt < 20 & nCount_RNA < 50000 & nCount_RNA > 200)",sep="")
  eval(parse(text = cmd_1)) 
}

#merge samples
cmd_1<-paste("second.batch.combined <- merge(",sample_summary$label[1],", y = c(",paste(as.character(sample_summary$label[2:length(sample_summary$label)]), collapse=", "),"), project = 'Pluristem.combined.set')",sep="")
eval(parse(text = cmd_1))
Idents(object = second.batch.combined) <- "orig.ident"
#VlnPlot(second.batch.combined, features = c("nFeature_RNA", "nCount_RNA","percent.mt"),group.by = "orig.ident", ncol = 3,pt.size = 0.01)
#VlnPlot(second.batch.combined, features = "nFeature_RNA",pt.size = 0.0) + ggtitle("Genes Detected Per Cell")+NoLegend()
vln_plot <- VlnPlot(second.batch.combined, features = "nFeature_RNA", pt.size = 0.01) +
  ggtitle("Genes Detected Per Cell")
vln_plot + theme(axis.text.x = element_text(angle = 90, hjust = 1))
vln_plot
VlnPlot(second.batch.combined, features = "nCount_RNA",pt.size = 0.0) + ggtitle("UMIs Per Cell")+NoLegend()
VlnPlot(object = second.batch.combined, features = "percent.mt",pt.size = 0.01) + ggtitle("Mitochondrial Content Per Cell")


plot1 <- FeatureScatter(second.batch.combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(second.batch.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

second.batch.combined<-SCTransform(second.batch.combined)
second.batch.combined<-RunPCA(second.batch.combined)
second.batch.combined<-RunUMAP(second.batch.combined,dims = 1:40)
DimPlot(second.batch.combined)
second.batch.combined<-FindNeighbors(second.batch.combined,dims = 1:40)
second.batch.combined<-FindClusters(second.batch.combined,resolution = 0.3)
second.batch.combined <- BuildClusterTree(second.batch.combined,reorder.numeric = T,reorder = T)

DimPlot(second.batch.combined, label = TRUE)
save(second.batch.combined, file = "in_vitro_MSCs.Rdata")

#############
#Further processing of the data 

setwd("") #set wd to folder with saved Rdata file 
load("in_vitro_MSCs.Rdata") 
vitro <- second.batch.combined


p1 <- DimPlot(vitro, group.by = "orig.ident") + NoLegend()  #check distribution of cells based on donor ID
p2 <- DimPlot(vitro, group.by = "seurat_clusters") + NoLegend() #check clusters
p3 <- DimPlot(vitro, group.by = "time") + NoLegend()
p4 <- DimPlot(vitro, group.by = "annotation") + NoLegend()
cowplot::plot_grid(p1, p2, p3, p4)


# select only P4
samples_to_subset <- c("Donor_22063_P4_ifng", "Donor_21320_P4_ifng", "Donor_20043_P4_ifng", "Donor_21069_P4_ifng", 
                       "Donor_22063_P4", "Donor_21320_P4", "Donor_20043_P4", "Donor_21069_P4")
cells_to_keep <- which(Idents(vitro) %in% samples_to_subset)
P4 <- subset(vitro, cells = cells_to_keep)
DimPlot(P4)
P4 <- SCTransform(P4)
P4 <- RunPCA(P4)
P4 <- RunUMAP(P4,dims = 1:20)
P4 <- FindNeighbors(P4, dims = 1:20)
P4 <- FindClusters(P4, resolution = 0.2)
P4 <- BuildClusterTree(P4,reorder.numeric = T,reorder = T)

DimPlot(P4, label = TRUE)

#Remove cluster 4, outlier
Idents(object = P4) <- "seurat_clusters"
samples_to_subset <- c("0", "1", "2", "3")
cells_to_keep <- which(Idents(P4) %in% samples_to_subset)
P4 <- subset(P4, cells = cells_to_keep)

saveRDS(P4, file ="P4.rds")
P4 <- readRDS("P4.rds")

p1 <- DimPlot(P4, group.by = "orig.ident") + NoLegend()  
p2 <- DimPlot(P4, group.by = "seurat_clusters") + NoLegend() 
p3 <- DimPlot(P4, group.by = "annotation") + NoLegend()

#Figure 2B
cowplot::plot_grid(p1, p2, p3, nrow = 3)

# Number per cells in each cluster
# Create a table of sample IDs (orig.ident) by cluster IDs (seurat_clusters)
cell_counts <- table(P4$orig.ident, P4$seurat_clusters)

# Convert to a data frame for easier viewing
cell_counts_df <- as.data.frame(cell_counts)
colnames(cell_counts_df) <- c("Sample", "Cluster", "Cell_Count")
write.csv(cell_counts_df, file = "cell_counts.csv")

# Figure S2A
ggplot(cell_counts_df, aes(x = Cluster, y = Cell_Count, fill = Sample)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.75), width = 0.7) +
  scale_x_discrete(expand = expansion(mult = c(0.1, 0.1))) + # Adds padding on x-axis
  labs(title = "Cell Counts per Cluster by Sample",
       x = "Seurat Cluster",
       y = "Cell Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Calculate DEGs for heatmap
Markers <- FindAllMarkers(P4, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

# Fig S2B
DoHeatmap(P4, features = top10$gene, size = 6)

# Fig S3A
p1 <- FeaturePlot(P4,"IFITM1", ncol = 1, pt.size = 3 + theme(legend.position = 'none', axis.text.x = element_text(size = 19)))
p2 <- FeaturePlot(P4,"IFI27", ncol = 1, pt.size = 3 + theme(legend.position = 'none', axis.text.x = element_text(size = 19)))
p3 <- FeaturePlot(P4,"HLA-DRA", ncol = 1, pt.size = 3 + theme(legend.position = 'none', axis.text.x = element_text(size = 19)))
p4 <- FeaturePlot(P4,"HLA-DMA", ncol = 1, pt.size = 3 + theme(legend.position = 'none', axis.text.x = element_text(size = 19)))
cowplot::plot_grid(p1, p2, p3, p4)


#Pathway analysis

# Fig S2C -> increased in C0
top_100_genes <- Markers %>%
  filter(cluster == 0) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 100)
head(top_100_genes, n = 15)
top_100_genes<- top_100_genes$gene
GO_results <- enrichGO(gene = top_100_genes, OrgDb = "org.Hs.eg.db", key = "SYMBOL", ont = "BP")
as.data.frame(GO_results)
fit <- plot(barplot(GO_results, showCategory = 10))

# Fig S2C -> increased in C1
top_100_genes <- Markers %>%
  filter(cluster == 1) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 100)
head(top_100_genes, n = 15)
top_100_genes<- top_100_genes$gene
GO_results <- enrichGO(gene = top_100_genes, OrgDb = "org.Hs.eg.db", key = "SYMBOL", ont = "BP")
as.data.frame(GO_results)
fit <- plot(barplot(GO_results, showCategory = 10))

# Fig S2C -> increased in C2
top_100_genes <- Markers %>%
  filter(cluster == 1) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 100)
head(top_100_genes, n = 15)
top_100_genes<- top_100_genes$gene
GO_results <- enrichGO(gene = top_100_genes, OrgDb = "org.Hs.eg.db", key = "SYMBOL", ont = "BP")
as.data.frame(GO_results)
fit <- plot(barplot(GO_results, showCategory = 10))

# Fig S2C -> increased in C3
top_100_genes <- Markers %>%
  filter(cluster == 1) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 100)
head(top_100_genes, n = 15)
top_100_genes<- top_100_genes$gene
GO_results <- enrichGO(gene = top_100_genes, OrgDb = "org.Hs.eg.db", key = "SYMBOL", ont = "BP")
as.data.frame(GO_results)
fit <- plot(barplot(GO_results, showCategory = 10))


# select only samples treated with IFNg
samples_to_subset <- c("Donor_22063_P4_ifng", "Donor_21320_P4_ifng", "Donor_20043_P4_ifng", "Donor_21069_P4_ifng")
cells_to_keep <- which(Idents(vitro) %in% samples_to_subset)
IFNG <- subset(vitro, cells = cells_to_keep)
DimPlot(IFNG)
IFNG <- SCTransform(IFNG)
IFNG<-RunPCA(IFNG)
IFNG<-RunUMAP(IFNG,dims = 1:40)
IFNG <- FindNeighbors(IFNG, dims = 1:40)
IFNG <- FindClusters(IFNG, resolution = 0.4)
IFNG<-BuildClusterTree(IFNG,reorder.numeric = T,reorder = T)

#Fig S3B
p1 <- DimPlot(IFNG, group.by = "seurat_clusters", pt.size = 1, label = TRUE, label.size = 8) + NoLegend()
p2 <- DimPlot(IFNG, group.by = "orig.ident", pt.size = 1) + NoLegend()
cowplot::plot_grid(p1, p2)

saveRDS(IFNG, file ="IFNG.rds")
IFNG <- readRDS("IFNG.rds")

# Figure 2C -> 069 vs 043
unique(IFNG$orig.ident)
Idents(IFNG) <- "orig.ident"
Markers_IFNg <- FindMarkers(IFNG, ident.1 = "Donor_21069_P4_ifng", ident.2 = "Donor_20043_P4_ifng")
write.csv(Markers_IFNg, file = "Marker_069_vs_043.csv")
Markers_IFNg$Significance <- "Not Significant"
Markers_IFNg$Significance[Markers_IFNg$p_val_adj < 0.01 & abs(Markers_IFNg$avg_log2FC) > 1] <- "Significant"


# Select top markers for each cluster comparison (positive log fold change)
top_pos_markers <- Markers_IFNg %>%
  top_n(n = 30, wt = avg_log2FC)
top_neg_markers <- Markers_IFNg %>%
  top_n(n = 30, wt = -avg_log2FC)

top_markers <- rbind(top_pos_markers, top_neg_markers)

# Create a volcano plot
ggplot(Markers_IFNg, aes(x = avg_log2FC, y = -log10(p_val_adj), color = Significance)) +
  geom_point(alpha = 0.8, size = 3) +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "gray")) +
  theme_minimal() +
  labs(title = "Volcano Plot of Differentially Expressed Genes",
       x = "log2 Fold Change",
       y = "-log10 Adjusted P-Value") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", size = 0.5) +
  theme(legend.position = "right") +
  geom_text(data = top_markers, aes(label = rownames(top_markers)), vjust = -0.5, hjust = 0.5, size = 3, color = "black")


# Figure 2C -> 069 vs 320
unique(IFNG$orig.ident)
Idents(IFNG) <- "orig.ident"
Markers_IFNg <- FindMarkers(IFNG, ident.1 = "Donor_21069_P4_ifng", ident.2 = "Donor_21320_P4_ifng")
write.csv(Markers_IFNg, file = "Marker_069_vs_320.csv")
Markers_IFNg$Significance <- "Not Significant"
Markers_IFNg$Significance[Markers_IFNg$p_val_adj < 0.01 & abs(Markers_IFNg$avg_log2FC) > 1] <- "Significant"


# Select top markers for each cluster comparison (positive log fold change)
top_pos_markers <- Markers_IFNg %>%
  top_n(n = 30, wt = avg_log2FC)
top_neg_markers <- Markers_IFNg %>%
  top_n(n = 30, wt = -avg_log2FC)

top_markers <- rbind(top_pos_markers, top_neg_markers)

# Create a volcano plot
ggplot(Markers_IFNg, aes(x = avg_log2FC, y = -log10(p_val_adj), color = Significance)) +
  geom_point(alpha = 0.8, size = 3) +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "gray")) +
  theme_minimal() +
  labs(title = "Volcano Plot of Differentially Expressed Genes",
       x = "log2 Fold Change",
       y = "-log10 Adjusted P-Value") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", size = 0.5) +
  theme(legend.position = "right") +
  geom_text(data = top_markers, aes(label = rownames(top_markers)), vjust = -0.5, hjust = 0.5, size = 3, color = "black")



# Figure 2D -> 063 vs 043
unique(IFNG$orig.ident)
Idents(IFNG) <- "orig.ident"
Markers_IFNg <- FindMarkers(IFNG, ident.1 = "Donor_22063_P4_ifng", ident.2 = "Donor_20043_P4_ifng")
write.csv(Markers_IFNg, file = "Marker_063_vs_043.csv")
Markers_IFNg$Significance <- "Not Significant"
Markers_IFNg$Significance[Markers_IFNg$p_val_adj < 0.01 & abs(Markers_IFNg$avg_log2FC) > 1] <- "Significant"


# Select top markers for each cluster comparison (positive log fold change)
top_pos_markers <- Markers_IFNg %>%
  top_n(n = 30, wt = avg_log2FC)
top_neg_markers <- Markers_IFNg %>%
  top_n(n = 30, wt = -avg_log2FC)

top_markers <- rbind(top_pos_markers, top_neg_markers)

# Create a volcano plot
ggplot(Markers_IFNg, aes(x = avg_log2FC, y = -log10(p_val_adj), color = Significance)) +
  geom_point(alpha = 0.8, size = 3) +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "gray")) +
  theme_minimal() +
  labs(title = "Volcano Plot of Differentially Expressed Genes",
       x = "log2 Fold Change",
       y = "-log10 Adjusted P-Value") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", size = 0.5) +
  theme(legend.position = "right") +
  geom_text(data = top_markers, aes(label = rownames(top_markers)), vjust = -0.5, hjust = 0.5, size = 3, color = "black")



# Figure 2D -> 063 vs 320
unique(IFNG$orig.ident)
Idents(IFNG) <- "orig.ident"
Markers_IFNg <- FindMarkers(IFNG, ident.1 = "Donor_22063_P4_ifng", ident.2 = "Donor_21320_P4_ifng")
write.csv(Markers_IFNg, file = "Marker_063_vs_320.csv")
Markers_IFNg$Significance <- "Not Significant"
Markers_IFNg$Significance[Markers_IFNg$p_val_adj < 0.01 & abs(Markers_IFNg$avg_log2FC) > 1] <- "Significant"


# Select top markers for each cluster comparison (positive log fold change)
top_pos_markers <- Markers_IFNg %>%
  top_n(n = 30, wt = avg_log2FC)
top_neg_markers <- Markers_IFNg %>%
  top_n(n = 30, wt = -avg_log2FC)

top_markers <- rbind(top_pos_markers, top_neg_markers)

# Create a volcano plot
ggplot(Markers_IFNg, aes(x = avg_log2FC, y = -log10(p_val_adj), color = Significance)) +
  geom_point(alpha = 0.8, size = 3) +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "gray")) +
  theme_minimal() +
  labs(title = "Volcano Plot of Differentially Expressed Genes",
       x = "log2 Fold Change",
       y = "-log10 Adjusted P-Value") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", size = 0.5) +
  theme(legend.position = "right") +
  geom_text(data = top_markers, aes(label = rownames(top_markers)), vjust = -0.5, hjust = 0.5, size = 3, color = "black")



# Figure S3C -> 069 vs 063
unique(IFNG$orig.ident)
Idents(IFNG) <- "orig.ident"
Markers_IFNg <- FindMarkers(IFNG, ident.1 = "Donor_21069_P4_ifng", ident.2 = "Donor_22063_P4_ifng")
write.csv(Markers_IFNg, file = "Marker_069_vs_063.csv")
Markers_IFNg$Significance <- "Not Significant"
Markers_IFNg$Significance[Markers_IFNg$p_val_adj < 0.01 & abs(Markers_IFNg$avg_log2FC) > 1] <- "Significant"


# Select top markers for each cluster comparison (positive log fold change)
top_pos_markers <- Markers_IFNg %>%
  top_n(n = 30, wt = avg_log2FC)
top_neg_markers <- Markers_IFNg %>%
  top_n(n = 30, wt = -avg_log2FC)

top_markers <- rbind(top_pos_markers, top_neg_markers)

# Create a volcano plot
ggplot(Markers_IFNg, aes(x = avg_log2FC, y = -log10(p_val_adj), color = Significance)) +
  geom_point(alpha = 0.8, size = 3) +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "gray")) +
  theme_minimal() +
  labs(title = "Volcano Plot of Differentially Expressed Genes",
       x = "log2 Fold Change",
       y = "-log10 Adjusted P-Value") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", size = 0.5) +
  theme(legend.position = "right") +
  geom_text(data = top_markers, aes(label = rownames(top_markers)), vjust = -0.5, hjust = 0.5, size = 3, color = "black")



# Figure S3C -> 320 vs 043
unique(IFNG$orig.ident)
Idents(IFNG) <- "orig.ident"
Markers_IFNg <- FindMarkers(IFNG, ident.1 = "Donor_21320_P4_ifng", ident.2 = "Donor_20043_P4_ifng")
write.csv(Markers_IFNg, file = "Marker_320_vs_043.csv")
Markers_IFNg$Significance <- "Not Significant"
Markers_IFNg$Significance[Markers_IFNg$p_val_adj < 0.01 & abs(Markers_IFNg$avg_log2FC) > 1] <- "Significant"


# Select top markers for each cluster comparison (positive log fold change)
top_pos_markers <- Markers_IFNg %>%
  top_n(n = 30, wt = avg_log2FC)
top_neg_markers <- Markers_IFNg %>%
  top_n(n = 30, wt = -avg_log2FC)

top_markers <- rbind(top_pos_markers, top_neg_markers)

# Create a volcano plot
ggplot(Markers_IFNg, aes(x = avg_log2FC, y = -log10(p_val_adj), color = Significance)) +
  geom_point(alpha = 0.8, size = 3) +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "gray")) +
  theme_minimal() +
  labs(title = "Volcano Plot of Differentially Expressed Genes",
       x = "log2 Fold Change",
       y = "-log10 Adjusted P-Value") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", size = 0.5) +
  theme(legend.position = "right") +
  geom_text(data = top_markers, aes(label = rownames(top_markers)), vjust = -0.5, hjust = 0.5, size = 3, color = "black")



#Fig S4A
p1 <- DimPlot(P4, label = TRUE) + NoLegend()
p2 <- DimPlot(P4, group.by = "orig.ident") + NoLegend()
p3 <- FeaturePlot(P4,"IDO1", ncol = 1, pt.size = 3 + theme(legend.position = 'none', axis.text.x = element_text(size = 19)))
cowplot::plot_grid(p1, p2, p3, ncol = 3)




