# scRNA seq analysis of publicly available data: GSE253355
# citation: Bandyopadhyay, S., Duffy, M. P., Ahn, K. J., Sussman, J. H., Pang, M., Smith, D., Duncan, G., Zhang, I., Huang, J., Lin, Y., Xiong, B., Imtiaz, T., Chen, C.-H., Thadi, A., Chen, C., Xu, J., Reichart, M., Martinez, Z., Diorio, C., . . . Tan, K. Mapping the cellular biogeography of human bone marrow niches using single-cell transcriptomics and proteomic imaging. Cell. 2024. https://doi.org/10.1016/j.cell.2024.04.013

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
library(CellChat)
library(patchwork)


# Download and unzip "GSE253355_Normal_Bone_Marrow_Atlas_Seurat_SB_v2.rds.gz"	from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253355

#load seurat object
setwd("")

BMN <- readRDS("C:\\Users\\Patrick Niekamp\\Desktop\\Experiments\\Exp. 57_Cell paper scRNA seq\\GSE253355_Normal_Bone_Marrow_Atlas_Seurat_SB_v2.rds")

#create cellchat object
cellchat <- createCellChat(object = BMN)
CellChatDB <- CellChatDB.human 
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
cellchat@DB <- CellChatDB.use
# Modify the cell order in the cellchat object
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
# Compute interaction network
cellchat <- computeCommunProb(cellchat, raw.use = FALSE)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)


saveRDS(cellchat, file = "C:\\Users\\Patrick Niekamp\\Desktop\\Experiments\\Exp. 57_Cell paper scRNA seq\\cellchat_results.rds")
saveRDS(cellchat_raw, file = "C:\\Users\\Patrick Niekamp\\Desktop\\Experiments\\Exp. 57_Cell paper scRNA seq\\cellchat_results_raw.rds")
cellchat <- readRDS("C:\\Users\\Patrick Niekamp\\Desktop\\Experiments\\Exp. 57_Cell paper scRNA seq\\cellchat_results.rds")



source <- c("Adipo-MSC", "THY1+ MSC", "Fibro-MSC","Osteo-Fibro MSC","Osteo-MSC","Osteoblast","VSMC", "AEC", "SEC" )
target <- c("GMP", "Late Myeloid", "Neutrophil", "Monocyte", "Cycling DCs", "pDC", "Macrophages", "Monocytes", "CLP")


# Fig. 3A
netVisual_individual(cellchat, signaling = "CSF",  remove.isolate = TRUE, pairLR.use = "CSF1_CSF1R", 
                     sources.use = sources.use, targets.use = target, vertex.receiver = vertex.receiver)


# Fig 3B (top)
clusters_subset <- c("Adipo-MSC", "THY1+ MSC", "Fibro-MSC","Osteo-Fibro MSC","Osteo-MSC","Osteoblast","VSMC", "AEC", "SEC") 
cells_to_keep <- which(Idents(BMN) %in% clusters_subset)
CSF1_exp <- subset(BMN, cells = cells_to_keep)
VlnPlot(CSF1_exp , features = "CSF1", pt.size = 0) 

# Fig 3B (bottom)
clusters_subset <- c("GMP", "Late Myeloid", "Neutrophil", "Monocyte", "Cycling DCs", "pDC", "Macrophages", "Monocytes", "CLP")
cells_to_keep <- which(Idents(BMN) %in% clusters_subset)
CSF1R_exp <- subset(BMN, cells = cells_to_keep)
VlnPlot(CSF1R_exp, features = "CSF1R", pt.size = 0) 
