rm(list = ls())
library(Seurat)
library(matchSCore2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pathview)
library(plyr)
library(magrittr)
library(ggsci)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(plyr)
library(future)
library(extrafont)
library(SCpubr)
library(openxlsx)
library(xlsx)

setwd("~/Documents/VNP_project/")

datafolder <- file.path(getwd(), "RDS_objects")
outfolder <- file.path(getwd(), "output")

if (!dir.exists(file.path(outfolder, "Myeloid"))) {
  dir.create(file.path(outfolder, "Myeloid"), recursive = T)
  dir.create(file.path(outfolder, "Myeloid", "DE-analysis"), recursive = T)  
  dir.create(file.path(outfolder, "Myeloid", "Excel"), recursive = T)
  dir.create(file.path(outfolder, "Myeloid", "Plots", "UMAP"), recursive = T)
  dir.create(file.path(outfolder, "Myeloid", "Plots", "Violin"), recursive = T)  
  dir.create(file.path(outfolder, "Myeloid", "Plots", "Density"), recursive = T)
  dir.create(file.path(outfolder, "Myeloid", "Plots", "DotPlot"), recursive = T)
  dir.create(file.path(outfolder, "Myeloid", "Plots", "Barplot"), recursive = T)
}

# LOAD DATA ---------

load(file.path(datafolder, "Clusters2_integrated_mito5.Rdata"))
Idents(Clusters) <- "seurat_clusters"

# PROCESS DATA ---------

coord <- as.data.frame(Clusters[["umap"]]@cell.embeddings)
cells <- row.names(coord[coord$UMAP_1 > (5),])
Clusters <- subset(Clusters, cells = cells, invert = T)

new.cluster.ids <- c(
  "Classical",
  "Intermediate",
  "Non-Classical",
  "Int Cytokine",
  "CD16+DC",
  "DC2 CD1C+",
  "Classical") 

names(new.cluster.ids) <- levels(Clusters)
Clusters <- RenameIdents(
  Clusters,
  new.cluster.ids)
Clusters[["Subclusters"]] <- Idents(Clusters)
Idents(Clusters) <- "Subclusters"
Clusters@active.assay <- "SCT"

save(Clusters, file = file.path(datafolder, "myeloid_annotated.Rdata"))

# PLOTS ---------

cluster_cols <- c(
  "Classical" = "#ECA809", # Marigold.
  "Intermediate" = "#043362", # Prussian Blue.
  "Non-Classical" = "#BC5210", # Burnt Orange.
  "Int Cytokine" = "#279185", # Celadon Green.
  "CD16+DC" = "#7EB356", # Bud Green.
  "DC2 CD1C+" = "#AC70FF" # Medium Purple.
)   

features = c(
  "S100A8","S100A12","LYZ","CD14","NEAT1","NAMPT","CSF3R","IL1B",
  "FCGR3A","MS4A7","MTSS1","TCF7L2","SIGLEC10","HLA-DRA","HLA-DPB1",
  "LGALS2","FCER1A","CLEC10A","CD1C","IFITM2","FCER1G","IFITM3","AIF1")
features = unique(features)

## UMAP ---------

pdf(file.path(outfolder, "Myeloid", "Plots", "UMAP", "Monocytes.pdf"), width = 4.675, height = 5.3)
do_DimPlot(Clusters, colors.use = cluster_cols, plot.title = "Monocytes")
dev.off()

## VIOLIN PLOT ---------

for (i in 1:length(features)) {
  pdf(file.path(outfolder, "Myeloid", "Plots", "Violin", paste0("Monocytes_VlnPlot_",features[i], ".pdf", sep = "")))
  print(do_VlnPlot(Clusters,features = features[i] ,colors.use = cluster_cols, group.by = "Subclusters"))
  dev.off()
}

## DENSITY PLOT ---------

for (i in 1:length(features)) {
  pdf(file.path(outfolder, "Myeloid", "Plots", "Density", paste0("Monocytes_Density_",features[i], ".pdf", sep = "")))
  print(do_NebulosaPlot(Clusters,features = features[i]))
  dev.off()
}

## DOT PLOT ---------

pdf(file.path(outfolder, "Myeloid", "Plots", "DotPlot", "Monocytes_DotPlot.pdf"), width = 6, height = 10)
do_DotPlot(Clusters, features = features, cols = cluster_cols, group.by = "Subclusters", flip = TRUE)
dev.off()

# GSEA ANALYSIS ---------

# Remove mitochondrial and ribosomal genes from downstream analysis

to.rm_rp <- rownames(Clusters@assays$SCT)[grep(("^RP"), rownames(Clusters@assays$SCT))]
to.rm_mt <- rownames(Clusters@assays$SCT)[grep(("^MT-"), rownames(Clusters@assays$SCT))]
to.rm <- c(to.rm_rp, to.rm_mt)

counts <- GetAssayData(Clusters, assay = "SCT")
counts <- counts[-(which(rownames(counts) %in% to.rm)),]
Clusters <- subset(Clusters, features = rownames(counts))

# Run GSEA for comparing Progressors ("Control") and VNPs ("VNP")

DE_Control <- list()
GO_DE_Control <- list()
DE_VNP <- list()
GO_DE_VNP <- list()

n <- length(levels(Clusters))

for (i in 1:n){
  
  Object <- subset(Clusters, idents = levels(Clusters)[i])
  Idents(Object) <- "group"
  
  DE_Control[[i]] <- FindMarkers(Object, ident.1 = "Control", ident.2 = "VNP", only.pos = TRUE, logfc.threshold = 0.1)
  DE_Control[[i]]  %<>% filter(p_val_adj < 0.05)
  names(DE_Control)[[i]] <- levels(Clusters)[i]
  ego <- enrichGO(rownames(DE_Control[[i]]), ont = "ALL", OrgDb = org.Hs.eg.db, keyType = "SYMBOL", minGSSize = 10)
  GO_DE_Control[[i]] <- as.data.frame(ego)
  names(GO_DE_Control)[[i]] <- levels(Clusters)[i]
  
  DE_VNP[[i]] <- FindMarkers(Object, ident.1 = "VNP", ident.2 = "Control", only.pos = TRUE, logfc.threshold = 0.1)
  DE_VNP[[i]]  %<>% filter(p_val_adj < 0.05)
  names(DE_VNP)[[i]] <- levels(Clusters)[i]
  ego <- enrichGO(rownames(DE_VNP[[i]]), ont = "ALL", OrgDb = org.Hs.eg.db, keyType = "SYMBOL", minGSSize = 10)
  GO_DE_VNP[[i]] <- as.data.frame(ego)
  names(GO_DE_VNP)[[i]] <- levels(Clusters)[i]

}

# Save GSEA tables

for (i in 1:n){
  write.xlsx(DE_Control[i], file = file = file.path(outfolder, "Myeloid", "DE-analysis", "DE_Control_Myeloid.xlsx"), sheetName = paste(i), append = T)
  write.xlsx(DE_VNP[i], file = file = file.path(outfolder, "Myeloid", "DE-analysis", "DE_VNP_Myeloid.xlsx"), sheetName = paste(i), append = T)
}

# Generate GO barplots per cell type and condition

for (i in 1:n){
  ego2 <- as.data.frame(GO_DE_Control[[i]])
  ego2 <- head(ego2, n = 20)
  p <- ggplot(ego2, aes(x = Description, y = -log10(p.adjust))) + 
              geom_bar(stat = "identity",
                position = "dodge", 
                fill="steelblue") + 
              ggtitle(paste(names(GO_DE_Control[i]), "Control", sep = " ")) +
              coord_flip() +
              theme_bw() +
              theme(panel.border = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank())
  
  pdf(file = file.path(outfolder, "Myeloid", "Plots", "Barplot", paste0("barplot_Control_", names(GO_DE_Control[i]), ".pdf")), width = 12, height = 3.5)
  print(p)
  dev.off()

  write.xlsx(GO_DE_Control[i], file = file.path(outfolder, "Myeloid", "DE-analysis", "GO_DE_Control_Myeloid.xlsx"), sheetName = paste(i), append = T)
}

for (i in 1:n){
  ego2 <- as.data.frame(GO_DE_VNP[[i]])
  ego2 <- head(ego2, n = 20)
  p <- ggplot(ego2, aes(x = Description, y = -log10(p.adjust))) +
              geom_bar(stat = "identity",
                position = "dodge",
                fill="steelblue") + 
              ggtitle(paste(names(GO_DE_VNP[i]), "VNP", sep = " ")) +
              coord_flip() +
              theme_bw()+ 
              theme(panel.border = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank()) 
    
  pdf(file = file.path(outfolder, "Myeloid", "Plots", "Barplot", paste0("barplot_VNP_", names(GO_DE_VNP[i]), ".pdf")), width = 6,height = 3.5)
  print(p)
  dev.off()
  write.xlsx(GO_DE_VNP[i], file =  file.path(outfolder, "Myeloid", "DE-analysis", "GO_DE_VNP_Myeloid.xlsx"), sheetName = paste(i), append = T)
}

# Save differentially expressed genes per cluster

Myeloid_markers <- FindAllMarkers(
  Clusters,
  assay = "SCT",
  only.pos = TRUE,
  logfc.threshold = 0.1)

Myeloid_markers  %<>% filter(p_val_adj < 0.05)

write.xlsx(Myeloid_markers, file = file.path(outfolder, "Myeloid", "DE-analysis", "Myeloid_AllMarkers.xlsx"))