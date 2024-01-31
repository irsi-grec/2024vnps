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

datafolder <- file.path(getwd(), "RDS_objects")
outfolder <- file.path(getwd(), "output")

if (!dir.exists(file.path(outfolder, "CD8_NK"))) {
  dir.create(file.path(outfolder, "CD8_NK"), recursive = T)
  dir.create(file.path(outfolder, "CD8_NK", "DE-analysis"), recursive = T)
  dir.create(file.path(outfolder, "CD8_NK", "Excel"), recursive = T)
  dir.create(file.path(outfolder, "CD8_NK", "Plots", "UMAP"), recursive = T)
  dir.create(file.path(outfolder, "CD8_NK", "Plots", "Violin"), recursive = T)  
  dir.create(file.path(outfolder, "CD8_NK", "Plots", "Density"), recursive = T)
  dir.create(file.path(outfolder, "CD8_NK", "Plots", "DotPlot"), recursive = T)
  dir.create(file.path(outfolder, "CD8_NK", "Plots", "Barplot"), recursive = T)
}

# LOAD DATA ---------

load(file.path(datafolder, "Clusters3_CD8-NK_integrated_mito5.Rdata"))
Idents(Clusters) <- "seurat_clusters"

# PROCESS DATA ---------

new.cluster.ids <- c(
  "Effector",
  "NK cytolysis",
  "Naive",
  "Terminal Effector",
  "Activated",
  "Central Memory",
  "Gamma Delta",
  "Activated",
  "Effector",
  "Transitional Mem",
  "NK",
  "Transitional Mem",
  "RBC",
  "NK response to cytokines",
  "Stem Cell Memory")

names(new.cluster.ids) <- levels(Clusters)
Clusters <- RenameIdents(
  Clusters,
  new.cluster.ids)
Clusters[["Subclusters"]] <- Idents(Clusters)
Idents(Clusters) <- "Subclusters"
Clusters@active.assay <- "SCT"

Clusters = subset(
  Clusters,
  idents = new.cluster.ids)

save(Clusters, file = file.path(datafolder, "CD8-NK_annotated.Rdata"))

# PLOTS ---------

cluster_cols <- c(
  "Effector" = "#ECA809", # Marigold.
  "NK cytolysis" = "#043362", # Prussian Blue.
  "Naive" = "#5F0F40", # Prussian Blue.
  "Terminal Effector" = "#009FF5", # Carolina Blue.
  "Activated" = "#BC5210", # Burnt Orange.
  "Central Memory" = "#279185", # Celadon Green.
  "Gamma Delta" = "#7EB356", # Bud Green.
  "Transitional Mem" = "#AC70FF", # Medium Purple.
  "NK" = "#63412C", # Van Dyke Brown.
  "NK response to cytokines" = "#D6D6D6", # Light grey.
  "Stem Cell Memory" = "#544B81") # Tyrian Purple.

features <- c(
  "CCR7","LEF1","SELL","FOXP1","KLRB1","TRGC1","TRGC2","HOPX","IL2RB","NFKBIA","KLRC1",
  "XCL2","IFITM3","GNLY","FCGR3A","KLRF1","KLRD1","KLRC2","KLRB1","GZMB","FCGR3A",
  "GZMA","CD63","FGFBP2","GNLY","ZEB2","GZMH","CCL4","CCL3","NKG7","S100A4","GZMH",
  "IL7R","ITGB1","ITGAL","ZEB2","CD74","HLA-DRB1","HLA-DRA","CCL5","IL7R", "LTB","CD27","FOS",
  "SATB1", "CD27","NOSIP")
features = unique(features)

## UMAP ---------

pdf(file.path(outfolder, "Plots", "UMAP", "CD8_NK.pdf"), width = 5, height = 5.7 )
do_DimPlot(Clusters, colors.use = cluster_cols, plot.title = "CD8 NK", legend.ncol = 3)
dev.off()

## VIOLIN PLOT ---------

for (i in 1:length(features)) {
  pdf(file.path(outfolder, "Plots", "Violin", paste0("CD8_NK_VlnPlot_",features[i], ".pdf", sep = "")))
  print(do_VlnPlot(Clusters,features = features[i] ,colors.use = cluster_cols, group.by = "Subclusters"))
  dev.off()
}

## DENSITY PLOT ---------

for (i in 1:length(features)) {
  pdf(file.path(outfolder, "Plots", "Density", paste0("CD8_NK_Density_",features[i], ".pdf", sep = "")))
  print(do_NebulosaPlot(Clusters,features = features[i]))
  dev.off()
}

## DOT PLOT ---------

pdf(file.path(outfolder, "Plots", "DotPlot", "CD8_NK_DotPlot.pdf"), width = 6, height = 10)
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
  write.xlsx(DE_Control[i], file = file = file.path(outfolder, "DE-analysis", "DE_Control_CD8_NK.xlsx"), sheetName = paste(i), append = T)
  write.xlsx(DE_VNP[i], file = file = file.path(outfolder, "DE-analysis", "DE_VNP_CD8_NK.xlsx"), sheetName = paste(i), append = T)
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
  
  pdf(file = file.path(outfolder, "Plots", "Barplot", paste0("barplot_Control_", names(GO_DE_Control[i]), ".pdf")), width = 12, height = 3.5)
  print(p)
  dev.off()
  
  write.xlsx(GO_DE_Control[i], file = file.path(outfolder, "DE-analysis", "GO_DE_Control_CD8_NK.xlsx"), sheetName = paste(i), append = T)
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
  
  pdf(file = file.path(outfolder, "Plots", "Barplot", paste0("barplot_VNP_", names(GO_DE_VNP[i]), ".pdf")), width = 6,height = 3.5)
  print(p)
  dev.off()
  write.xlsx(GO_DE_VNP[i], file =  file.path(outfolder, "DE-analysis", "GO_DE_VNP_CD8.xlsx"), sheetName = paste(i), append = T)
}

# Save differentially expressed genes per cluster

CD8_markers <- FindAllMarkers(
  Clusters,
  assay = "SCT",
  only.pos = TRUE,
  logfc.threshold = 0.1)

CD8_markers  %<>% filter(p_val_adj < 0.05)

write.xlsx(CD8_markers, file = file.path(outfolder, "DE-analysis", "CD8-NK_AllMarkers.xlsx"))

