rm (list = ls())
library(Seurat)
library(sctransform)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(dplyr)
library(openxlsx)

library(future)
plan("multiprocess", workers = 2)
options(future.globals.maxSize= 10485760000)

setwd("~/Documents/VNP_project/")

datafolder <- file.path(getwd(), "RDS_objects")
outfolder <- file.path(getwd(), "output")

# LOAD DATA ---------

load(file.path(datafolder, "Clusters_integrated_mito5.Rdata"))

All_cells <-  SubsetData(
  object = Clusters,
  ident.use = c("0", "1", "2", "3", "4", "5", "6", "7", "9", "11", "12", "13", "18", "20", "21", "22", "23", "24"))

rm(Clusters)
gc()

# PROCESS DATA ---------

#1. Generate a list with two datasets: c(ds1, ds2)

H.list <- SplitObject(All_cells, split.by = "ID")

#2. SCTransform each dataset (ds1 y ds2)

H.list <- lapply(X = H.list, FUN = function(x) {
  x <- SCTransform(
    x,
    vars.to.regress = c("percent.mt"),
    verbose = FALSE)
})

#3. SelectIntegrationFeatures to extract anchor features (Note: the higher the number of features, the better will be the results)

H.features <- SelectIntegrationFeatures(
  object.list = H.list,
  nfeatures = 2000)

#4. PrepSCTIntegration on the list and anchor features
#https://satijalab.org/seurat/v3.0/future_vignette.html

H.list <- PrepSCTIntegration(
  object.list = H.list,
  anchor.features = H.features,
  verbose = FALSE)

#5. FindIntegrationAnchors (method = SCT, important)

H.anchors <- FindIntegrationAnchors(
  object.list = H.list,
  normalization.method = "SCT",
  anchor.features = H.features)

#6. IntegrateDataset (method = SCT)

H.integrated <- IntegrateData(
  anchorset = H.anchors,
  normalization.method = "SCT",
  dims = 1:30,
  verbose = FALSE)

#7. t-SNE and Clustering

PCA <- RunPCA(
  H.integrated,
  npcs = 30,
  verbose = FALSE)

UMAP <- RunUMAP(
  PCA,
  reduction = "pca",
  dims = 1:20)

Neighbors <- FindNeighbors(
  UMAP,
  reduction = "pca",
  dims = 1:20)

Clusters <- FindClusters(
  Neighbors,
  resolution = 0.5)

# SAVE DATA --------- 

save(Clusters, file = file.path(datafolder, "Clusters3_integrated_mito5.Rdata"))

# FIND ALL MARKERS ---------

Clusters@active.assay <- "SCT"

All.markers_norm <- FindAllMarkers(
  Clusters,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25)

write.xlsx(All.markers_norm, file = file.path(outfolder, "FindAllMarkers", "All.markers_Clusters3_Integration_mito5.xlsx"))
saveRDS(All.markers_norm, file = file.path(outfolder, "FindAllMarkers", "All.markers_Clusters3_Integration_mito5.rds"))