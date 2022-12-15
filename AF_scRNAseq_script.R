Sys.time()
set.seed(365)

library(tidyverse)
library(dplyr)
library(ggtext)
library(patchwork)
library(Matrix)
library(umap)
library(DESeq2)
library(reshape2)
library(ggplot2)
library(Rtsne)
library(ggtext)
library(patchwork)
library(viridis)
library(ragg)
library(Seurat)
library(scales)
library(fgsea)
library(ggridges)
library(SingleR)

sessionInfo()

num_DEGS_to_use = 2000

## Set paths
  base_directory = ".../"
  matrix_base_directory_batch1 = paste0(base_directory, "Batch_1/Processed Data/")
  matrix_base_directory_batch2 = paste0(base_directory, "Batch_2/Processed Data/")
  matrix_base_directory_batch3 = paste0(base_directory, "Batch_3/Processed Data/")
  matrix_base_directory_batch4 = paste0(base_directory, "Batch_4/Processed Data/")
working_directory = paste0(base_directory, "R_scripts/Batch2/scRNAseq_analysis/")

gc()
############################################## Import raw SC data
lapply(c(1:10), function(i){
  
  if (i==1) {raw_sample_name <- "687"
  matrix_directory = paste0(matrix_base_directory_batch1, "687/687/outs/filtered_feature_bc_matrix/")}
  if (i==2) {raw_sample_name <- "689"
  matrix_directory = paste0(matrix_base_directory_batch1, "689/689/outs/filtered_feature_bc_matrix/")}
  if (i==3) {raw_sample_name <- "696"
  matrix_directory = paste0(matrix_base_directory_batch1, "696/696/outs/filtered_feature_bc_matrix/")}
  if (i==4) {raw_sample_name <- "Combo"
  matrix_directory = paste0(matrix_base_directory_batch1, "Combo/Combo/outs/filtered_feature_bc_matrix/")}
  if (i==5) {raw_sample_name <- "713"
  matrix_directory = paste0(matrix_base_directory_batch2, "713/713/outs/filtered_feature_bc_matrix/")}
  if (i==6) {raw_sample_name <- "719"
  matrix_directory = paste0(matrix_base_directory_batch2, "719/719/outs/filtered_feature_bc_matrix/")}
  if (i==7) {raw_sample_name <- "725"
  matrix_directory = paste0(matrix_base_directory_batch2, "725/725/outs/filtered_feature_bc_matrix/")}
  if (i==8) {raw_sample_name <- "741"
  matrix_directory = paste0(matrix_base_directory_batch3, "741/741/outs/filtered_feature_bc_matrix/")}
  if (i==9) {raw_sample_name <- "758"
  matrix_directory = paste0(matrix_base_directory_batch4, "758/758/outs/filtered_feature_bc_matrix/")}
  if (i==10) {raw_sample_name <- "744"
  matrix_directory = paste0(matrix_base_directory_batch4, "744/744/outs/filtered_feature_bc_matrix/")}
  
  barcode.path <- paste0(matrix_directory, "barcodes.tsv.gz")
  features.path <- paste0(matrix_directory, "features.tsv.gz")
  matrix.path <- paste0(matrix_directory, "matrix.mtx.gz")
  
  matrix <- readMM(file = matrix.path)
  feature.names = read.delim(features.path,
                             header = FALSE,
                             stringsAsFactors = FALSE)
  barcode.names = read.delim(barcode.path,
                             header = FALSE,
                             stringsAsFactors = FALSE)
  
  colnames(matrix) = barcode.names$V1
  rownames(matrix) = feature.names$V2
  
  min_cells <- 3
  number_of_PCAs <<- 20
  
  seurat_object <- CreateSeuratObject(counts = matrix, project = "scRNAseq_1stBatch", min.cells=min_cells, min.features = 1)
  
  setwd(paste0(working_directory, "Quality_Control", sep="/", raw_sample_name))
  
  seurat_object@meta.data[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "MT-")[,1]
  seurat_object@meta.data[["HSP"]] <- PercentageFeatureSet(seurat_object, pattern = "HSP")[,1]
  seurat_object@meta.data[["Ribo1"]] <- (PercentageFeatureSet(seurat_object, pattern = "RPL")[,1])
  seurat_object@meta.data[["Ribo2"]] <- (PercentageFeatureSet(seurat_object, pattern = "RPS")[,1])
  seurat_object@meta.data[["HB"]] <- PercentageFeatureSet(seurat_object, pattern = "HB")[,1]
  
  #define genes that drive lots of technical variance & remove from sample
  hkGeneREGEX=c("DNAJ", "EIF", "RPL", "RPS", "RPN1", "POLR", "SNX", "HSP", "H1FX", 
                "H2AF", "PRKA", "NDUF", "ATP", "PSM", "UBA", "UBE", "USP", "TXN")
  coreExcludeGenes = unlist(unique(c(grep('\\.[0-9]+',rownames(seurat_object),value=TRUE), #Poorly characterised
                                     grep('MALAT1',rownames(seurat_object),value=TRUE), #Contamination or highly expressed poorly characterised
                                     grep('^MT-',rownames(seurat_object),value=TRUE), #Mitochondria
                                     lapply(1:length(hkGeneREGEX), function(e){
                                       grep(hkGeneREGEX[e], rownames(seurat_object), value = TRUE)
                                     })
  ))
  )
  filtered_genes <- c(rownames(seurat_object), coreExcludeGenes)
  filtered_genes <- filtered_genes[!duplicated(filtered_genes)]
  seurat_object <- subset(seurat_object, features = filtered_genes)
  
  ### Measure doublet probability (mapped upon combination of objects) - https://bioconductor.org/books/3.14/OSCA.advanced/doublet-detection.html#doublet-simulation
    genes <- as.data.frame(rownames(seurat_object@assays[["RNA"]]@data))
    colnames(genes) <- "gene_short_name"
    rownames(genes) <- genes[,1]
    doublet_detection <- SeuratWrappers::as.cell_data_set(seurat_object)
    
    doublet_density <- scDblFinder::computeDoubletDensity(doublet_detection)
    summary(doublet_density)
    seurat_object[["DoubletScore"]] <- doublet_density
  
  # Extract into seperate matrix
  metadata <- seurat_object@meta.data
    metadata <- metadata[,c("nCount_RNA", "nFeature_RNA", "percent.mt", "DoubletScore")]
  
  # Feature QC Plots
  pdf("FeatureScatterplots_preFiltering.pdf")
    print(ggplot() + geom_point(data = metadata, aes(x = nCount_RNA, y = nFeature_RNA, color = DoubletScore)))
  dev.off()
  
  # Violin QC Plots
  pdf("ViolinPlots_preFiltering.pdf")
    print(VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "DoubletScore"), ncol = 2))
  dev.off()
  
  ### Filter data 
  seurat_object <- subset(seurat_object, subset = nFeature_RNA > 150 & nFeature_RNA < 8000 & percent.mt < 10)
  
  
  # Re-extract metadata
    metadata <- seurat_object@meta.data
    metadata <- metadata[,c("nCount_RNA", "nFeature_RNA", "percent.mt", "DoubletScore")]
  
  write.table(metadata, file = paste0("metadata", raw_sample_name, "_postfiltering.txt"), sep = "\t",
              row.names = TRUE, col.names = TRUE)
  
  # Feature QC Plots
  pdf("FeatureScatterplots_postFiltering.pdf")
    print(ggplot() + geom_point(data = metadata, aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt)))
    print(ggplot() + geom_point(data = metadata, aes(x = nCount_RNA, y = nFeature_RNA, color = DoubletScore)))
  dev.off()
  
  # Violin QC Plots
  pdf("ViolinPlots_postFiltering.pdf")
    print(VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "DoubletScore"), ncol = 2))
  dev.off()
  
  metadata$Barcode <- rownames(metadata)
  metadata$Barcode <- paste0(raw_sample_name, "_", metadata$Barcode)
  
  ## Prepare individual data for integration
  seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat_object <- FindVariableFeatures(object = seurat_object, selection.method="vst", nfeatures = num_DEGS_to_use)
  DefaultAssay(seurat_object) <- "RNA"
  
  ## Save data to global environment and rename
  assign((paste0("seurat_object", sep="_", raw_sample_name)), seurat_object, envir=.GlobalEnv)
  assign((paste0("metadata", sep="_", raw_sample_name)), metadata, envir=.GlobalEnv)
  message(paste0(raw_sample_name, " complete"))
})

metadata <- rbind(metadata_687, metadata_689, metadata_696, metadata_713, metadata_Combo, 
                  metadata_719, metadata_725, metadata_741, metadata_758, metadata_744)

  setwd(paste0(working_directory, "Quality_Control"))
  write.table(metadata, file = "metadata_combined.txt", sep = "\t",
              row.names = TRUE, col.names = TRUE)

rownames(metadata) <- metadata$Barcode



#### Split combo lane
newID <- read.table(file=paste0(base_directory, "R_scripts/Input_Data/donor_ids_ComboLane21.tsv"), header = T)
newID <- split(newID, f=newID$best_singlet)

# Combo = 680 (10k cells) + 696 (16k cells) 695 (2.8k cells)
barcodes680 <- newID$donor0
barcodes680 <- barcodes680$cell
barcodes695 <- newID$donor1
barcodes695 <- barcodes695$cell
barcodes696 <- newID$donor2
barcodes696 <- barcodes696$cell

## Extrat data from combo object
seurat_object_680 <- subset(seurat_object_Combo, cells = barcodes680)
seurat_object_695 <- subset(seurat_object_Combo, cells = barcodes695)
seurat_object_696new <- subset(seurat_object_Combo, cells = barcodes696)

seurat_object_696 <- merge(seurat_object_696, seurat_object_696new)

#### Combine my data
seurat_object_687@meta.data[["orig_sample"]] <- "32PCW_687"
seurat_object_689@meta.data[["orig_sample"]] <- "17PCW_689"
seurat_object_696@meta.data[["orig_sample"]] <- "30PCW_696"
seurat_object_680@meta.data[["orig_sample"]] <- "17PCW_680"
seurat_object_695@meta.data[["orig_sample"]] <- "13PCW_695"
seurat_object_713@meta.data[["orig_sample"]] <- "14PCW_713"
seurat_object_719@meta.data[["orig_sample"]] <- "17PCW_719"
seurat_object_725@meta.data[["orig_sample"]] <- "19PCW_725"
seurat_object_741@meta.data[["orig_sample"]] <- "28PCW_741"
seurat_object_758@meta.data[["orig_sample"]] <- "28PCW_758"


#### Split 744 lane - Zombie?
newID <- read.table(file=paste0(base_directory, "R_scripts/Input_Data/donor_ids_744.tsv"), header = T)
newID <- split(newID, f=newID$donor_id)

# 744 = 696(1) (343 cells) + 744 (358 cells) 689 (243 cells)
barcodes744_696 <- newID$donor0
barcodes744_696 <- barcodes744_696$cell
barcodes744_744 <- newID$donor1
barcodes744_744 <- barcodes744_744$cell
barcodes744_689 <- newID$donor2
barcodes744_689 <- barcodes744_689$cell

## Extrat data from combo object
seurat_object_744 <- subset(seurat_object_744, cells = barcodes744_744)
seurat_object_744@meta.data[["orig_sample"]] <- "25PCW_744"

#### Batch correction & integration
Batch1_seurat_object <- merge(seurat_object_687, seurat_object_689,  
                              merge.data = TRUE, add.cell.ids = c("687", "689"))
Batch2_seurat_object <- merge(seurat_object_696, y=c(seurat_object_680, seurat_object_695),  
                              merge.data = TRUE, add.cell.ids = c("696", "680", "695"))
Batch3_seurat_object <- seurat_object_713
Batch3_seurat_object <- RenameCells(Batch3_seurat_object, add.cell.id = "713")
Batch4_seurat_object <- seurat_object_719
Batch4_seurat_object <- RenameCells(Batch4_seurat_object, add.cell.id = "719")
Batch5_seurat_object <- seurat_object_725
Batch5_seurat_object <- RenameCells(Batch5_seurat_object, add.cell.id = "725")
Batch6_seurat_object <- seurat_object_741 # excl seurat_object_743
Batch6_seurat_object <- RenameCells(Batch6_seurat_object, add.cell.id = "741")
Batch7_seurat_object <- merge(seurat_object_744, seurat_object_758,  
                              merge.data = TRUE, add.cell.ids = c("744", "758"))

seurat_merged <- merge(Batch1_seurat_object, y=c(Batch2_seurat_object, Batch3_seurat_object, 
                                                   Batch4_seurat_object, Batch5_seurat_object,
                                                   Batch6_seurat_object, Batch7_seurat_object), 
                         merge.data = TRUE, project = "AF_SC")
dim(seurat_merged)
number_of_cells_analysed=length(colnames(seurat_merged))

seurat_combined_list <- list()
seurat_combined_list[["Batch1"]] <- Batch1_seurat_object
seurat_combined_list[["Batch2"]] <- Batch2_seurat_object
seurat_combined_list[["Batch3"]] <- Batch3_seurat_object
seurat_combined_list[["Batch4"]] <- Batch4_seurat_object
seurat_combined_list[["Batch5"]] <- Batch5_seurat_object
seurat_combined_list[["Batch6"]] <- Batch6_seurat_object
seurat_combined_list[["Batch7"]] <- Batch7_seurat_object

  features <- SelectIntegrationFeatures(object.list = seurat_combined_list)
  anchors <- FindIntegrationAnchors(object.list = seurat_combined_list, anchor.features = features)
  
  saveRDS(anchors, paste0(working_directory, "full_SC_anchors_", ceiling((number_of_cells_analysed)/1000),"k.RDS"))
  # anchors <- readRDS(paste0(working_directory, "full_SC_anchors_", ceiling((number_of_cells_analysed)/1000),"k.RDS"))
  
  seurat_combined <- IntegrateData(anchorset = anchors)
  DefaultAssay(seurat_combined) <- "integrated"
  rm(seurat_combined_list, anchors)


rownames(metadata) <- metadata$Barcode
seurat_object <- seurat_combined


saveRDS(seurat_object, paste0(working_directory, "CombiningDataSetsseurat_object_normalised_", ceiling(number_of_cells_analysed/1000),"k.rds"))
seurat_object_normalised = seurat_object


#### QC
setwd(paste0(working_directory, "Quality_Control"))
seurat_object@meta.data[["DoubletScore"]] <- metadata[colnames(seurat_object),"DoubletScore"]

# Feature QC Plots
pdf("combined_FeatureScatterplots.pdf")
ggplot() + geom_point(data = metadata, aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt))
ggplot() + geom_point(data = metadata, aes(x = nCount_RNA, y = nFeature_RNA, color = DoubletScore))
dev.off()

# Violin QC Plots
pdf("combined_ViolinPlots.pdf")
VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "DoubletScore"), ncol = 2)
dev.off()







setwd(paste0(working_directory, "Quality_Control"))
# Feature QC Plots
pdf("post-normalisation_combined_FeatureScatterplots.pdf")
ggplot() + geom_point(data = metadata, aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt))
dev.off()

# Violin QC Plots
pdf("post-normalisation_combined_ViolinPlots.pdf")
VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()




#### Checking for female gender cells
## Checking single cell co-expression of main female gene = XIST
rm(female_genes)
female_genes <- seurat_object@assays[["RNA"]]@counts["XIST",]
female_genes <- as.data.frame(female_genes)
colnames(female_genes) <- "XIST"

XIST_pos <- filter(female_genes, female_genes$XIST!=0)
XIST_neg <- filter(female_genes, female_genes$XIST==0)

write.table(female_genes, file = "Maternal_cells.txt", sep = "\t",
            row.names = TRUE, col.names = TRUE)


#### Data summarisation
## Select top 2,000 Highly Variable Genes (HGV)
seurat_object <- FindVariableFeatures(object = seurat_object, selection.method="vst", nfeatures = num_DEGS_to_use)
top20 <- head(VariableFeatures(seurat_object), 20)

setwd(working_directory)
pdf("VariableFeaturesPlot.pdf")
plot1HGV <- VariableFeaturePlot(seurat_object)
plot1HGV
LabelPoints(plot = plot1HGV, points = top20, repel = TRUE)
dev.off()





#### Assessing no. of PCAs required for Dimensionality Reduction - Parameters
## Save list of all genes detected
all.genes <- rownames(seurat_object)

setwd(working_directory)
write.table(all.genes, file = "Seurat_COMPLETE_Gene_list.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE)

## Scaling expression before PCAs
seurat_object <- ScaleData(seurat_object, features = all.genes, do.center = TRUE, do.scale = TRUE)

## Run PCA
seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object),npcs = 50)

# Plots
setwd(working_directory)
pdf("Dimensionality_plots.pdf", width=30, height=30)
VizDimLoadings(seurat_object, dims = 1:20, reduction = "pca")
DimPlot(seurat_object, reduction = "pca")
DimHeatmap(seurat_object, dims = 1:20, cells = 500, balanced = TRUE)
ElbowPlot(seurat_object,ndims=50)
dev.off()

## From elbow plot, what is correct number of PCAs?
number_of_PCAs = 20






#### Clustering Analysis - Parameters
setwd(working_directory)

## Dimensionality reduction
seurat_object <- RunUMAP(seurat_object, dims = 1:number_of_PCAs)
seurat_object <- RunTSNE(object = seurat_object, check_duplicates = FALSE)

# Set the resolution parameter based on modularity optimization, the parameter is correlated to the scale of observing communities.
# The higher the resolution parameter, the larger the number of smaller communities.
seurat_object <- FindNeighbors(seurat_object, dims = 1:number_of_PCAs)
seurat_object <- FindClusters(seurat_object, resolution = 0.95)

head(Idents(seurat_object), 5)

setwd(working_directory)
## Save plots
pdf("Clustering_plots.pdf")
DimPlot(seurat_object, reduction = "umap")
DimPlot(seurat_object, reduction = "tsne")
dev.off()

saveRDS(seurat_object, paste0(working_directory, "combined_complete_seurat_object_", ceiling(number_of_cells_analysed/1000),"k.rds"))
seurat_object <- readRDS(paste0(working_directory, "combined_complete_seurat_object_", ceiling(number_of_cells_analysed/1000),"k.rds"))
metadata <- seurat_object@meta.data





#### UMAP QC Plots
UMAPdata=as.data.frame(seurat_object@meta.data)
UMAPdata$UMAP_1 <- seurat_object@reductions[["umap"]]@cell.embeddings[,"UMAP_1"]
UMAPdata$UMAP_2 <- seurat_object@reductions[["umap"]]@cell.embeddings[,"UMAP_2"]
UMAPdata$Cluster=seurat_object@meta.data$seurat_clusters
UMAPdata$Barcode <- rownames(UMAPdata)
seurat_metadata_clusters <- UMAPdata


setwd(paste0(working_directory, "Quality_Control"))
pdf("post_filtering_UMAP_QC.pdf")
ggplot() + geom_point(data = UMAPdata, aes(x = UMAP_1, y = UMAP_2, color = percent.mt))
ggplot() + geom_point(data = UMAPdata, aes(x = UMAP_1, y = UMAP_2, color = nFeature_RNA))
ggplot() + geom_point(data = UMAPdata, aes(x = UMAP_1, y = UMAP_2, color = nCount_RNA))
dev.off()

pdf("doublet_scoring_UMAP_QC.pdf")
ggplot() + geom_point(data = UMAPdata, aes(x = UMAP_1, y = UMAP_2, color = DoubletScore))
dev.off()

## Set trimesters
second_trimester_samples <- c("13PCW_695", "14PCW_713", "14PCW_743_cKIT", 
                              "30PCW_696", "17PCW_680", "17PCW_689", "17PCW_719", 
                              "19PCW_725")
lapply(1:length(second_trimester_samples), function(y){
UMAPdata[(UMAPdata$orig_sample==second_trimester_samples[y]),"Trimester"] <<- "Second"
})
second_trimester_samples <- c("25PCW_744", "28PCW_758", "28PCW_741", 
                              "32PCW_687")
lapply(1:length(second_trimester_samples), function(y){
  UMAPdata[(UMAPdata$orig_sample==second_trimester_samples[y]),"Trimester"] <<- "Third"
})

barcodes680 <- (UMAPdata[(UMAPdata$orig_sample=="17PCW_680"),"Barcode"])
barcodes695 <- (UMAPdata[(UMAPdata$orig_sample=="13PCW_695"),"Barcode"])
barcodes696 <- (UMAPdata[(UMAPdata$orig_sample=="30PCW_696"),"Barcode"])
barcodes687 <- (UMAPdata[(UMAPdata$orig_sample=="32PCW_687"),"Barcode"])
barcodes689 <- (UMAPdata[(UMAPdata$orig_sample=="17PCW_689"),"Barcode"])
barcodes713 <- (UMAPdata[(UMAPdata$orig_sample=="14PCW_713"),"Barcode"])
barcodes719 <- (UMAPdata[(UMAPdata$orig_sample=="17PCW_719"),"Barcode"])
barcodes725 <- (UMAPdata[(UMAPdata$orig_sample=="19PCW_725"),"Barcode"])
barcodes741 <- (UMAPdata[(UMAPdata$orig_sample=="28PCW_741"),"Barcode"])
barcodes744 <- (UMAPdata[(UMAPdata$orig_sample=="25PCW_744"),"Barcode"])
barcodes758 <- (UMAPdata[(UMAPdata$orig_sample=="28PCW_758"),"Barcode"])

UMAPdata$week <- as.numeric(substr(UMAPdata$orig_sample, 1,2))+2

colour_plot <- ggplot(UMAPdata, aes(x = UMAP_1, y = UMAP_2, color=orig_sample)) +
  scale_color_viridis(option="turbo", discrete=TRUE)
my_colours <- unique((ggplot_build(colour_plot))[["data"]][[1]][["colour"]])

setwd(paste0(working_directory, "Quality_Control"))
pdf("Batch_Differences_UMAP_QC.pdf", width=9, height=7)
ggplot() + geom_point(data = UMAPdata, size = 0.4, shape = 16,  
                      aes(x = UMAP_1, y = UMAP_2, color = orig_sample)) + 
  scale_color_manual(values=my_colours, name="Sample") +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme_light(base_family = "sans",base_line_size = 0.5,base_rect_size = 0.1,base_size = 14)

ggplot() + geom_point(data = UMAPdata, size = 0.4, shape = 16,  
                      aes(x = UMAP_1, y = UMAP_2, color = Trimester)) + 
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme_light(base_family = "sans",base_line_size = 0.5,base_rect_size = 0.1,base_size = 14)

ggplot() + geom_point(data = UMAPdata, size = 0.4, shape = 16,  
                      aes(x = UMAP_1, y = UMAP_2, color = as.numeric(week))) + 
  scale_color_gradient2(mid = "orange", high="darkorchid4", name="Gestational Age") +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme_light(base_family = "sans",base_line_size = 0.5,base_rect_size = 0.1,base_size = 14)

DimPlot(object = seurat_object, cells.highlight = barcodes695, cols.highlight = "red", cols = "gray", order = TRUE, pt.size = 0.1, sizes.highlight = 0.1) + ggtitle("695 13PCW cells")
DimPlot(object = seurat_object, cells.highlight = barcodes713, cols.highlight = "red", cols = "gray", order = TRUE, pt.size = 0.1, sizes.highlight = 0.1) + ggtitle("713 14PCW cells")
DimPlot(object = seurat_object, cells.highlight = barcodes719, cols.highlight = "red", cols = "gray", order = TRUE, pt.size = 0.1, sizes.highlight = 0.1) + ggtitle("719 17PCW cells")
DimPlot(object = seurat_object, cells.highlight = barcodes680, cols.highlight = "red", cols = "gray", order = TRUE, pt.size = 0.1, sizes.highlight = 0.1) + ggtitle("680 17PCW cells")
DimPlot(object = seurat_object, cells.highlight = barcodes689, cols.highlight = "red", cols = "gray", order = TRUE, pt.size = 0.1, sizes.highlight = 0.1) + ggtitle("689 17PCW cells")
DimPlot(object = seurat_object, cells.highlight = barcodes725, cols.highlight = "red", cols = "gray", order = TRUE, pt.size = 0.1, sizes.highlight = 0.1) + ggtitle("725 19PCW cells")
DimPlot(object = seurat_object, cells.highlight = barcodes741, cols.highlight = "red", cols = "gray", order = TRUE, pt.size = 0.1, sizes.highlight = 0.1) + ggtitle("741 28PCW cells")
DimPlot(object = seurat_object, cells.highlight = barcodes696, cols.highlight = "red", cols = "gray", order = TRUE, pt.size = 0.1, sizes.highlight = 0.1) + ggtitle("696 30PCW cells")
DimPlot(object = seurat_object, cells.highlight = barcodes687, cols.highlight = "red", cols = "gray", order = TRUE, pt.size = 0.1, sizes.highlight = 0.1) + ggtitle("687 32PCW cells")
DimPlot(object = seurat_object, cells.highlight = barcodes744, cols.highlight = "red", cols = "gray", order = TRUE, pt.size = 0.1, sizes.highlight = 0.1) + ggtitle("744 25PCW cells")
DimPlot(object = seurat_object, cells.highlight = barcodes758, cols.highlight = "red", cols = "gray", order = TRUE, pt.size = 0.1, sizes.highlight = 0.1) + ggtitle("758 28PCW cells")
dev.off()





setwd(paste0(working_directory, "Quality_Control"))
write.table(metadata, file = "metadata_UMAP_combined.txt", sep = "\t",
            row.names = TRUE, col.names = TRUE)





## Standardized theme & ggplot layouts
UMAP_ggplot <- function(input, title, save_as_PDF) {
  
  ## Extract colours from the viridis rainbow colours scale package
  colour_plot <- ggplot(input, aes(x = UMAP_1, y = UMAP_2, color=(factor(Cluster)))) +
    scale_color_viridis(option="turbo", discrete=TRUE)
  my_colours <- unique((ggplot_build(colour_plot))[["data"]][[1]][["colour"]])
  my_colours <- colorspace::desaturate(my_colours, 0.01)
  
  plot <- ggplot(input, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(
      aes(color=as.factor(Cluster)),
      size=1.05,
      alpha=0.95,
      shape=16,  
      stroke=0) +
    guides(
      colour = guide_legend(override.aes = list(size=5, stroke=0)))  +  #size of legend points
    scale_color_manual(
      values = my_colours,
      name="Clusters:") +
    labs(
      x="UMAP 1",
      y="UMAP 2",
      title=title) +
    theme_light(
      base_family = "sans",
      base_line_size = 0.2,
      base_rect_size = 0.2
    ) +
    geom_label(
      data = input %>%                     # filter what data to label
        group_by(Cluster) %>%
        filter(row_number()==ceiling(n()/2)), # find middle point for each cluster
      aes(label=as.factor(Cluster), color=as.factor(Cluster)),
      fill="grey70",
      size=3,
      fontface = "bold",
      nudge_x = -1,
      nudge_y = -2,
      show.legend = FALSE,
      label.padding = unit(0.15, "lines"),
      label.size = 1,
    )
  print(plot)                       

  ## Save plot
  if (save_as_PDF == "y") {  
    pdf(paste0(working_directory, "/", title, "_", ceiling(number_of_cells_analysed/1000), "k.pdf"), height=8, width=9)
    print(plot)
    dev.off()
  }
  if (save_as_PDF == "n") {print('save_as_PDF="y" to save as PDF')}
}
# input(UMAPdata), "title", save_as_PDF("y"/"n"),
UMAP_ggplot(UMAPdata, "Standard_Clustering", "y")




## Label major cell types using SingleR
  SCE <- Seurat::as.SingleCellExperiment(seurat_object)
  
  ref <- celldex::HumanPrimaryCellAtlasData()
  pred <- SingleR::SingleR(test = SCE, ref = ref, assay.type.test=1,
                           labels = ref$label.main)
  pred <- as.data.frame(pred)
  pred$Barcode <- rownames(pred)
  
  setwd(working_directory)
  write.csv(pred, "singleR_output.csv")
  seurat_SingleR_stemcell <- subset(seurat_object, cells=pred[pred$labels=="Tissue_stem_cells","Barcode"])
  saveRDS(seurat_SingleR_stemcell, "seurat_SingleR_stemcell.rds")

setwd(working_directory)
pred <- read.csv("singleR_output.csv")
colnames(pred)[colnames(pred) == 'X'] <- 'Barcode'

pred <- pred[,c("Barcode", "pruned.labels")]

pred <- merge(pred, seurat_metadata_clusters, by="Barcode")

pred <- pred[!is.na(pred$pruned.labels),]
pred$epi_label <- "Unlabelled"
pred[pred$pruned.labels=="Epithelial_cells","epi_label"] <- "Epithelial_cells"

plot <- ggplot(pred, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(
    aes(color=epi_label),
    size=1.05,
    alpha=0.95,
    shape=16,                # shape 16 means we can ignore 'stroke'
    stroke=0) +
  guides(
    colour = guide_legend(override.aes = list(size=5, stroke=0)))  +  #size of legend points
  scale_color_manual(
    values = c("darkorange2", "gray80"),
    name="Cell Labels") +
  labs(
    x="UMAP 1",
    y="UMAP 2",
    title="Cell Labelling with SingleR") +
  theme_light(
    base_family = "sans",
    base_line_size = 0.2,
    base_rect_size = 0.2
  )
pdf(paste0(working_directory, "SingleR_cellLabelling_epithelialOnly.pdf"), height=8, width=9)
print(plot)
dev.off()










#### Epithelial barcodes - Select epithelial clusters by hand
epithelial_clusters=c(11,8,4,9,1,10,6,0,2,15,24,12,3,7,19,21) 
epithelial_seurat_object <- subset(seurat_object, idents = epithelial_clusters)
epithelial_SC_barcodes <- colnames(epithelial_seurat_object)
setwd(working_directory)
write.csv(epithelial_SC_barcodes, "epithelial_SC_barcodes_myriad.csv")



#### ggplot Violin plot of Lung Genes in epithelial cluster only when local and all when cluster
seurat_object_epithelia <- subset(seurat_object, cells=epithelial_SC_barcodes)

seurat <- seurat_object_epithelia
genes_to_analyse=as.data.frame(read.csv(paste0(base_directory, "R_scripts/Input_Data/**** Please see publication for gene list ****")), header=TRUE)
genes=genes_to_analyse$lung_markers

genes_all <- genes
genes <- genes[genes %in% rownames(seurat_object_epithelia@assays[["RNA"]]@data)]
data <- as.data.frame(seurat@assays[["RNA"]]@data[genes[1],])

for (x in c(2:(length(genes)))) {
  data[,x] <- seurat@assays[["RNA"]]@data[genes[x],]
}
colnames(data) <- genes
data$Barcode <- rownames(data)

setwd(working_directory)
write.csv(data[rowSums(data[,c(1:length(colnames(data))-1)])>0,], "Lung_gene_data_Oct22.csv")

UMAPlong <- gather_(data, "Gene", "Expression", genes)
UMAPlong[UMAPlong$Expression==0,] <- NA
UMAPlong <- UMAPlong[complete.cases(UMAPlong), ]

data_summary <- function(x) {m <- mean(x)
ymin <- m-sd(x)
ymax <- m+sd(x)
return(c(y=m,ymin=ymin,ymax=ymax))
}

setwd(working_directory)
tiff("Expression_of_Lung_progentior_markers_in_the_epithelial_AF_SC_cluster.tiff", units="in", width=10, height=3, res=150, compression = 'lzw')
ggplot(UMAPlong, aes(x=Gene, y=Expression, color=Gene, fill=Gene)) +
  geom_violin() +
  stat_summary(fun.data=data_summary, color="black") +
  theme_light(base_family = "Arial",base_line_size = 0.5,base_rect_size = 0.1,base_size = 20) +
  theme(text=element_text(color="black", size=14, face="bold"),
                          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_x_discrete(limits = genes_all) +
  ylim(c(0,6)) +
  ylab("Lung Marker Expression") + xlab(NULL) + NoLegend() +
  ggtitle(paste0("Epithelial clusters ", sum(epithelial_clusters)))
dev.off()

## Lung, expression of NKX2.1
lung_barcodes <- data[rowSums(data[,c(1:length(colnames(data))-1)])>0,]
genes=c("SOX2","SOX9","ETV4","ETV5","ID2","GATA6","MSX2","NKX2.1","FOXA2")
genes <- genes[genes %in% colnames(lung_barcodes)]
lung_barcodes_binary <- lung_barcodes
rownames(lung_barcodes_binary) <- lung_barcodes_binary$Barcode
lung_barcodes_binary <- lung_barcodes_binary[,genes]
lung_barcodes_binary[lung_barcodes_binary>0] <- 1
lung_barcodes_binary$sum <- rowSums(lung_barcodes_binary)
lung_barcodes_binary <- lung_barcodes_binary[lung_barcodes_binary$sum>=2,]
lung_barcodes <- rownames(lung_barcodes_binary)

barcodes=lung_barcodes
barcodes <- barcodes[!duplicated(barcodes)]
df <- seurat[barcodes,]

pdf(paste0(directory, "/Lung_coexpression.pdf"), height=8, width=9)
ggplot() + 
  geom_point(data=seurat, aes(x=UMAP_1, y=UMAP_2), color="gray80", size=3, shape=16) + 
  geom_point(data=df, aes(x=UMAP_1, y=UMAP_2), color="black", stroke=0.1, fill="red", shape=21, size=3) +
  theme_light(
    base_family = "sans",
    base_line_size = 0.2,
    base_rect_size = 0.2)
dev.off()




#### ggplot Violin plot of kidney Genes in epithelial cluster only when local and all when cluster
seurat_object_epithelia <- subset(seurat_object, cells=epithelial_SC_barcodes)

seurat <- seurat_object_epithelia
genes_to_analyse=as.data.frame(read.csv(paste0(base_directory, "R_scripts/Input_Data/**** Please see publication for gene list ****")), header=TRUE)
genes=genes_to_analyse$kidney_markers

genes_all <- genes
genes <- genes[genes %in% rownames(seurat_object_epithelia@assays[["RNA"]]@data)]
data <- as.data.frame(seurat@assays[["RNA"]]@data[genes[1],])

for (x in c(2:(length(genes)))) {
  data[,x] <- seurat@assays[["RNA"]]@data[genes[x],]
}
colnames(data) <- genes
data$Barcode <- rownames(data)

setwd(working_directory)
write.csv(data[rowSums(data[,c(1:length(colnames(data))-1)])>0,], "Kidney_gene_data_Oct22.csv")

UMAPlong <- gather_(data, "Gene", "Expression", genes)
UMAPlong[UMAPlong$Expression==0,] <- NA
UMAPlong <- UMAPlong[complete.cases(UMAPlong), ]

data_summary <- function(x) {m <- mean(x)
ymin <- m-sd(x)
ymax <- m+sd(x)
return(c(y=m,ymin=ymin,ymax=ymax))
}

setwd(working_directory)
tiff("Expression_of_Kidney_progentior_markers_in_the_epithelial_AF_SC_cluster.tiff", units="in", width=10, height=3, res=150, compression = 'lzw')
ggplot(UMAPlong, aes(x=Gene, y=Expression, color=Gene, fill=Gene)) +
  geom_violin() +
  stat_summary(fun.data=data_summary, color="black") +
  theme_light(base_family = "Arial",base_line_size = 0.5,base_rect_size = 0.1,base_size = 20) +
  theme(text=element_text(color="black", size=14, face="bold"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_x_discrete(limits = genes_all) +
  ylim(c(0,6)) +
  ylab("Kidney Marker Expression") + xlab(NULL) + NoLegend() +
  ggtitle(ggtitle(paste0("Epithelial clusters ", sum(epithelial_clusters))))
dev.off()


## Kidney, expression of PAX8 or LHX1
genes <- c("PAX8", "LHX1", "PAX2", "SIX2", "WT1")
kidney_barcodes <- data[rowSums(data[,c(1:length(colnames(data))-1)])>0,]
kidney_barcodes_binary <- kidney_barcodes
kidney_barcodes_binary <- kidney_barcodes_binary[!duplicated(kidney_barcodes_binary$Barcode),]
rownames(kidney_barcodes_binary) <- kidney_barcodes_binary$Barcode
kidney_barcodes_binary <- kidney_barcodes_binary[,genes]
kidney_barcodes_binary[kidney_barcodes_binary>0] <- 1
kidney_barcodes_binary$sum <- rowSums(kidney_barcodes_binary)
kidney_barcodes_binary <- kidney_barcodes_binary[kidney_barcodes_binary$sum>=2,]
kidney_barcodes <- rownames(kidney_barcodes_binary)
kidney_barcodes <- kidney_barcodes$Barcode

barcodes=kidney_barcodes
barcodes <- barcodes[!duplicated(barcodes)]
df <- seurat[barcodes,]

pdf(paste0(directory, "/Kidney_coexpression.pdf"), height=8, width=9)
ggplot() + 
  geom_point(data=seurat, aes(x=UMAP_1, y=UMAP_2), color="gray80", size=3, shape=16) + 
  geom_point(data=df, aes(x=UMAP_1, y=UMAP_2), color="black", stroke=0.1, fill="red", shape=21, size=3) +
  theme_light(
    base_family = "sans",
    base_line_size = 0.2,
    base_rect_size = 0.2)
dev.off()





#### ggplot Violin plot of intestinal Genes in epithelial cluster only when local and all when cluster
seurat_object_epithelia <- subset(seurat_object, cells=epithelial_SC_barcodes)

seurat <- seurat_object_epithelia

genes_to_analyse=as.data.frame(read.csv(paste0(base_directory, "R_scripts/Input_Data/**** Please see publication for gene list ****")), header=TRUE)
genes=genes_to_analyse$intestine_markers

genes_all <- genes
genes <- genes[genes %in% rownames(seurat_object_epithelia@assays[["RNA"]]@data)]
data <- as.data.frame(seurat@assays[["RNA"]]@data[genes[1],])

for (x in c(2:(length(genes)))) {
  data[,x] <- seurat@assays[["RNA"]]@data[genes[x],]
}
colnames(data) <- genes
data$Barcode <- rownames(data)

setwd(working_directory)
write.csv(data[rowSums(data[,c(1:length(colnames(data))-1)])>0,], "Intestinal_gene_data_Oct22.csv")

UMAPlong <- gather_(data, "Gene", "Expression", genes)
UMAPlong[UMAPlong$Expression==0,] <- NA
UMAPlong <- UMAPlong[complete.cases(UMAPlong), ]

data_summary <- function(x) {m <- mean(x)
ymin <- m-sd(x)
ymax <- m+sd(x)
return(c(y=m,ymin=ymin,ymax=ymax))
}

setwd(working_directory)
tiff("Expression_of_intestinal_progentior_markers_in_the_epithelial_AF_SC_cluster.tiff", units="in", width=10, height=3, res=150, compression = 'lzw')
ggplot(UMAPlong, aes(x=Gene, y=Expression, color=Gene, fill=Gene)) +
  geom_violin() +
  stat_summary(fun.data=data_summary, color="black") +
  theme_light(base_family = "Arial",base_line_size = 0.5,base_rect_size = 0.1,base_size = 20) +
  theme(text=element_text(color="black", size=14, face="bold"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_x_discrete(limits = genes_all) +
  ylim(c(0,6)) +
  ylab("Intestine Marker Expression") + xlab(NULL) + NoLegend() +
  ggtitle(paste0("Epithelial clusters ", sum(epithelial_clusters)))
dev.off()



intestine_barcodes <- data[rowSums(data[,c(1:length(colnames(data))-1)])>0,]
# SOX9 + one of these genes
genes=c("ASCL2","ATOH1", "DCLK1", "NGN3","MUC2", "MUC5AC", "FABP1", "FABP2", "ALPI", "LRIG1")
genes <- genes[genes %in% colnames(intestine_barcodes)]
intestine_barcodes_binary <- intestine_barcodes
rownames(intestine_barcodes_binary) <- intestine_barcodes_binary$Barcode
intestine_barcodes_binary <- intestine_barcodes_binary[intestine_barcodes_binary$SOX9>0,]
intestine_barcodes_binary <- intestine_barcodes_binary[,genes]
intestine_barcodes_binary[intestine_barcodes_binary>0] <- 1
intestine_barcodes_binary$sum <- rowSums(intestine_barcodes_binary)
intestine_barcodes_binary <- intestine_barcodes_binary[intestine_barcodes_binary$sum>=1,]
intestine_barcodes <- rownames(intestine_barcodes_binary)

barcodes=intestine_barcodes
barcodes <- barcodes[!duplicated(barcodes)]
df <- seurat[barcodes,]

pdf(paste0(directory, "/Intestine_coexpression.pdf"), height=8, width=9)
ggplot() + 
  geom_point(data=seurat, aes(x=UMAP_1, y=UMAP_2), color="gray80", size=3, shape=16) + 
  geom_point(data=df, aes(x=UMAP_1, y=UMAP_2), color="black", stroke=0.1, fill="red", shape=21, size=3) +
  theme_light(
    base_family = "sans",
    base_line_size = 0.2,
    base_rect_size = 0.2)
dev.off()





####### UMAP expression of epithelial genes
gene <- rev(c("CDH1", "EPCAM", "KRT8", "KRT10", "KRT17", "KRT19"))
gene <- gene[gene %in% rownames(seurat_object_epithelia@assays[["RNA"]]@data)]

clustering=lapply(1:length(gene), function(x){
  
  gene_data <- as.data.frame(seurat_object_epithelia@assays[["RNA"]]@data[gene[x],])
  colnames(gene_data) <- "data"
  gene_data$Barcode <- rownames(gene_data)
  gene_data <- merge(UMAPdata, gene_data, by="Barcode")
  rownames(gene_data) <- gene_data$Barcode
  Gene_pos <- (filter(gene_data, gene_data$data!=0))
  Gene_neg <- (filter(gene_data, gene_data$data==0))
  
  setwd(paste0(working_directory, "UMAP_gene_highlights"))
  tiff(paste0(gene[x], '_expression_scRNAseq.tiff'), units="in", width=5, height=4, res=150, compression = 'lzw')
  plot <- ggplot() + geom_point(data = Gene_neg, size = 0.4, shape = 16,  aes(x = UMAP_1, y = UMAP_2, col = data)) +
    geom_point(data = Gene_pos, size = 0.4, shape = 16,  aes(x = UMAP_1, y = UMAP_2, col = data)) +
    labs(x="UMAP1", y="UMAP2",
         title=(paste(gene[x], " Expression in scRNAseq"))) +
    scale_color_gradient2(mid = "gray90", high="#FF8C00", name=gene[x]) + theme(panel.background = element_rect(fill = 'white', colour = 'black')) +
    theme_light(base_family = "sans",base_line_size = 0.5,base_rect_size = 0.1,base_size = 14)
  print(plot)
  dev.off()
})


#### Epithelial marker violin
genes <- rev(c("EPCAM", "CDH1", "KRT8", "KRT10", "KRT17", "KRT19"))
genes_all <- genes
genes <- genes[genes %in% rownames(seurat_object_epithelia@assays[["RNA"]]@data)]
data <- as.data.frame(seurat_object_epithelia@assays[["RNA"]]@data[genes[1],])

for (x in c(2:(length(genes)))) {
  data[,x] <- seurat_object_epithelia@assays[["RNA"]]@data[genes[x],]
}
colnames(data) <- genes
data$Barcode <- rownames(data)

UMAPlong <- gather_(data, "Gene", "Expression", genes)
UMAPlong[UMAPlong$Expression==0,] <- NA
UMAPlong <- UMAPlong[complete.cases(UMAPlong), ]

data_summary <- function(x) {m <- mean(x)
ymin <- m-sd(x)
ymax <- m+sd(x)
return(c(y=m,ymin=ymin,ymax=ymax))
}

setwd(working_directory)
tiff("Expression_of_Epithelial_markers_in_all_of_the_AF_SC.tiff", units="in", width=5, height=9, res=150, compression = 'lzw')
ggplot(UMAPlong, aes(x=Gene, y=Expression)) +
  geom_violin(color="black") +
  stat_summary(fun.data=data_summary) +
  theme_light(base_family = "Arial",base_line_size = 0.5,base_rect_size = 0.1,base_size = 20) +
  theme(text=element_text(color="black", size=14, face="bold")) +
  scale_x_discrete(limits = genes_all) +
  ylim(c(0,ceiling(max(data[,c(1:length(genes))])))) +
  ylab("Epithelial Marker Expression") + xlab(NULL) + NoLegend() +
  ggtitle("All clusters") +
  coord_flip()
dev.off()


