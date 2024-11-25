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
library(fgsea)
library(ggridges)
library(scales)
library(ggrepel)
library(SingleR)

sessionInfo()

## Set paths
  base_directory = "XXXX"
  cellbender_directory = paste0(base_directory, "cellbender_outputs/")

working_directory = paste0(base_directory, "AF_scRNAseq_analysis/")

input_path = paste0(base_directory, "Input_Data/")

gc()
############################################## Import raw SC data
dir.create(working_directory)
dir.create(paste0(working_directory, "Quality_Control"))
lapply(c(1:11), function(i){
  if (i==1) raw_sample_name <- "687"
  if (i==2) raw_sample_name <- "689"
  if (i==3) raw_sample_name <- "696"
  if (i==4) raw_sample_name <- "Combo"
  if (i==5) raw_sample_name <- "713"
  if (i==6) raw_sample_name <- "719"
  if (i==7) raw_sample_name <- "725"
  if (i==8) raw_sample_name <- "741"
  if (i==9) raw_sample_name <- "758a"
  if (i==10) raw_sample_name <- "744"
  if (i==11) raw_sample_name <- "774"
  
  message(paste0(raw_sample_name, " beginning"))
  
    matrix <- Read10X_h5(paste0(cellbender_directory, raw_sample_name, "/CellBender_", raw_sample_name ,".h5"),
                         use.names = TRUE, unique.features = TRUE)
    cellbender_barcodes <- read.csv(paste0(cellbender_directory, raw_sample_name, "/CellBender_", raw_sample_name ,"_cell_barcodes.csv"), header=FALSE)
    cellbender_barcodes$V1 <- paste(raw_sample_name, cellbender_barcodes$V1, sep="-")
    all_barcodes <- colnames(matrix)
    all_barcodes <- paste(raw_sample_name, all_barcodes, sep="-")
    cellbender_empty_barcodes <- all_barcodes[!all_barcodes %in% cellbender_barcodes$V1]
    assign((paste0("cellbender_empty_barcodes", sep="_", raw_sample_name)), cellbender_empty_barcodes, envir=.GlobalEnv)

  seurat_object <- CreateSeuratObject(counts = matrix, project = "scRNAseq", min.cells=0, min.features = 1)
  
  seurat_object <- RenameCells(seurat_object, add.cell.id = raw_sample_name)
  seurat_object <- RenameCells(seurat_object, new.names = gsub("_", "-", colnames(seurat_object)))
  
  dir.create(paste0(working_directory, "Quality_Control", sep="/", raw_sample_name))
  setwd(paste0(working_directory, "Quality_Control", sep="/", raw_sample_name))
  
  seurat_object@meta.data[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "MT-")[,1]
  seurat_object@meta.data[["HSP"]] <- PercentageFeatureSet(seurat_object, pattern = "HSP")[,1]
  seurat_object@meta.data[["Ribo1"]] <- PercentageFeatureSet(seurat_object, pattern = "RPL")[,1]
  seurat_object@meta.data[["Ribo2"]] <- PercentageFeatureSet(seurat_object, pattern = "RPS")[,1]
  seurat_object@meta.data[["HB"]] <- PercentageFeatureSet(seurat_object, pattern = "HB")[,1]
  
  # Define genes that drive lots of technical variance & remove from sample
  hkGeneREGEX=c("DNAJ", "EIF", "RPL", "RPS", "RPN1", "POLR", "SNX", "HSP", "H1FX",
                "H2AF", "PRKA", "NDUF", "ATP", "PSM", "UBA", "UBE", "USP", "TXN")
  coreExcludeGenes = unlist(unique(c(grep('\\.[0-9]+',rownames(seurat_object),value=TRUE), #Poorly characterised
                                     grep('MALAT1',rownames(seurat_object),value=TRUE), #Contamination or highly expressed & poorly characterised
                                     grep('^MT-',rownames(seurat_object),value=TRUE), #Mitochondria
                                     lapply(1:length(hkGeneREGEX), function(e){
                                       grep(hkGeneREGEX[e], rownames(seurat_object), value = TRUE)
                                     })
  ))
  )
  filtered_genes <- c(rownames(seurat_object), coreExcludeGenes)
  filtered_genes <- filtered_genes[!(duplicated(filtered_genes) | duplicated(filtered_genes, fromLast = TRUE))]
  seurat_object <- subset(seurat_object, features = filtered_genes)
  
  # Extract into seperate matrix
  metadata <- seurat_object@meta.data
  metadata$Barcode <- rownames(metadata)
  
  # Feature QC Plots
  pdf("FeatureScatterplots_preFiltering.pdf")
  print(ggplot() + geom_point(data = metadata, aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt)))
  print(ggplot() + geom_point(data = metadata, aes(x = nCount_RNA, y = nFeature_RNA, color = souporcell_doublet_status)))
  dev.off()
  
  # Violin QC Plots
  pdf("ViolinPlots_preFiltering.pdf")
  print(VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 2))
  dev.off()
  
  pdf("log_nCounts_preFiltering.pdf")
  print(ggplot() + geom_histogram(data = metadata, aes(x = log(nCount_RNA))))
  dev.off()
  
  seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)

  assign((paste0("seurat_object", sep="_", raw_sample_name)), seurat_object, envir=.GlobalEnv)
  assign((paste0("metadata", sep="_", raw_sample_name)), seurat_object@meta.data, envir=.GlobalEnv)
  assign((paste0("seurat_object_barcodes", sep="_", raw_sample_name)), colnames(seurat_object), envir=.GlobalEnv)
  message(paste0(raw_sample_name, " complete"))
})

#### Split combo lane
newID <- read.table(file=paste0(base_directory, "R_scripts/Input_Data/donor_ids_ComboLane21.tsv"), header = T)
newID$cell <- paste("Combo", newID$cell, sep="-")
newID <- split(newID, f=newID$best_singlet)

barcodes680 <- newID$donor1
barcodes680 <- barcodes680$cell
barcodes695 <- newID$donor2
barcodes695 <- barcodes695$cell
barcodes696 <- newID$donor0
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
seurat_object_758a@meta.data[["orig_sample"]] <- "28PCW_758a"
seurat_object_774@meta.data[["orig_sample"]] <- "20PCW_774"

seurat_object_687@meta.data[["PCW"]] <- 32
seurat_object_689@meta.data[["PCW"]] <- 17
seurat_object_696@meta.data[["PCW"]] <- 30
seurat_object_680@meta.data[["PCW"]] <- 17
seurat_object_695@meta.data[["PCW"]] <- 13
seurat_object_713@meta.data[["PCW"]] <- 14
seurat_object_719@meta.data[["PCW"]] <- 17
seurat_object_725@meta.data[["PCW"]] <- 19
seurat_object_741@meta.data[["PCW"]] <- 28
seurat_object_758a@meta.data[["PCW"]] <- 28
seurat_object_774@meta.data[["PCW"]] <- 20

seurat_object_687@meta.data[["Trimester"]] <- "Third"
seurat_object_689@meta.data[["Trimester"]] <- "Second"
seurat_object_696@meta.data[["Trimester"]] <- "Third"
seurat_object_680@meta.data[["Trimester"]] <- "Second"
seurat_object_695@meta.data[["Trimester"]] <- "Second"
seurat_object_713@meta.data[["Trimester"]] <- "Second"
seurat_object_719@meta.data[["Trimester"]] <- "Second"
seurat_object_725@meta.data[["Trimester"]] <- "Second"
seurat_object_741@meta.data[["Trimester"]] <- "Third"
seurat_object_758a@meta.data[["Trimester"]] <- "Third"
seurat_object_774@meta.data[["Trimester"]] <- "Second"

#### Select appropriate 744 cells
newID <- read.table(file=paste0(base_directory, "R_scripts/Input_Data/donor_ids_744.tsv"), header = T)
newID$cell <- paste("744", newID$cell, sep="-")
newID <- split(newID, f=newID$donor_id)

barcodes744_744 <- newID$donor1
barcodes744_744 <- barcodes744_744$cell

seurat_object_744 <- subset(seurat_object_744, cells = barcodes744_744)

seurat_object_744@meta.data[["orig_sample"]] <- "25PCW_744"
seurat_object_744@meta.data[["PCW"]] <- 25
seurat_object_744@meta.data[["Trimester"]] <- "Third"

seurat_combined_list <- list()
seurat_combined_list[["Batch1"]] <- seurat_object_687
seurat_combined_list[["Batch2"]] <- seurat_object_689
seurat_combined_list[["Batch3"]] <- seurat_object_696
seurat_combined_list[["Batch4"]] <- seurat_object_680
seurat_combined_list[["Batch5"]] <- seurat_object_695
seurat_combined_list[["Batch6"]] <- seurat_object_713
seurat_combined_list[["Batch7"]] <- seurat_object_719
seurat_combined_list[["Batch8"]] <- seurat_object_725
seurat_combined_list[["Batch10"]] <- seurat_object_741
seurat_combined_list[["Batch12"]]<- seurat_object_744
seurat_combined_list[["Batch13"]] <- seurat_object_758a
seurat_combined_list[["Batch18"]] <- seurat_object_774

features <- SelectIntegrationFeatures(object.list = seurat_combined_list, nfeatures = 2000)
anchors <- FindIntegrationAnchors(object.list = seurat_combined_list, anchor.features = features)
AF_seurat_object <- IntegrateData(anchorset = anchors)
DefaultAssay(AF_seurat_object) <- "integrated"
num_cells_analysed=length(colnames(AF_seurat_object))

# Filter
  AF_seurat_object <- subset(AF_seurat_object, subset = nFeature_RNA > 150)


# Select epithelial cells
AF_seurat_object_processed <- AF_seurat_object
AF_seurat_object_processed <- FindVariableFeatures(object = AF_seurat_object_processed, selection.method="vst", nfeatures = num_DEGS_to_use)
AF_seurat_object_processed <- ScaleData(AF_seurat_object_processed, features = rownames(AF_seurat_object_processed), do.center = TRUE, do.scale = TRUE)
AF_seurat_object_processed <- RunPCA(AF_seurat_object_processed, features = VariableFeatures(object = AF_seurat_object_processed), npcs = 50)
number_of_PCAs = 20
AF_seurat_object_processed <- RunUMAP(AF_seurat_object_processed, dims = 1:number_of_PCAs)
AF_seurat_object_processed <- RunTSNE(object = AF_seurat_object_processed, check_duplicates = FALSE)
AF_seurat_object_processed <- FindNeighbors(AF_seurat_object_processed, dims = 1:number_of_PCAs)
AF_seurat_object_processed <- FindClusters(AF_seurat_object_processed, resolution = 0.2)

setwd(working_directory)
AF_UMAPdata=as.data.frame(AF_seurat_object_processed@meta.data)
AF_UMAPdata$UMAP_1 <- AF_seurat_object_processed@reductions[["umap"]]@cell.embeddings[,"UMAP_1"]
AF_UMAPdata$UMAP_2 <- AF_seurat_object_processed@reductions[["umap"]]@cell.embeddings[,"UMAP_2"]
AF_UMAPdata$Barcode <- rownames(AF_UMAPdata)

## Remove low quality cluster
low_quality_clusters=XX

Idents(AF_seurat_object_processed) <- AF_seurat_object_processed$seurat_clusters
AF_seurat_object_processed <- subset(x = AF_seurat_object_processed, idents = low_quality_clusters, invert = TRUE)

AF_seurat_object_processed <- FindVariableFeatures(object = AF_seurat_object_processed, selection.method="vst", nfeatures = num_DEGS_to_use)
AF_seurat_object_processed <- ScaleData(AF_seurat_object_processed, features = rownames(AF_seurat_object_processed), do.center = TRUE, do.scale = TRUE)
AF_seurat_object_processed <- RunPCA(AF_seurat_object_processed, features = VariableFeatures(object = AF_seurat_object_processed), npcs = 50)
number_of_PCAs = 20
AF_seurat_object_processed <- RunUMAP(AF_seurat_object_processed, dims = 1:number_of_PCAs)
AF_seurat_object_processed <- RunTSNE(object = AF_seurat_object_processed, check_duplicates = FALSE)
AF_seurat_object_processed <- FindNeighbors(AF_seurat_object_processed, dims = 1:number_of_PCAs)
AF_seurat_object_processed <- FindClusters(AF_seurat_object_processed, resolution = 0.2)


#### Plot expression of specific epithelial genes
gene <- rev(c("CDH1", "EPCAM", "KRT8", "KRT10", "KRT17", "KRT13"))

gene <- gene[gene %in% rownames(AF_seurat_object_processed@assays[["RNA"]]@data)]

dir.create(paste0(working_directory, "UMAP_gene_highlights_allCells"))

AF_UMAPdata=as.data.frame(AF_seurat_object_processed@meta.data)
AF_UMAPdata$UMAP_1 <- AF_seurat_object_processed@reductions[["umap"]]@cell.embeddings[,"UMAP_1"]
AF_UMAPdata$UMAP_2 <- AF_seurat_object_processed@reductions[["umap"]]@cell.embeddings[,"UMAP_2"]
AF_UMAPdata$Barcode <- rownames(AF_UMAPdata)

clustering=lapply(1:length(gene), function(x){
  
  gene_data <- as.data.frame(AF_seurat_object_processed@assays[["RNA"]]@data[gene[x],])
  
  colnames(gene_data) <- "data"
  gene_data$Barcode <- rownames(gene_data)
  gene_data <- merge(AF_UMAPdata, gene_data, by="Barcode")
  rownames(gene_data) <- gene_data$Barcode
  Gene_pos <- (filter(gene_data, gene_data$data>=3))
  Gene_neg <- (filter(gene_data, gene_data$data<3))
  
  setwd(paste0(working_directory, "UMAP_gene_highlights_allCells"))
  tiff(paste0(gene[x], '_expression_scRNAseq.tiff'), units="in", width=5, height=4, res=150, compression = 'lzw')
  plot <- ggplot() +
    geom_point(data = gene_data, size = 0.4, shape = 16,  aes(x = UMAP_1, y = UMAP_2), col = "gray90") +
    geom_point(data = Gene_neg[Gene_neg$data>2,], size = 0.4, shape = 16,  aes(x = UMAP_1, y = UMAP_2, col = data)) +
    geom_point(data = Gene_pos[Gene_pos$data>2,], size = 0.4, shape = 16,  aes(x = UMAP_1, y = UMAP_2, col = data)) +
    labs(x="UMAP1", y="UMAP2",
         title=(paste(gene[x], " Expression in scRNAseq"))) +
    scale_color_gradient2(mid = "gray90", high="#FF8C00", name=gene[x]) + theme(panel.background = element_rect(fill = 'white', colour = 'black')) +
    theme_light(base_family = "sans",base_line_size = 0.5,base_rect_size = 0.1,base_size = 14) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(plot)
  dev.off()
})




### Label epithelial cells using SingleR
  SCE <- Seurat::as.SingleCellExperiment(AF_seurat_object_processed)
  
  ref <- celldex::HumanPrimaryCellAtlasData()
  pred_HCA <- SingleR::SingleR(test = SCE, ref = ref, assay.type.test=1,
                               labels = ref$label.main)
  pred_HCA <- as.data.frame(pred_HCA)
  pred_HCA$Barcode <- rownames(pred_HCA)
  
  setwd(working_directory)
  write.csv(pred_HCA, "singleR_output_HCA.csv")

setwd(working_directory)
pred_HCA <- read.csv("singleR_output_HCA.csv")
colnames(pred_HCA)[colnames(pred_HCA) == 'X'] <- 'Barcode'

pred_HCA <- pred_HCA[,c("Barcode", "pruned.labels")]

AF_UMAPdata=as.data.frame(AF_seurat_object_processed@meta.data)
AF_UMAPdata$UMAP_1 <- AF_seurat_object_processed@reductions[["umap"]]@cell.embeddings[,"UMAP_1"]
AF_UMAPdata$UMAP_2 <- AF_seurat_object_processed@reductions[["umap"]]@cell.embeddings[,"UMAP_2"]
AF_UMAPdata$Barcode <- rownames(AF_UMAPdata)

pred_HCA <- pred_HCA[,c("Barcode", "pruned.labels")]
pred_HCA <- merge(pred_HCA, AF_UMAPdata, by="Barcode")
pred_HCA <- pred_HCA[!is.na(pred_HCA$pruned.labels),]

## Extract colours from the viridis rainbow colours scale package
colour_plot <- ggplot(pred_HCA, aes(x = UMAP_1, y = UMAP_2, color=(factor(pruned.labels)))) +
  scale_color_viridis(option="turbo", discrete=TRUE)
my_colours <- unique((ggplot_build(colour_plot))[["data"]][[1]][["colour"]])

pred_HCA$epi_label <- "Unlabelled"
pred_HCA[pred_HCA$pruned.labels=="Epithelial_cells","epi_label"] <- "Epithelial_cells"

plot <- ggplot(pred_HCA, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(
    aes(color=epi_label),
    size=1.05,
    alpha=0.95,
    shape=16,             
    stroke=0) +
  guides(
    colour = guide_legend(override.aes = list(size=5, stroke=0)))  + 
  scale_color_manual(
    values = c("darkorange2", "gray80"),
    name="Cell Labels") +
  labs(x="UMAP 1", y="UMAP 2", title="Cell Labelling with SingleR",
       caption = paste0("Made by M. Beesley on ", Sys.Date())) +
  theme_light(base_family = "sans",base_line_size = 0.2,base_rect_size = 0.2) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

pdf(paste0(working_directory, "SingleR_cellLabelling_epithelialOnly.pdf"), height=8, width=9)
print(plot)
dev.off()

AF_seurat_object_processed_allCells <- AF_seurat_object_processed



AF_UMAPdata=as.data.frame(AF_seurat_object_processed@meta.data)
AF_UMAPdata$UMAP_1 <- AF_seurat_object_processed@reductions[["umap"]]@cell.embeddings[,"UMAP_1"]
AF_UMAPdata$UMAP_2 <- AF_seurat_object_processed@reductions[["umap"]]@cell.embeddings[,"UMAP_2"]
AF_UMAPdata$Barcode <- rownames(AF_UMAPdata)
AF_UMAPdata <- merge(AF_UMAPdata, B4_metadata, by="Barcode", all.x=TRUE)
AF_UMAPdata_backup <- AF_UMAPdata
write.csv(AF_UMAPdata, paste0(working_directory, "AF_UMAPdata_V1.csv"))

AF_UMAPdata_allCells <- AF_UMAPdata

setwd(working_directory)
pdf("Batch_Differences_UMAP_QC.pdf", height=6, width=8)
colour_plot <- ggplot(AF_UMAPdata, aes(x = UMAP_1, y = UMAP_2, color=orig_sample)) +
  scale_color_viridis(option="turbo", discrete=TRUE)
my_colours <- unique((ggplot_build(colour_plot))[["data"]][[1]][["colour"]])
print(ggplot() + geom_point(data = AF_UMAPdata, size = 0.4, shape = 16,
                            aes(x = UMAP_1, y = UMAP_2, color = orig_sample)) +
        scale_color_manual(values=my_colours, name="Sample") +
        guides(colour = guide_legend(override.aes = list(size=5))) +
        theme_light(base_family = "sans",base_line_size = 0.5,base_rect_size = 0.1,base_size = 14) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
dev.off()

setwd(working_directory)
pdf("Gestational_age_split_AF.pdf", width=9, height=7)
print(ggplot() + geom_point(data = AF_UMAPdata, size = 0.4, shape = 16,
                            aes(x = UMAP_1, y = UMAP_2, color = Trimester)) +
        guides(colour = guide_legend(override.aes = list(size=5))) +
        theme_light(base_family = "sans",base_line_size = 0.5,base_rect_size = 0.1,base_size = 14) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

print(ggplot() + geom_point(data = AF_UMAPdata, size = 0.4, shape = 16,
                            aes(x = UMAP_1, y = UMAP_2, color = PCW)) +
        scale_color_continuous(breaks=c(13,20,25,32)) +
        theme_light(base_family = "sans",base_line_size = 0.5,base_rect_size = 0.1,base_size = 14) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
dev.off()

#### UMAP QC Plots
setwd(paste0(working_directory, "Quality_Control"))
pdf("UMAP_QC.pdf", width=9, height=7)
ggplot() + geom_point(data = AF_UMAPdata, size = 0.4, shape = 16,
                      aes(x = UMAP_1, y = UMAP_2, color = as.numeric(percent.mt))) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme_light(base_family = "sans",base_line_size = 0.5,base_rect_size = 0.1,base_size = 14)

ggplot() + geom_point(data = AF_UMAPdata, size = 0.4, shape = 16,
                      aes(x = UMAP_1, y = UMAP_2, color = as.numeric(nFeature_RNA))) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme_light(base_family = "sans",base_line_size = 0.5,base_rect_size = 0.1,base_size = 14)

ggplot() + geom_point(data = AF_UMAPdata, size = 0.4, shape = 16,
                      aes(x = UMAP_1, y = UMAP_2, color = as.numeric(nCount_RNA))) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme_light(base_family = "sans",base_line_size = 0.5,base_rect_size = 0.1,base_size = 14)
dev.off()



## Select epithelial clusters
epithelial_clusters=c(x,x,x,x)
resolution=xx

Idents(AF_seurat_object_processed) <- AF_seurat_object_processed$seurat_clusters
AF_seurat_object_processed <- subset(x = AF_seurat_object_processed, idents = epithelial_clusters, invert = FALSE)
AF_seurat_object_processed <- FindVariableFeatures(object = AF_seurat_object_processed, selection.method="vst", nfeatures = num_DEGS_to_use)
AF_seurat_object_processed <- ScaleData(AF_seurat_object_processed, features = rownames(AF_seurat_object_processed), do.center = TRUE, do.scale = TRUE)
AF_seurat_object_processed <- RunPCA(AF_seurat_object_processed, features = VariableFeatures(object = AF_seurat_object_processed), npcs = 50)
number_of_PCAs = 20
AF_seurat_object_processed <- RunUMAP(AF_seurat_object_processed, dims = 1:number_of_PCAs)
AF_seurat_object_processed <- RunTSNE(object = AF_seurat_object_processed, check_duplicates = FALSE)
AF_seurat_object_processed <- FindNeighbors(AF_seurat_object_processed, dims = 1:number_of_PCAs)
AF_seurat_object_processed <- FindClusters(AF_seurat_object_processed, resolution = resolution)

AF_seurat_object_processed_preAnalysis <- AF_seurat_object_processed





#### Analysis of epithelial
DefaultAssay(AF_seurat_object_processed) <- "RNA"

AF_UMAPdata=as.data.frame(AF_seurat_object_processed@meta.data)
AF_UMAPdata$UMAP_1 <- AF_seurat_object_processed@reductions[["umap"]]@cell.embeddings[,"UMAP_1"]
AF_UMAPdata$UMAP_2 <- AF_seurat_object_processed@reductions[["umap"]]@cell.embeddings[,"UMAP_2"]
AF_UMAPdata$seurat_clusters=AF_seurat_object_processed@meta.data$seurat_clusters
AF_UMAPdata$Barcode <- rownames(AF_UMAPdata)
AF_UMAPdata <- merge(AF_UMAPdata, B4_metadata, by="Barcode", all.x=TRUE)
AF_UMAPdata_backup <- AF_UMAPdata


## Extract colours from the viridis rainbow colours scale package
colour_plot <- ggplot(AF_UMAPdata, aes(x = UMAP_1, y = UMAP_2, color=(factor(seurat_clusters)))) +
  scale_color_viridis(option="turbo", discrete=TRUE)
my_colours <- unique((ggplot_build(colour_plot))[["data"]][[1]][["colour"]])
my_colours <- colorspace::desaturate(my_colours, 0.01)

plot <- ggplot(AF_UMAPdata, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(
    aes(color=as.factor(seurat_clusters)),
    size=1.05,
    alpha=0.95,
    shape=16,
    stroke=0) +
  guides(
    colour = guide_legend(override.aes = list(size=5, stroke=0)))  +
  scale_color_manual(
    values = my_colours,
    name="Clusters:") +
  labs(
    x="UMAP 1",
    y="UMAP 2",
    title="Standard Clustering",
    caption = paste0("Made by M. Beesley on ", Sys.Date())
  ) +
  theme_light(
    base_family = "sans",
    base_line_size = 0.2,
    base_rect_size = 0.2
  )
pdf(paste0(working_directory, "/", "Clustering", "_", ceiling(num_cells_analysed/1000), "k.pdf"), height=8, width=9)
print(plot)
dev.off()





top_pathways_to_examine=1
probability_cutoff=40

pathways <- escape::getGeneSets(library = "C8")
pathway_names <- as.data.frame((names(pathways)))
pathway_names <- pathway_names[,1]

# Filter for relevant tissues
pathway_names_filtered <- c(pathway_names[grepl('KIDNEY', pathway_names)], pathway_names[grepl('LUNG', pathway_names)],
                            pathway_names[grepl('ESOPHAGE', pathway_names)], pathway_names[grepl('STOMACH', pathway_names)],
                            pathway_names[grepl('OLFACTORY', pathway_names)], pathway_names[grepl('GASTRIC', pathway_names)],
                            pathway_names[grepl('EPITHEL', pathway_names)], pathway_names[grepl('LIVER', pathway_names)],
                            pathway_names[grepl('INTESTIN', pathway_names)])

# Filter for fetal and progenitor
pathway_names_filtered <- c(pathway_names_filtered[grepl('FET', pathway_names_filtered)],
                            pathway_names_filtered[grepl('PROGENITOR', pathway_names_filtered)],
                            pathway_names_filtered[grepl('STEM', pathway_names_filtered)])
# Filter out immune
pathway_names_filtered <- pathway_names_filtered[!grepl('IMMUNE', pathway_names_filtered)]
pathway_names_filtered <- pathway_names_filtered[!grepl('BLOOD', pathway_names_filtered)]
pathway_names_filtered <- pathway_names_filtered[!grepl('LYMP', pathway_names_filtered)]
pathway_names_filtered <- pathway_names_filtered[!grepl('MEYL', pathway_names_filtered)]
pathway_names_filtered <- pathway_names_filtered[!grepl('KARY', pathway_names_filtered)]
pathway_names_filtered <- pathway_names_filtered[!grepl('ERYT', pathway_names_filtered)]
pathway_names_filtered <- pathway_names_filtered[!grepl('ENS', pathway_names_filtered)]

pathway_names_filtered <- pathway_names_filtered[!duplicated(pathway_names_filtered)]
refined_pathways <- pathways[pathway_names_filtered]


GSEA_C8_label <- escape::enrichIt(obj = AF_seurat_object_processed,
                                  gene.sets = refined_pathways,
                                  groups = 1000, cores = 12)

GSEA_C8_label <- as.data.frame(GSEA_C8_label)
GSEA_C8_label <- t(GSEA_C8_label)

pathways_tissue <- as.data.frame(rownames(GSEA_C8_label))
colnames(pathways_tissue) <- "Pathway"

pathways_tissue[grepl("EMBRYO", pathways_tissue[,1], fixed = TRUE), "GO_Tissue"] <- "Embryo"
pathways_tissue[grepl("FETAL", pathways_tissue[,1], fixed = TRUE), "GO_Tissue"] <- "Fetal"
pathways_tissue[grepl("EPITHEL", pathways_tissue[,1], fixed = TRUE), "GO_Tissue"] <- "Epithelial"
pathways_tissue[grepl("SQUAMOUS", pathways_tissue[,1], fixed = TRUE), "GO_Tissue"] <- "Squamous Epithelial"
pathways_tissue[grepl("CILIATED", pathways_tissue[,1], fixed = TRUE), "GO_Tissue"] <- "Ciliated Epithelial"
pathways_tissue[grepl("PROGENITOR", pathways_tissue[,1], fixed = TRUE), "GO_Tissue"] <- "Progenitor"
pathways_tissue[grepl("CORNEAL", pathways_tissue[,1], fixed = TRUE), "GO_Tissue"] <- "Cornea"
pathways_tissue[grepl("RETINA", pathways_tissue[,1], fixed = TRUE), "GO_Tissue"] <- "Retina"
pathways_tissue[grepl("EYE", pathways_tissue[,1], fixed = TRUE), "GO_Tissue"] <- "Eye"
pathways_tissue[grepl("OVARY", pathways_tissue[,1], fixed = TRUE), "GO_Tissue"] <- "Ovary"
pathways_tissue[grepl("LIVER", pathways_tissue[,1], fixed = TRUE), "GO_Tissue"] <- "Liver"
pathways_tissue[grepl("OLFACTORY", pathways_tissue[,1], fixed = TRUE), "GO_Tissue"] <- "Olfactory"
pathways_tissue[grepl("KIDNEY", pathways_tissue[,1], fixed = TRUE), "GO_Tissue"] <- "Kidney"
pathways_tissue[grepl("BLADD", pathways_tissue[,1], fixed = TRUE), "GO_Tissue"] <- "Bladder"
pathways_tissue[grepl("URET", pathways_tissue[,1], fixed = TRUE), "GO_Tissue"] <- "Ureter"
pathways_tissue[grepl("GASTRIC", pathways_tissue[,1], fixed = TRUE), "GO_Tissue"] <- "GI"
pathways_tissue[grepl("DUODEN", pathways_tissue[,1], fixed = TRUE), "GO_Tissue"] <- "Duodenum"
pathways_tissue[grepl("INTESTIN", pathways_tissue[,1], fixed = TRUE), "GO_Tissue"] <- "GI"
pathways_tissue[grepl("STOMACH", pathways_tissue[,1], fixed = TRUE), "GO_Tissue"] <- "GI"
pathways_tissue[grepl("ESOPHAGE", pathways_tissue[,1], fixed = TRUE), "GO_Tissue"] <- "Esophageous"
pathways_tissue[grepl("LUNG", pathways_tissue[,1], fixed = TRUE), "GO_Tissue"] <- "Lung"
pathways_tissue[grepl("BRONCHIOLAR", pathways_tissue[,1], fixed = TRUE), "GO_Tissue"] <- "Lung"
pathways_tissue[grepl("RESPIRA", pathways_tissue[,1], fixed = TRUE), "GO_Tissue"] <- "Lung"
pathways_tissue[grepl("PLACENTA", pathways_tissue[,1], fixed = TRUE), "GO_Tissue"] <- "Placenta"
pathways_tissue[grepl("IMMUNE", pathways_tissue[,1], fixed = TRUE), "GO_Tissue"] <- "Immune"
pathways_tissue[grepl("NEURON", pathways_tissue[,1], fixed = TRUE), "GO_Tissue"] <- "Immune"
pathways_tissue[grepl("LYMPHOCYTE", pathways_tissue[,1], fixed = TRUE), "GO_Tissue"] <- "Immune"
pathways_tissue[grepl("MACROPHAGE", pathways_tissue[,1], fixed = TRUE), "GO_Tissue"] <- "Immune"

maximum <- data.frame()
lapply(c(1:nrow(pathways_tissue)), function(g){
  message(g)
  invisible(maximum[g,"maximum"] <- max(GSEA_C8_label[g,]))
  invisible(assign("maximum", maximum, envir=.GlobalEnv))
})

maximum$Pathway <- pathways_tissue$Pathway
maximum$GO_Tissue <- pathways_tissue$GO_Tissue


cell_database <- as.data.frame(GSEA_C8_label[,1])
cell_database$Pathway <- pathways_tissue$Pathway
cell_database$GO_Tissue <- pathways_tissue$GO_Tissue
colnames(cell_database) <- c("Probability", "Pathway", "GO_Tissue")
cell_database <- cell_database[order(cell_database$Probability, decreasing=TRUE),]
cell_database <- head(cell_database, top_pathways_to_examine)
data <- as.data.frame(table(cell_database[,"GO_Tissue"]))
data <- data[order(data$Freq, decreasing=TRUE),]
cell_database_output <- as.data.frame(data[1,1])
cell_database_output$Probability <- head(cell_database$Probability, 1)
cell_database_output$Barcode <- colnames(GSEA_C8_label)[1]
colnames(cell_database_output) <- c("GO_Tissue","Probability", "Barcode")

lapply(c(2:length(colnames(GSEA_C8_label))), function(g){
  message(g)
  invisible(cell <- as.data.frame(GSEA_C8_label[,g]))
  invisible(cell$Pathway <- pathways_tissue$Pathway)
  invisible(cell$GO_Tissue <- pathways_tissue$GO_Tissue)
  invisible(colnames(cell) <- c("Probability", "Pathway", "GO_Tissue"))
  invisible(cell <- cell[order(cell$Probability, decreasing=TRUE),])
  invisible(cell <- head(cell, top_pathways_to_examine))
  invisible(data <- as.data.frame(table(cell[,"GO_Tissue"])))
  invisible(data <- data[order(data$Freq, decreasing=TRUE),])
  invisible(cell_output <- as.data.frame(data[1,1]))
  invisible(cell_output$Probability <- head(cell$Probability, 1))
  invisible(colnames(cell_output) <- c("GO_Tissue","Probability"))
  invisible(cell_output$Barcode <- colnames(GSEA_C8_label)[g])
  invisible(cell_database_output <- rbind(cell_database_output,cell_output))
  invisible(assign("cell_database_output", cell_database_output, envir=.GlobalEnv))
})

cell_database_output$GO_Tissue <- as.character(cell_database_output$GO_Tissue)
cell_database_output[cell_database_output$Probability<(probability_cutoff/100), "GO_Tissue"] <- "Unknown"
cell_database_output[is.na(cell_database_output$GO_Tissue), "GO_Tissue"] <- "Unknown"

cell_database_output$Barcode <- gsub("X", "", cell_database_output$Barcode)
cell_database_output$Barcode <- gsub("\\.", "-", cell_database_output$Barcode)

cell_database_output <- cell_database_output[!is.na(cell_database_output$Barcode),]
AF_UMAPdata_GO <- AF_UMAPdata[,c("Barcode", "seurat_clusters", "UMAP_1", "UMAP_2", "orig_sample", "Trimester", "PCW", "nCount_RNA", "nFeature_RNA")]
AF_UMAPdata_GO = merge(AF_UMAPdata_GO, cell_database_output, by="Barcode", all.x=TRUE)

  ggplot() +
    geom_point(data=AF_UMAPdata_backup, aes(x = UMAP_1, y = UMAP_2),
               color="gray80", size=1.05, alpha=1, shape=16, stroke=0) +
    geom_point(data=AF_UMAPdata_GO[AF_UMAPdata_GO$GO_Tissue=="Kidney",], aes(x = UMAP_1, y = UMAP_2),
               color="darkorange4", size=1.05, alpha=1, shape=16, stroke=0) +
    labs(x="UMAP 1", y="UMAP 2", title=paste0("AF epithelia GO labelled, ", table(AF_UMAPdata_GO$GO_Tissue=="Kidney")[["TRUE"]],  " ", "Kidney", " cells")) +
    theme_light(base_family = "sans", base_line_size = 0.2, base_rect_size = 0.2) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  ggplot() +
    geom_point(data=AF_UMAPdata_backup, aes(x = UMAP_1, y = UMAP_2),
               color="gray80", size=1.05, alpha=1, shape=16, stroke=0) +
    geom_point(data=AF_UMAPdata_GO[AF_UMAPdata_GO$GO_Tissue=="Lung",], aes(x = UMAP_1, y = UMAP_2),
               color="darkorange4", size=1.05, alpha=1, shape=16, stroke=0) +
    labs(x="UMAP 1", y="UMAP 2", title=paste0("AF epithelia GO labelled, ", table(AF_UMAPdata_GO$GO_Tissue=="Lung")[["TRUE"]],  " ", "Lung", " cells")) +
    theme_light(base_family = "sans", base_line_size = 0.2, base_rect_size = 0.2) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  ggplot() +
    geom_point(data=AF_UMAPdata_backup, aes(x = UMAP_1, y = UMAP_2),
               color="gray80", size=1.05, alpha=1, shape=16, stroke=0) +
    geom_point(data=AF_UMAPdata_GO[AF_UMAPdata_GO$GO_Tissue=="GI",], aes(x = UMAP_1, y = UMAP_2),
               color="darkorange4", size=1.05, alpha=1, shape=16, stroke=0) +
    labs(x="UMAP 1", y="UMAP 2", title=paste0("AF epithelia GO labelled, ", table(AF_UMAPdata_GO$GO_Tissue=="GI")[["TRUE"]],  " ", "GI", " cells")) +
    theme_light(base_family = "sans", base_line_size = 0.2, base_rect_size = 0.2) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())






epithelial_Kidney_barcodes = cell_database_output[cell_database_output$GO_Tissue=="Kidney","Barcode"]
epithelial_lung_barcodes = cell_database_output[cell_database_output$GO_Tissue=="Lung","Barcode"]
epithelial_GI_barcodes = c(cell_database_output[cell_database_output$GO_Tissue=="GI","Barcode"])

#### Lung
epithelial_lung_barcodes <- epithelial_lung_barcodes[!duplicated(epithelial_lung_barcodes)]
seurat_object_epithelia <- subset(AF_seurat_object_processed, cells=epithelial_lung_barcodes)

genes <- rev(c("NKX2-1", "SOX9", "ETV4", "ETV5", "GATA6", "ID2"))
genes_all <- genes
genes <- genes[genes %in% rownames(seurat_object_epithelia@assays[["RNA"]]@data)]

## Reprocess just Lung cells
Idents(seurat_object_epithelia) <- seurat_object_epithelia$seurat_clusters
seurat_object_epithelia <- FindVariableFeatures(object = seurat_object_epithelia, selection.method="vst", nfeatures = num_DEGS_to_use)
seurat_object_epithelia <- ScaleData(seurat_object_epithelia, features = rownames(seurat_object_epithelia), do.center = TRUE, do.scale = TRUE)
seurat_object_epithelia <- RunPCA(seurat_object_epithelia, features = VariableFeatures(object = seurat_object_epithelia), npcs = 50)
number_of_PCAs = 20
seurat_object_epithelia <- RunUMAP(seurat_object_epithelia, dims = 1:number_of_PCAs)
seurat_object_epithelia <- RunTSNE(object = seurat_object_epithelia, check_duplicates = FALSE)
seurat_object_epithelia <- FindNeighbors(seurat_object_epithelia, dims = 1:number_of_PCAs)
seurat_object_epithelia <- FindClusters(seurat_object_epithelia, resolution = 0.99)

AF_UMAPdata=as.data.frame(seurat_object_epithelia@meta.data)
AF_UMAPdata$UMAP_1 <- seurat_object_epithelia@reductions[["umap"]]@cell.embeddings[,"UMAP_1"]
AF_UMAPdata$UMAP_2 <- seurat_object_epithelia@reductions[["umap"]]@cell.embeddings[,"UMAP_2"]
AF_UMAPdata$Barcode <- rownames(AF_UMAPdata)

seurat_object_epithelia = AddModuleScore(seurat_object_epithelia, features=list(genes), name="lung")
lung_data <- as.data.frame(seurat_object_epithelia$lung1)
lung_data$Barcode <- rownames(lung_data)
colnames(lung_data) <- c("Lung", "Barcode")
AF_UMAPdata_MS <- merge(AF_UMAPdata, lung_data, by="Barcode")

tiff("UMAP of Lung progenitor gene scores in the Lung AF cells.tiff")
print(ggplot() +
        geom_point(data = AF_UMAPdata_MS, size = 0.4, shape = 16, color="gray80",
                   aes(x = UMAP_1, y = UMAP_2)) +
        geom_point(data = AF_UMAPdata_MS[AF_UMAPdata_MS$Lung>0,], size = 0.4, shape = 16,
                   aes(x = UMAP_1, y = UMAP_2), color = "darkred") +
        guides(colour = guide_legend(override.aes = list(size=5, stroke=0))) +
        theme_light(base_family = "sans", base_line_size = 0.5, base_rect_size = 0.1, base_size = 14) +
        ggtitle(paste0("ModuleScore of Lung markers")) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
dev.off()

tiff("Violin plot of Lung progenitor gene scores in the Lung AF cells.tiff")
print(ggplot() + geom_violin(data=AF_UMAPdata_MS_GO, aes(x="Lung", y=Lung)) +
        geom_jitter(data=AF_UMAPdata_MS_GO, aes(x="Lung", y=Lung), size=0.3) +
        theme_light(base_family = "sans", base_line_size = 0.5, base_rect_size = 0.1, base_size = 14) +
        ggtitle("Module score of all lung cells for lung progenitor markers") +
        xlab(NULL) + ylab("MS") +
        geom_hline(yintercept=0.01, linetype="dashed", color="red", linewidth=1))
dev.off()





#### Kidney
epithelial_Kidney_barcodes <- epithelial_Kidney_barcodes[!duplicated(epithelial_Kidney_barcodes)]
seurat_object_epithelia <- subset(AF_seurat_object_processed, cells=epithelial_Kidney_barcodes)

seurat_object_epithelia_Kidney <- seurat_object_epithelia

seurat <- AF_seurat_object_processed
genes <- rev(c("PAX2", "PAX8", "LHX1", "JAG1", "SIX2", "RET", "HNF4A", "GATA3", "POU3F3",
               "WT1"))
genes_all <- genes
genes <- genes[genes %in% rownames(seurat_object_epithelia@assays[["RNA"]]@data)]

## Reprocess just Kidney cells
Idents(seurat_object_epithelia) <- seurat_object_epithelia$seurat_clusters
seurat_object_epithelia <- FindVariableFeatures(object = seurat_object_epithelia, selection.method="vst", nfeatures = num_DEGS_to_use)
seurat_object_epithelia <- ScaleData(seurat_object_epithelia, features = rownames(seurat_object_epithelia), do.center = TRUE, do.scale = TRUE)
seurat_object_epithelia <- RunPCA(seurat_object_epithelia, features = VariableFeatures(object = seurat_object_epithelia), npcs = 50)
number_of_PCAs = 20
seurat_object_epithelia <- RunUMAP(seurat_object_epithelia, dims = 1:number_of_PCAs)
seurat_object_epithelia <- RunTSNE(object = seurat_object_epithelia, check_duplicates = FALSE)
seurat_object_epithelia <- FindNeighbors(seurat_object_epithelia, dims = 1:number_of_PCAs)
seurat_object_epithelia <- FindClusters(seurat_object_epithelia, resolution = 0.99)

AF_UMAPdata=as.data.frame(seurat_object_epithelia@meta.data)
AF_UMAPdata$UMAP_1 <- seurat_object_epithelia@reductions[["umap"]]@cell.embeddings[,"UMAP_1"]
AF_UMAPdata$UMAP_2 <- seurat_object_epithelia@reductions[["umap"]]@cell.embeddings[,"UMAP_2"]
AF_UMAPdata$Barcode <- rownames(AF_UMAPdata)

seurat_object_epithelia = AddModuleScore(seurat_object_epithelia, features=list(genes), name="Kidney")
Kidney_data <- as.data.frame(seurat_object_epithelia$Kidney1)
Kidney_data$Barcode <- rownames(Kidney_data)
colnames(Kidney_data) <- c("Kidney", "Barcode")
AF_UMAPdata_MS <- merge(AF_UMAPdata, Kidney_data, by="Barcode")

pdf("MS_of_Kidney_markers_in_Kidney_AF_SC.pdf", height=5, width=4)
print(ggplot() +
        geom_point(data = AF_UMAPdata_MS, size = 0.4, shape = 16, color="gray80",
                   aes(x = UMAP_1, y = UMAP_2)) +
        geom_point(data = AF_UMAPdata_MS[AF_UMAPdata_MS$Kidney>0,], size = 0.4, shape = 16,
                   aes(x = UMAP_1, y = UMAP_2), color = "darkred") +
        guides(colour = guide_legend(override.aes = list(size=5, stroke=0))) +
        theme_light(base_family = "sans", base_line_size = 0.5, base_rect_size = 0.1, base_size = 14) +
        ggtitle(paste0("ModuleScore of Kidney markers")) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
dev.off()



#### GI
epithelial_GI_barcodes <- epithelial_GI_barcodes[!duplicated(epithelial_GI_barcodes)]
seurat_object_epithelia <- subset(AF_seurat_object_processed, cells=epithelial_GI_barcodes)

seurat <- seurat_object_epithelia
genes <- (c("LGR5", "OLFM4", "LRIG1", "CDX2", "CD44", "LYZ", "SMOC2", "PROCR"))
genes_all <- genes
genes <- genes[genes %in% rownames(seurat@assays[["RNA"]]@data)]

## Reprocess just GI cells
Idents(seurat_object_epithelia) <- seurat_object_epithelia$seurat_clusters
seurat_object_epithelia <- FindVariableFeatures(object = seurat_object_epithelia, selection.method="vst", nfeatures = num_DEGS_to_use)
seurat_object_epithelia <- ScaleData(seurat_object_epithelia, features = rownames(seurat_object_epithelia), do.center = TRUE, do.scale = TRUE)
seurat_object_epithelia <- RunPCA(seurat_object_epithelia, features = VariableFeatures(object = seurat_object_epithelia), npcs = 50)
number_of_PCAs = 20
seurat_object_epithelia <- RunUMAP(seurat_object_epithelia, dims = 1:number_of_PCAs)
seurat_object_epithelia <- RunTSNE(object = seurat_object_epithelia, check_duplicates = FALSE)
seurat_object_epithelia <- FindNeighbors(seurat_object_epithelia, dims = 1:number_of_PCAs)
seurat_object_epithelia <- FindClusters(seurat_object_epithelia, resolution = 0.99)

AF_UMAPdata=as.data.frame(seurat_object_epithelia@meta.data)
AF_UMAPdata$UMAP_1 <- seurat_object_epithelia@reductions[["umap"]]@cell.embeddings[,"UMAP_1"]
AF_UMAPdata$UMAP_2 <- seurat_object_epithelia@reductions[["umap"]]@cell.embeddings[,"UMAP_2"]
AF_UMAPdata$Barcode <- rownames(AF_UMAPdata)

seurat_object_epithelia = AddModuleScore(seurat_object_epithelia, features=list(genes), name="GI")
GI_data <- as.data.frame(seurat_object_epithelia$GI1)
GI_data$Barcode <- rownames(GI_data)
colnames(GI_data) <- c("GI", "Barcode")
AF_UMAPdata_MS <- merge(AF_UMAPdata, GI_data, by="Barcode")

pdf("MS_of_GI_markers_in_GI_AF_SC.pdf", height=5, width=4)
print(ggplot() +
        geom_point(data = AF_UMAPdata_MS, size = 0.4, shape = 16, color="gray80",
                   aes(x = UMAP_1, y = UMAP_2)) +
        geom_point(data = AF_UMAPdata_MS[AF_UMAPdata_MS$GI>0,], size = 0.4, shape = 16,
                   aes(x = UMAP_1, y = UMAP_2), color = "darkred") +
        guides(colour = guide_legend(override.aes = list(size=5, stroke=0))) +
        theme_light(base_family = "sans", base_line_size = 0.5, base_rect_size = 0.1, base_size = 14) +
        ggtitle(paste0("ModuleScore of GI markers")) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
dev.off()


