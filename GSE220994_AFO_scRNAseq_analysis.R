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
    library(DevKidCC)

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
lapply(c(1:6), function(i){
  if (i==1) {raw_sample_name <- "OU1"}
  if (i==2) {raw_sample_name <- "OU2"}
  if (i==3) {raw_sample_name <- "OD3"}
  if (i==4) {raw_sample_name <- "OD4"}
  if (i==5) {raw_sample_name <- "OD5"}
  if (i==6) {raw_sample_name <- "OD6"}

  matrix_directory = paste0(matrix_base_directory, raw_sample_name, "/", raw_sample_name, "/outs/filtered_feature_bc_matrix/")
  message(paste0(raw_sample_name, " beginning"))

  seurat_object <- CreateSeuratObject(Read10X(matrix_directory))
  
  seurat_object <- RenameCells(seurat_object, add.cell.id = raw_sample_name)
  seurat_object <- RenameCells(seurat_object, new.names = gsub("_", "-", colnames(seurat_object)))

  dir.create(paste0(working_directory, "Quality_Control", sep="/", raw_sample_name))
  setwd(paste0(working_directory, "Quality_Control", sep="/", raw_sample_name))

  seurat_object@meta.data[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "MT-")[,1]
  seurat_object@meta.data[["HSP"]] <- PercentageFeatureSet(seurat_object, pattern = "HSP")[,1]
  seurat_object@meta.data[["Ribo1"]] <- PercentageFeatureSet(seurat_object, pattern = "RPL")[,1]
  seurat_object@meta.data[["Ribo2"]] <- PercentageFeatureSet(seurat_object, pattern = "RPS")[,1]
  seurat_object@meta.data[["HB"]] <- PercentageFeatureSet(seurat_object, pattern = "HB")[,1]

  #define genes that drive lots of technical variance & remove from sample
  hkGeneREGEX=c("DNAJ", "EIF", "RPL", "RPS", "RPN1", "POLR", "SNX", "HSP", "H1FX",
                "H2org", "PRKA", "NDUF", "ATP", "PSM", "UBA", "UBE", "USP", "TXN")
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

  metadata <- seurat_object@meta.data
  metadata$Barcode <- rownames(metadata)

  seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)

    seurat_object <- subset(seurat_object, subset = nFeature_RNA > 150)

  barcodes = colnames(seurat_object)

  assign((paste0("seurat_object", sep="_", raw_sample_name)), seurat_object, envir=.GlobalEnv)
  assign((paste0("metadata", sep="_", raw_sample_name)), seurat_object@meta.data, envir=.GlobalEnv)
  assign((paste0("seurat_object_barcodes", sep="_", raw_sample_name)), barcodes, envir=.GlobalEnv)
  message(paste0(raw_sample_name, " complete"))
})

 #### Combine my data
 seurat_object_OU1@meta.data[["orig_sample"]] <- "OU1"
 seurat_object_OU2@meta.data[["orig_sample"]] <- "OU2"
 seurat_object_OD3@meta.data[["orig_sample"]] <- "OD3"
 seurat_object_OD4@meta.data[["orig_sample"]] <- "OD4"
 seurat_object_OD5@meta.data[["orig_sample"]] <- "OD5"
 seurat_object_OD6@meta.data[["orig_sample"]] <- "OD6"

 seurat_object_OU1@meta.data[["diff"]] <- "Undifferentiated"
 seurat_object_OU2@meta.data[["diff"]] <- "Undifferentiated"
 seurat_object_OD3@meta.data[["diff"]] <- "Differentiated"
 seurat_object_OD4@meta.data[["diff"]] <- "Differentiated"
 seurat_object_OD5@meta.data[["diff"]] <- "Differentiated"
 seurat_object_OD6@meta.data[["diff"]] <- "Differentiated"

 seurat_combined_list <- list()
 seurat_combined_list[["OU1"]] <- seurat_object_OU1
 seurat_combined_list[["OU2"]] <- seurat_object_OU2
 seurat_combined_list[["OD3"]] <- seurat_object_OD3
 seurat_combined_list[["OD4"]] <- seurat_object_OD4
 seurat_combined_list[["OD5"]] <- seurat_object_OD5
 seurat_combined_list[["OD6"]] <- seurat_object_OD6

 features <- SelectIntegrationFeatures(object.list = seurat_combined_list, nfeatures = num_DEGS_to_use)
 anchors <- FindIntegrationAnchors(object.list = seurat_combined_list, anchor.features = features)
 org_seurat_object <- IntegrateData(anchorset = anchors)
 DefaultAssay(org_seurat_object) <- "integrated"
 message("Integration complete")
 num_cells_analysed=length(colnames(org_seurat_object))
 
 org_seurat_object_processed <- org_seurat_object

 
# Remove doublets
donor_ids_O123456 <- read.table(file=paste0(input_path, "GSE220994_Organoids_scRNAseq_extended_metadata.csv"), header = T)

org_seurat_object_processed <- subset(org_seurat_object_processed, cells=donor_ids_O123456[donor_ids_O123456$patient!="doublet", "cell"])

# Select epithelial cells
DefaultAssay(org_seurat_object_processed) <- "integrated"
org_seurat_object_processed <- FindVariableFeatures(object = org_seurat_object_processed, selection.method="vst", nfeatures = num_DEGS_to_use)
org_seurat_object_processed <- ScaleData(org_seurat_object_processed, features = rownames(org_seurat_object_processed), do.center = TRUE, do.scale = TRUE)
org_seurat_object_processed <- RunPCA(org_seurat_object_processed, features = VariableFeatures(object = org_seurat_object_processed), npcs = 50)
number_of_PCAs = 20
org_seurat_object_processed <- RunUMAP(org_seurat_object_processed, dims = 1:number_of_PCAs)
org_seurat_object_processed <- RunTSNE(object = org_seurat_object_processed, check_duplicates = FALSE)
org_seurat_object_processed <- FindNeighbors(org_seurat_object_processed, dims = 1:number_of_PCAs)
org_seurat_object_processed <- FindClusters(org_seurat_object_processed, resolution = 0.7)

DefaultAssay(org_seurat_object_processed) <- "RNA"

donor_ids_O123456$patient <- "Unknown"
donor_ids_O123456[donor_ids_O123456$donor_id=="b993", "patient"] <- "993"
donor_ids_O123456[donor_ids_O123456$donor_id=="b10CDH", "patient"] <- "10CDH"
donor_ids_O123456[donor_ids_O123456$donor_id=="b6SB", "patient"] <- "6SB"
donor_ids_O123456[donor_ids_O123456$donor_id=="b724", "patient"] <- "724"
donor_ids_O123456[donor_ids_O123456$donor_id=="b730", "patient"] <- "730"
donor_ids_O123456[donor_ids_O123456$donor_id=="b760", "patient"] <- "760"
donor_ids_O123456[donor_ids_O123456$donor_id=="b744", "patient"] <- "744"
donor_ids_O123456[donor_ids_O123456$donor_id=="b773", "patient"] <- "773"
donor_ids_O123456[donor_ids_O123456$donor_id=="b777", "patient"] <- "777"
donor_ids_O123456[donor_ids_O123456$donor_id=="b768", "patient"] <- "768"
donor_ids_O123456[donor_ids_O123456$donor_id=="b742", "patient"] <- "742"
donor_ids_O123456[donor_ids_O123456$donor_id=="b727", "patient"] <- "727"

barcodes_993 <- donor_ids_O123456[donor_ids_O123456$patient=="993", "cell"]
barcodes_10CDH <- donor_ids_O123456[donor_ids_O123456$patient=="10CDH", "cell"]
barcodes_6SB <- donor_ids_O123456[donor_ids_O123456$patient=="6SB", "cell"]
barcodes_724 <- donor_ids_O123456[donor_ids_O123456$patient=="724", "cell"]
barcodes_727 <- donor_ids_O123456[donor_ids_O123456$patient=="727", "cell"]
barcodes_730 <- donor_ids_O123456[donor_ids_O123456$patient=="730", "cell"]
barcodes_742 <- donor_ids_O123456[donor_ids_O123456$patient=="742", "cell"]
barcodes_744 <- donor_ids_O123456[donor_ids_O123456$patient=="744", "cell"]
barcodes_760 <- donor_ids_O123456[donor_ids_O123456$patient=="760", "cell"]
barcodes_768 <- donor_ids_O123456[donor_ids_O123456$patient=="768", "cell"]
barcodes_773 <- donor_ids_O123456[donor_ids_O123456$patient=="773", "cell"]
barcodes_777 <- donor_ids_O123456[donor_ids_O123456$patient=="777", "cell"]

#OU1 = 993, 10CDH, 724, 727, 730, 744, 760, 773
barcodes_OU1_993 <- barcodes_993[substr(barcodes_993, 1,3)=="OU1"]
barcodes_OU1_10CDH <- barcodes_10CDH[substr(barcodes_10CDH, 1,3)=="OU1"]
barcodes_OU1_724 <- barcodes_724[substr(barcodes_724, 1,3)=="OU1"]
barcodes_OU1_727 <- barcodes_727[substr(barcodes_727, 1,3)=="OU1"]
barcodes_OU1_730 <- barcodes_730[substr(barcodes_730, 1,3)=="OU1"]
barcodes_OU1_744 <- barcodes_744[substr(barcodes_744, 1,3)=="OU1"]
barcodes_OU1_760 <- barcodes_760[substr(barcodes_760, 1,3)=="OU1"]
barcodes_OU1_773 <- barcodes_773[substr(barcodes_773, 1,3)=="OU1"]

#OU2 = 993, 6SB, 730, 742, 760, 768 (few cells), 773, 777
barcodes_OU2_993 <- barcodes_993[substr(barcodes_993, 1,3)=="OU2"]
barcodes_OU2_6SB <- barcodes_6SB[substr(barcodes_6SB, 1,3)=="OU2"]
barcodes_OU2_730 <- barcodes_730[substr(barcodes_730, 1,3)=="OU2"]
barcodes_OU2_742 <- barcodes_742[substr(barcodes_742, 1,3)=="OU2"]
barcodes_OU2_760 <- barcodes_760[substr(barcodes_760, 1,3)=="OU2"]
barcodes_OU2_768 <- barcodes_768[substr(barcodes_768, 1,3)=="OU2"]
barcodes_OU2_773 <- barcodes_773[substr(barcodes_773, 1,3)=="OU2"]
barcodes_OU2_777 <- barcodes_777[substr(barcodes_777, 1,3)=="OU2"]

#OD3 = 993, 6SB, 724, 730, 760, 773, 777, 10SB_AF
barcodes_OD3_993 <- barcodes_993[substr(barcodes_993, 1,3)=="OD3"]
barcodes_OD3_6SB <- barcodes_6SB[substr(barcodes_6SB, 1,3)=="OD3"]
barcodes_OD3_724 <- barcodes_724[substr(barcodes_724, 1,3)=="OD3"]
barcodes_OD3_727 <- barcodes_727[substr(barcodes_727, 1,3)=="OD3"]
barcodes_OD3_730 <- barcodes_730[substr(barcodes_730, 1,3)=="OD3"]
barcodes_OD3_760 <- barcodes_760[substr(barcodes_760, 1,3)=="OD3"]
barcodes_OD3_773 <- barcodes_773[substr(barcodes_773, 1,3)=="OD3"]
barcodes_OD3_777 <- barcodes_777[substr(barcodes_777, 1,3)=="OD3"]
barcodes_OD3_10SB_AF <- barcodes_10SB_AF[substr(barcodes_10SB_AF, 1,3)=="OD3"]

#OD4 = 993, 10CDH, 727, 730, 760, 773, 10SB_AF
barcodes_OD4_993 <- barcodes_993[substr(barcodes_993, 1,3)=="OD4"]
barcodes_OD4_10CDH <- barcodes_10CDH[substr(barcodes_10CDH, 1,3)=="OD4"]
barcodes_OD4_730 <- barcodes_730[substr(barcodes_730, 1,3)=="OD4"]
barcodes_OD4_760 <- barcodes_760[substr(barcodes_760, 1,3)=="OD4"]
barcodes_OD4_773 <- barcodes_773[substr(barcodes_773, 1,3)=="OD4"]
barcodes_OD4_10SB_AF <- barcodes_10SB_AF[substr(barcodes_10SB_AF, 1,3)=="OD4"]
barcodes_OD4_727 <- barcodes_727[substr(barcodes_727, 1,3)=="OD4"]

#OD5 = 6SB, 727, 730, 744, 760, 773, 10SB_AF
barcodes_OD5_6SB <- barcodes_6SB[substr(barcodes_6SB, 1,3)=="OD5"]
barcodes_OD5_727 <- barcodes_727[substr(barcodes_727, 1,3)=="OD5"]
barcodes_OD5_730 <- barcodes_730[substr(barcodes_730, 1,3)=="OD5"]
barcodes_OD5_744 <- barcodes_744[substr(barcodes_744, 1,3)=="OD5"]
barcodes_OD5_760 <- barcodes_760[substr(barcodes_760, 1,3)=="OD5"]
barcodes_OD5_773 <- barcodes_773[substr(barcodes_773, 1,3)=="OD5"]
barcodes_OD5_10SB_AF <- barcodes_10SB_AF[substr(barcodes_10SB_AF, 1,3)=="OD5"]

#OD6 = 10CDH, 727, 730, 744, 760, 777, 10SB_AF
barcodes_OD6_10CDH <- barcodes_10CDH[substr(barcodes_10CDH, 1,3)=="OD6"]
barcodes_OD6_727 <- barcodes_727[substr(barcodes_727, 1,3)=="OD6"]
barcodes_OD6_730 <- barcodes_730[substr(barcodes_730, 1,3)=="OD6"]
barcodes_OD6_744 <- barcodes_744[substr(barcodes_744, 1,3)=="OD6"]
barcodes_OD6_760 <- barcodes_760[substr(barcodes_760, 1,3)=="OD6"]
barcodes_OD6_777 <- barcodes_777[substr(barcodes_777, 1,3)=="OD6"]
barcodes_OD6_10SB_AF <- barcodes_10SB_AF[substr(barcodes_10SB_AF, 1,3)=="OD6"]

# Label individual organoid lines
metadata <- org_seurat_object_processed@meta.data
metadata$Barcode <- rownames(metadata)
metadata[,"organoid_line"] <- "Unassigned"
metadata[barcodes_OU1_993,"organoid_line"] <- "993-SiAFO"
metadata[barcodes_OU1_724,"organoid_line"] <- "724-KAFO"
metadata[barcodes_OU1_727,"organoid_line"] <- "727-LAFO"
metadata[barcodes_OU1_730,"organoid_line"] <- "730-cdhLTFO"
metadata[barcodes_OU1_744,"organoid_line"] <- "744-KAFO"
metadata[barcodes_OU1_760,"organoid_line"] <- "760-cdhLAFO"
metadata[barcodes_OU1_773,"organoid_line"] <- "773-LAFO"
metadata[barcodes_OU1_10CDH,"organoid_line"] <- "10CDH-cdhLTFO"

metadata[barcodes_OU2_993,"organoid_line"] <- "993-SiAFO"
metadata[barcodes_OU2_730,"organoid_line"] <- "730-cdhLTFO"
metadata[barcodes_OU2_742,"organoid_line"] <- "742-KAFO"
metadata[barcodes_OU2_760,"organoid_line"] <- "760-cdhLAFO"
metadata[barcodes_OU2_768,"organoid_line"] <- "768-cdhLAFO"
metadata[barcodes_OU2_773,"organoid_line"] <- "773-SiAFO"
metadata[barcodes_OU2_777,"organoid_line"] <- "777-LAFO"
metadata[barcodes_OU2_6SB,"organoid_line"] <- "6SB-KAFO"

metadata[barcodes_OD3_993,"organoid_line"] <- "993-SiAFO"
metadata[barcodes_OD3_724,"organoid_line"] <- "724-KAFO"
metadata[barcodes_OD3_730,"organoid_line"] <- "730-CDH_proximal_LTFO"
metadata[barcodes_OD3_760,"organoid_line"] <- "760-CDH_proximal_LAFO"
metadata[barcodes_OD3_773,"organoid_line"] <- "773-SiAFO"
metadata[barcodes_OD3_777,"organoid_line"] <- "777-proximal_LAFO"
metadata[barcodes_OD3_6SB,"organoid_line"] <- "6SB-KAFO"
metadata[barcodes_OD3_6SB,"diff"] <- "Undifferentiated"

metadata[barcodes_OD4_993,"organoid_line"] <- "993-SiAFO"
metadata[barcodes_OD4_727,"organoid_line"] <- "727-distal_LAFO"
metadata[barcodes_OD4_730,"organoid_line"] <- "730-CDH_distal_LTFO"
metadata[barcodes_OD4_760,"organoid_line"] <- "760-CDH_distal_LAFO"
metadata[barcodes_OD4_773,"organoid_line"] <- "773-proximal_LAFO"
metadata[barcodes_OD4_10CDH,"organoid_line"] <- "10CDH-CDH_distal_LTFO"

metadata[barcodes_OD5_727,"organoid_line"] <- "727-proximal_LAFO"
metadata[barcodes_OD5_730,"organoid_line"] <- "730-CDH_proximal_LTFO"
metadata[barcodes_OD5_744,"organoid_line"] <- "744-KAFO"
metadata[barcodes_OD5_760,"organoid_line"] <- "760-CDH_proximal_LAFO"
metadata[barcodes_OD5_773,"organoid_line"] <- "773-distal_LAFO"
metadata[barcodes_OD5_6SB,"organoid_line"] <- "6SB-KAFO"

metadata[barcodes_OD6_727,"organoid_line"] <- "727-KAFO"
metadata[barcodes_OD6_727,"diff"] <- "Undifferentiated"
metadata[barcodes_OD6_730,"organoid_line"] <- "730-CDH_distal_LTFO"
metadata[barcodes_OD6_744,"organoid_line"] <- "744-LAFO"
metadata[barcodes_OD6_744,"diff"] <- "Undifferentiated"
metadata[barcodes_OD6_760,"organoid_line"] <- "760-CDH_distal_LAFO"
metadata[barcodes_OD6_777,"organoid_line"] <- "777-distal_LAFO"
metadata[barcodes_OD6_10CDH,"organoid_line"] <- "10CDH-CDH_proximal_LTFO"

metadata$patient <- substr(metadata$organoid_line, 1, 4)
metadata$patient <- gsub("-", "", metadata$patient)
metadata$patient <- gsub("10CD", "10CDH", metadata$patient)
metadata$patient <- gsub("Unas", "Unassigned", metadata$patient)
metadata$identity <- substr(metadata$organoid_line, 5, nchar(metadata$organoid_line))
metadata$identity <- gsub("H-", "", metadata$identity)
metadata$identity <- gsub("signed", "Unassigned", metadata$identity)
metadata$lineage <- substr(metadata$identity, (nchar(metadata$identity)-3), nchar(metadata$identity))
metadata$lineage <- gsub("gned", "Unassigned", metadata$lineage)
metadata$lineage <- gsub("iAFO", "SiAFO", metadata$lineage)

org_seurat_object_processed@meta.data <- metadata[colnames(org_seurat_object_processed),]

org_seurat_object_processed_preAnalysis <- org_seurat_object_processed


org_UMAPdata=as.data.frame(org_seurat_object_processed@meta.data)
org_UMAPdata$UMAP_1 <- org_seurat_object_processed@reductions[["umap"]]@cell.embeddings[,"UMAP_1"]
org_UMAPdata$UMAP_2 <- org_seurat_object_processed@reductions[["umap"]]@cell.embeddings[,"UMAP_2"]
org_UMAPdata$Cluster=org_seurat_object_processed@meta.data$seurat_clusters
org_UMAPdata$Barcode <- rownames(org_UMAPdata)
org_UMAPdata_backup <- org_UMAPdata


# Select LAFO
Idents(org_seurat_object_processed) <- org_seurat_object_processed$identity
org_seurat_object_lung <- subset(x = org_seurat_object_processed, idents = c("proximal_LAFO", "distal_LAFO", "LAFO"), invert = FALSE)
saveRDS(org_seurat_object_lung, paste0(working_directory, "../LAFO_seurat_object_", ceiling(length(colnames(org_seurat_object_lung))/1000),"k.rds"))

# Select LAFO & CDH
Idents(org_seurat_object_processed) <- org_seurat_object_processed$identity
org_seurat_object_lung <- subset(x = org_seurat_object_processed, idents = c("proximal_LAFO", "distal_LAFO", "LAFO",
                                                                             "CDH_distal_LAFO", "CDH_distal_LTFO", "CDH_proximal_LAFO",
                                                                             "CDH_proximal_LTFO", "cdhLAFO", "cdhLTFO"), invert = FALSE)
saveRDS(org_seurat_object_lung, paste0(working_directory, "../LAFO_CDH_seurat_object_", ceiling(length(colnames(org_seurat_object_lung))/1000),"k.rds"))

# Select KAFO
Idents(org_seurat_object_processed) <- org_seurat_object_processed$identity
org_seurat_object_kidney <- subset(x = org_seurat_object_processed, idents = "KAFO", invert = FALSE)
saveRDS(org_seurat_object_kidney, paste0(working_directory, "../KAFO_seurat_object_", ceiling(length(colnames(org_seurat_object_kidney))/1000),"k.rds"))

# Select SiAFO
Idents(org_seurat_object_processed) <- org_seurat_object_processed$identity
org_seurat_object_GI <- subset(x = org_seurat_object_processed, idents = "SiAFO", invert = FALSE)
saveRDS(org_seurat_object_GI, paste0(working_directory, "../SiAFO_seurat_object_", ceiling(length(colnames(org_seurat_object_GI))/1000),"k.rds"))

