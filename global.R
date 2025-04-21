
###### contain all over backendcode here ####
library(shiny)
library(shinydashboard)
library(shinyjs)
library(tools)
library(Seurat)
library(shinybusy)
library(dplyr)
library(DT)
library(shinydashboardPlus)
library(ggplot2)
library(shinybusy)
library(glue)
library(markdown)
library(ggthemes)
library(grid)
library(png)

#when we upload the file it will be kept in temp folder
## 1. check if it is rds extenstion
## 2. check if it is Seurat objcet
# Read in file and perform validation.
load_seurat_obj <- function(path){
  errors <- c()
  # check file extension
  if (!tolower(tools::file_ext(path)) == "rds") { # ignores case
    errors <- c(errors, "Invalid file we neeed rds.") #add eror in error vector
    return(errors)
  }
  
  # try to read in file
  tryCatch(
    {
      obj <- readRDS(path)
    },
    error = function(e) {
      errors <- c(errors, "Invalid rds file.")
      return(errors)
    }
  )
  
  # Validate obj is a seurat object inherit function is built-in in R
  if (!inherits(obj, "Seurat")){
    errors <- c(errors, "The file is not a seurat object")
    return(errors)
  }
  
  return(obj)
}


create_metadata_UMAP <- function(obj, col){
  if (col %in% c("nCount_RNA", "nFeature_RNA", "percent.mt")){
    feature_value <- obj@meta.data[, col]
    
    # Extract UMAP coordinates and rename columns to UMAP_1 and UMAP_2
    umap_coords <- as.data.frame(obj@reductions$umap@cell.embeddings)
    colnames(umap_coords)[1:2] <- c("UMAP_1", "UMAP_2")
    
    # Combine with feature value
    col_df <- data.frame(umap_coords, feature_value = feature_value)
    
    # Plot
    umap <- ggplot(data = col_df) +
      geom_point(mapping = aes(UMAP_1, UMAP_2, color = log10(feature_value)), size = 0.01) +
      scale_colour_gradientn(colours = rainbow(7))
  } else if (col %in% colnames(obj@meta.data)) {
    umap <- DimPlot(obj, pt.size = .1, label = F, label.size = 4, group.by = col, reduction = "umap")
  } else {
    umap <- ggplot() +
      theme_void() +
      geom_text(aes(x = 0.5, y = 0.5, label = "collumn doesn't exist"), size = 20, color = "gray73", fontface = "bold") +
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  }
  return(umap)
}


create_feature_plot <- function(obj, gene) {
  if (gene %in% rownames(obj)) {
    FP <- FeaturePlot(obj, features = gene, pt.size = 0.001, combine = FALSE)
  } else {
    FP <- ggplot() + 
      theme_void() + 
      geom_text(aes(x = 0.5, y = 0.5, label = "Gene doesn't exist"), size = 20, color = "gray73", fontface = "bold") +
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  }
  return(FP)
 }

 
create_metadata_TSNE <- function(obj, col, dims = 1:10) {
  # Automatically run t-SNE if missing
  if (!("tsne" %in% names(obj@reductions))) {
    message("t-SNE reduction not found. Running RunTSNE() with dims = ", paste(dims, collapse = ", "))
    obj <- RunTSNE(obj, dims = dims)
  }
  
  if (col %in% c("nCount_RNA", "nFeature_RNA", "percent.mt")) {
    feature_value <- obj@meta.data[, col]
    
    # Extract t-SNE coordinates and rename columns
    tsne_coords <- as.data.frame(obj@reductions$tsne@cell.embeddings)
    colnames(tsne_coords)[1:2] <- c("TSNE_1", "TSNE_2")
    
    # Combine with metadata feature
    col_df <- data.frame(tsne_coords, feature_value = feature_value)
    
    # Plot
    tsne <- ggplot(data = col_df) +
      geom_point(mapping = aes(TSNE_1, TSNE_2, color = log10(feature_value)), size = 0.01) +
      scale_colour_gradientn(colours = rainbow(7))
    
  } else if (col %in% colnames(obj@meta.data)) {
    tsne <- DimPlot(obj, pt.size = .1, label = FALSE, label.size = 4, group.by = col, reduction = "tsne")
  } else {
    tsne <- ggplot() +
      theme_void() +
      geom_text(aes(x = 0.5, y = 0.5, label = "col doesn't exist"), size = 20, color = "gray73", fontface = "bold") +
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  }
  
  return(tsne)
}


