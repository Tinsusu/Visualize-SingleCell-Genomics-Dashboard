
###### contain all over backendcode here ####
library(shiny)
library(shinydashboard)
library(shinyjs)
library(tools)

#when we upload the file it will be kept in temp folder

# Read in file and perform validation.
load_seurat_obj <- function(path){
  errors <- c()
  # check file extension
  if (!tolower(tools::file_ext(path)) == "rds") { # ignores case
    errors <- c(errors, "Invalid rds file.")
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
  
  # Validate obj is a seurat object
  if (!inherits(obj, "Seurat")){
    errors <- c(errors, "File is not a seurat object")
    return(errors)
  }
  
  return(obj)
}

 


