# We will use Seurat to read in 10X sequencing data and to preprocess them. The output data will be saved in .csv files. 
# 
# Required packages: Seurat(3.2.2), dplyr(1.0.2), data.table(1.13.2)

myworking_directory <- "~/input/expression/"
setwd(myworking_directory)

library(Seurat)
library(dplyr)
library(data.table)

# set up parameters
path.folder = 'raw/'

removeCellCycle = FALSE

if (removeCellCycle){
  condition <- "_cc_regressed_dif"
} else {
  condition <- ""
}

data.names = c('mmB_D5','mmB_DesLO_D11','mmB_DesLO_D17')
cutoff_RNA_ls = c(200,200,1250)
cutoff_RNA_us = c(9000,9000,7500)
cutoff_percents = c(6.5, 7.5, 30)

for(i in 1:3) {
  data.name = data.names[i]
  cutoff_RNA_l = cutoff_RNA_ls[i]
  cutoff_RNA_u = cutoff_RNA_us[i]
  cutoff_percent = cutoff_percents[i]
  
  print(paste(data.name, cutoff_RNA_l, cutoff_RNA_u, cutoff_percent))
  
  # read in parameters 
  project.name = paste0('10X_',data.name)
  
  # read in 10X data   
  current.data = Read10X(data.dir = paste0(path.folder, data.name), gene.column = 2)
  
  # check number of cells & genes 
  print(dim(current.data)) 
  
  # create objects & filter out low-quality genes and cells 
  current = CreateSeuratObject(counts = current.data, min.cells = 0, min.features = 200, project = project.name)
  
  # check number of cells & genes 
  print(dim(current[["RNA"]]@data))
  
  # remove low-quality cells 
  current[["percent.mt"]] <- PercentageFeatureSet(object = current, pattern = "^MT-")
  VlnPlot(object = current, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  
  plot1 <- FeatureScatter(object = current, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(object = current, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  CombinePlots(plots = list(plot1, plot2))
  
  # subset data 
  current <- subset(x = current, subset = nFeature_RNA > cutoff_RNA_l & nFeature_RNA < cutoff_RNA_u & percent.mt < cutoff_percent)
  
  # normalize data 
  current <- NormalizeData(object = current)  # default normalization:TPM (with scaling factor = 1e4) and then log1p()  
  
  # check number of cells & genes 
  print(dim(current[["RNA"]]@data))
  
  # find most variable genes 
  current <- FindVariableFeatures(object = current)  # default returning 2000 features
  
  if (removeCellCycle) {
    # regress out cell-cycle effect 
    s.genes <- cc.genes$s.genes
    g2m.genes <- cc.genes$g2m.genes
  
    # regress out variance between deviding & non-deviding (apply to ALL genes)
    current <- CellCycleScoring(object = current, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
    current$CC.Difference <- current$S.Score - current$G2M.Score
    current <- ScaleData(current, vars.to.regress = "CC.Difference", features = rownames(current))
    
    # check number of cells & genes 
    print(dim(current[["RNA"]]@scale.data))
    
    data_to_save <- current[['RNA']]@scale.data
    
  } else {
    
    data_to_save <- current[['RNA']]@data
    
  }
  
  data_to_write_out <- as.data.frame(as.matrix(data_to_save))
  fwrite(x = data_to_write_out, file = paste0(data.name, condition, ".csv"))
}
