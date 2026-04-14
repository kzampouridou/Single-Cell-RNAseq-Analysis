####Intergration GEX/BCR###

library(Seurat)
library(scRepertoire)
library(readxl)

# 1. LOAD DATA
# Check if we have the annotated file     
seurat_obj <- readRDS("/Sample.annotated.rds")
class(seurat_obj)

# Check if we have the output of scRepertoire
bcr_df <- read_excel("/combined.Sample_BCR.clustered.xlsx")
bcr_list <- list(bcr_df)
names(bcr_list) <- "E22150_BCR"

# 2. BARCODE ALIGMENT CHECK
# Seurat barcodes must match BCR barcodes exactly

message("Seurat Barcode Example: ", head(Cells(seurat_obj), 1))
message("BCR Barcode Example: ", head(bcr_list[[1]]$barcode, 1))

# Rename the Seurat/BCR cells to add the prefix matching the BCR/Seurat data if needed
# This turns "AAAC..." into "E22150_BCR_AAAC..."
seurat_obj <- RenameCells(seurat_obj, add.cell.id = "Sample_BCR")

# 3. INTERGRATION
seurat_integrated <- combineExpression(
  bcr_list, 
  seurat_obj, 
  cloneCall = "strict", 
  chain = "both",
  proportion = FALSE,
  cloneSize = c(Single = 1, Small = 5, Medium = 20, Large = 100)
)

# 4. ANALYSIS & VISUALIZATION
# only B-cells will have values in these new columns
head(seurat_integrated@meta.data[, c("CTstrict", "cloneSize")])

# plotting the clonal lineages
p1 <- DimPlot(seurat_integrated, 
              group.by = "cloneSize", 
              reduction = "ref.umap") +ggtitle("BCR Clonal Expansion")

p1


