
# STEP 1: LOAD DATA-------------------------------------------------------------------------------------

library(tidyverse)
library(writexl)
library(Seurat)
library(devtools)
library(presto)

# --- SETTINGS ---
min_cells <- 5          # A gene must appear in at least 5 cells to be kept


# --- LOAD THE .h5 FILE ---
# Read10X_h5() reads the compressed CellRanger output file
# It returns a matrix: rows = genes, columns = cell barcodes

raw.data <- Read10X_h5(
  "sample_name_filtered_feature_bc_matrix.h5"
)

# --- CREATE THE SEURAT OBJECT ---
raw.so <- CreateSeuratObject(
  counts    = raw.data,    # the raw count matrix from above
  project   = "sample_name",       # sample name — will appear in plots
  min.cells = min_cells       # filter: keep genes seen in ≥5 cells
)

# --- CHECK THE OBJECT ---
raw.so



# STEP 2: QUALITY CONTROL------------------------------------------------------------------------------------


# --- CALCULATE QC METRICS ---

# Percentage of reads from mitochondrial genes
# High mt% = dying/damaged cell (its cytoplasm leaked but mitochondria stayed)
raw.so[["percent.mt"]] <- PercentageFeatureSet(raw.so, pattern = "^MT-")

# Percentage of reads from ribosomal genes
# Useful to track but we don't usually filter on this
raw.so[["percent.rb"]] <- PercentageFeatureSet(raw.so, pattern = "^RP[SL]")

# --- INSPECT THE METADATA TABLE ---
head(raw.so@meta.data, 10)   # first 10 cells
tail(raw.so@meta.data, 10)   # last 10 cells

# --- VISUALIZE QC METRICS ---
p1 <- VlnPlot(raw.so,
              features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
              ncol = 3,
              pt.size = 0.1)   # show individual cells as dots

p2 <- VlnPlot(raw.so,
              features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
              ncol = 3,
              pt.size = 0)     # no dots — cleaner view

p1 / p2   # display both plots stacked

ggsave("sample_name_dots_violin_plot.png", plot = p1, bg = "white", width = 12, height = 6)
ggsave("sample_name_no_dots_violin_plot.png", plot = p2, bg = "white", width = 12, height = 6)

# --- CALCULATE QUANTILES (to help choose filter thresholds) ---

# Mitochondrial %
quantile(raw.so@meta.data$percent.mt, seq(0, 1, 0.1))    # every 10%
quantile(raw.so@meta.data$percent.mt, seq(0.9, 1, 0.01)) # top 10%, fine-grained

# RNA counts (UMIs)
quantile(raw.so@meta.data$nCount_RNA, seq(0, 1, 0.05))   # every 5%
quantile(raw.so@meta.data$nCount_RNA, seq(0.99, 1, 0.001)) # top 1%

# Gene counts
quantile(raw.so@meta.data$nFeature_RNA, seq(0, 1, 0.05))    # every 5%
quantile(raw.so@meta.data$nFeature_RNA, seq(0, 0.05, 0.001)) # bottom 5%


# STEP 3: FILTER CELLS & SAVE----------------------------------------------------------------------------------------


# --- DEFINE THRESHOLDS ---
percent.mt_max   <- 10     # remove cells with >15% mitochondrial reads (dead/damaged)
nCount_RNA_min   <- 600    # remove cells with <600 UMIs (empty droplets)
nCount_RNA_max   <- 13000  # remove cells with >13000 UMIs (doublets)
nFeature_RNA_min <- 300    # remove cells with <300 genes (empty droplets)

# --- COUNT CELLS BEFORE FILTERING ---
nCells_before <- nrow(raw.so@meta.data)

# --- APPLY FILTERS ---
sample_name.so <- subset(raw.so, 
                    percent.mt  <= percent.mt_max   &
                      nCount_RNA  >= nCount_RNA_min   &
                      nCount_RNA  <= nCount_RNA_max   &
                      nFeature_RNA >= nFeature_RNA_min)

# --- COUNT CELLS AFTER FILTERING ---
nCells_after <- nrow(sample_name.so@meta.data)

# --- REPORT HOW MANY CELLS WERE REMOVED ---
message("Cells BEFORE filtering: ", nCells_before)
message("Cells AFTER  filtering: ", nCells_after)
message("Cells removed: ",          nCells_before - nCells_after,
        " (", round((nCells_before - nCells_after) / nCells_before * 100, 1), "%)")

# --- QUICK CHECK ---
head(sample_name.so@meta.data, 10)
sample_name.so

# --- SAVE THE FILTERED OBJECT ---
saveRDS(sample_name.so, file = "sample_name.so_qc_passed.rds")



# STEP 4: SCTRANSFORM NORMALIZATION---------------------------------------------------------------------------------------


# This single function does three things:
# 1. Normalizes the data (fixes differences in sequencing depth)
# 2. Scales the data (prepares it for PCA)
# 3. Finds Variable Features (the 3000 most interesting genes)

sample_name.so <- SCTransform(sample_name.so, 
                         vars.to.regress = "percent.mt", 
                         verbose = FALSE)

# --- VISUALIZE VARIABLE FEATURES ---
# These are the genes that vary the most between your cells
top10 <- head(VariableFeatures(sample_name.so), 10)
topgenes <- VariableFeaturePlot(sample_name.so)
topgeneslabeled <- LabelPoints(plot = topgenes, points = top10, repel = TRUE)

topgeneslabeled

ggsave("sample_name_hvg_plot.png", plot = topgeneslabeled, bg = "white", width = 10, height = 7)

# --- SAVE PROGRESS ---
saveRDS(sample_name.so, file = "sample_name_SCT_qc_passed.so.rds")


# STEP 5: PCA & ELBOW PLOT----------------------------------------------------------------------


sample_name.so <- RunPCA(sample_name.so, verbose = FALSE)

# Visualize the first few PCs and their top genes
p1 <- VizDimLoadings(sample_name.so, dims = 1:2, reduction = "pca")
ggsave("sample_name_PCA_Loadings.png", plot = p1, width = 10, height = 7, bg = "white")

p2 <- DimPlot(sample_name.so, reduction = "pca") + NoLegend()
ggsave("sample_name_PCA_Plot.png", plot = p2, width = 7, height = 7, bg = "white")

# This is the most important plot for the next step:
p3 <- ElbowPlot(sample_name.so, ndims = 50)
ggsave("sample_name_ElbowPlot.png", plot = p3, width = 8, height = 5, bg = "white")



# STEP 6: UMAP & CLUSTERING----------------------------------------------------------------------------


# 1. Find Neighbors: Builds a "social network" of cells
# We use the first PCs based on our Elbow Plot, need to change dims
sample_name.so <- FindNeighbors(sample_name.so, k.param = 20, dims = 1:25)

# 2. Find Clusters: Groups the cells
# Resolution 0.5 is a good starting point (higher = more clusters)
sample_name.so <- FindClusters(sample_name.so, resolution = 0.4, algorithm = 3) # algorithm 3 = Leiden

# 3. Run UMAP: Creates the 2D visualization
sample_name.so <- RunUMAP(sample_name.so, dims = 1:25, n.neighbors = 20, min.dist = 0.3, seed.use = 42)

# 4. Visualize and Save
p_umap <- DimPlot(sample_name.so, reduction = "umap", label = TRUE) + 
  ggtitle("UMAP Clustering - Sample sample_name")

p_umap

ggsave("sample_name_UMAP1_Clusters.png", plot = p_umap, width = 8, height = 7, bg = "white")



####SAVE RDS for ANNOTATION!!!!----------

saveRDS(sample_name.so, "sample_name_clustered.rds")


# STEP 7: FIND MARKER GENES --------------------------------------------------------------------------------------------------


# --- FIND MARKER GENES FOR EVERY CLUSTER ---
# For each cluster, finds genes that are MORE expressed there
# compared to ALL other clusters combined

sample_name.markers <- FindAllMarkers(
  sample_name.so,
  only.pos       = TRUE,   # only keep UPREGULATED markers (positive markers)
  min.pct        = 0.1,    # gene must be detected in ≥10% of cells in the cluster
  logfc.threshold = 0.25   # minimum log2 fold change (expression difference)
)

# --- QUICK LOOK AT THE RESULTS ---
sample_name.markers %>% as_tibble()

#--- Extract top markers per cluster---
# --- TOP 10 MARKERS PER CLUSTER (by fold change) ---
top10_per_cluster <- sample_name.markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10) %>%
  ungroup()

top10_per_cluster


# --- TOP MARKER PER CLUSTER (quick overview) ---
sample_name.markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 1)

# --- LOOK UP A SPECIFIC CLUSTER ---
sample_name.markers %>%
  as_tibble() %>%
  filter(cluster == 0) %>%   # change the number to whichever cluster you want
  head(20)

# --- LOOK UP A SPECIFIC GENE ---
# Useful to check CLL marker genes!
sample_name.markers %>%
  as_tibble() %>%
  filter(gene == "CD5")      # try also: CD19, CXCR4, LEF1, ROR1, MS4A1


#--- Save the results ---

# --- SAVE AS .rds (to reload in R later) ---
saveRDS(sample_name.markers, file = "sample_name.markers.rds")

# --- SAVE AS EXCEL (to open in Excel) ---
# All markers
write_xlsx(sample_name.markers, "sample_name_AllMarkers.xlsx")

# Top 10 per cluster only
write_xlsx(top10_per_cluster, "sample_name_Top10_per_Cluster.xlsx")



# STEP 8: VISUALIZATION------------------------------------------------------------------------------------------------------------------------------------------------------------



###### markers_heatmap ######

top10 <- sample_name.markers %>%
  filter(avg_log2FC > 1) %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10) %>%
  ungroup()

# Change to RNA assay, normalization and scaling
DefaultAssay(sample_name.so) <- "RNA"

sample_name.so <- NormalizeData(sample_name.so)

sample_name.so <- ScaleData(sample_name.so, features = top10$gene)

markers_heatmap <- DoHeatmap(sample_name.so, 
                             features = top10$gene, 
                             group.by = "seurat_clusters") + NoLegend()

markers_heatmap

ggsave("sample_name_Top_Markers_Heatmap.png", plot = markers_heatmap, width = 14, height = 10, bg = "white")


# See exactly where a classic B-cell or T-cell  marker lights up
FeaturePlot(sample_name.so, features = c("----", "----"))

# Look at the distribution of those same markers across clusters 0-8
VlnPlot(sample_name.so, features = c("----", "---"))

#Proliferating center

CD5_CXCR4<- FeaturePlot(sample_name.so, 
                        features = c("CD5", "CXCR4"), 
                        blend = TRUE, 
                        pt.size = 0.3)
ggsave("sample_name_CD5_CXCR.png", plot= CD5_CXCR4, bg = "white",width = 10, height = 7)



####################################### STEP 9: ANNOTATION FIRST ########################################



# Switch to the RNA assay for better biological visualization 
DefaultAssay(sample_name.so.annotated) <- "RNA"
sample_name.so.annotated <- NormalizeData(sample_name.so.annotated)

# Plot key CLL markers
# CD5 is the most important (T-cell marker normally, but expressed in CLL B-cells)
FeaturePlot(sample_name.so.annotated, 
            features = c("CD5", "MS4A1", "LEF1", "ROR1"), 
            ncol = 2, 
            reduction = "umap")



# Step 10: Find the cell----------------------------------------------------------------------------------------------------------------------------



# Find the barcode for the cell 
find_cell <- WhichCells(sample_name.so.annotated, expression = predicted.celltype.l2 == "Treg")

# Print the barcode to your console
find_cell

# Look at the UMIs, Total Genes, and Mitochondrial % 
sample_name.so.annotated@meta.data[find_cell, c("nCount_RNA", "nFeature_RNA", "percent.mt")]

# Spot it on UMAP
# Draw the UMAP highlighting ONLY our target cell
DimPlot(
  sample_name.so.annotated, 
  reduction = "umap", 
  cells.highlight = find_cell, 
  cols.highlight = "red", 
  sizes.highlight = 3  # Make it extra big so you can't miss it!
) + 
  NoLegend() +
  ggtitle("Location of the cell")

# The genes expressed in find_cell compared to the dataset, log2fc>1

# 1. Set the active identities to your fine-grained labels
Idents(sample_name.so.annotated) <- "predicted.celltype.l2"

# 2. Run the comparison using the existing Treg label!
targeted_comparison <- FindMarkers(
  sample_name.so.annotated,
  ident.1 = "Treg",                                            # Your single cell
  ident.2 = c("CD4 Naive", "CD4 TCM", "CD4 CTL", "CD4 TEM"),   # Grouping the normal CD4s together
  logfc.threshold = 0.5,
  min.cells.group = 1
)

# 3. View the results to see the biological noise
head(targeted_comparison, 10)



# Step 11: Find markers in Azimuth---------------------------------------------------------------------------------------------------------------------------------



# Find markers at all clusters
# Ensure Seurat is focused on the Level 2 annotations
Idents(sample_name.so.annotated) <- "predicted.celltype.l2"

# Run FindAllMarkers across the entire annotated dataset
azimuth_l2_markers <- FindAllMarkers(
  sample_name.so.annotated,
  only.pos = TRUE,
  min.pct = 0.1,
  logfc.threshold = 0.25
)

# Find markers at one specific cluster
# Find the distinguishing genes for a single cell type
#specific_celltype_markers <- FindMarkers(
 # sample_name.so.annotated,
  #ident.1 = "B intermediate",  
  #only.pos = TRUE,             # Only show genes that are turned UP in this cell type
  #min.pct = 0.1,
  #logfc.threshold = 0.25)


# Extract and view the Top 10 genes for EVERY cell type
top10_azimuth_markers <- azimuth_l2_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10) %>%
  ungroup()

# Print the list to the console
top10_azimuth_markers

# Save this highly valuable biological list to Excel!
write_xlsx(top10_azimuth_markers, "sample_name_Azimuth_L2_Top10_Markers.xlsx")
