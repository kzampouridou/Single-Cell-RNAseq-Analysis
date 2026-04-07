library(tidyverse)
library(Seurat)
library(devtools)
library(writexl)
library(dplyr)

#----------------------------------------------------------------------

#Read the hf5 files from cell ranger for each sample
min_cells <- 5

Sample.data <- Read10X_h5("C:/filtered_feature_bc_matrix.h5")
Sample.so <- CreateSeuratObject(counts = Sample.data, 
                                project = "Sample", 
                                min.cells = min_cells)
Sample.so


#-------------------------------------------------------------------------------------------------------------


#Start the QC by filtering for mitochondrial genes
Sample.so[["percent.mt"]] <- PercentageFeatureSet(Sample.so, pattern = "^MT-")
Sample.so[["percent.rb"]] <- PercentageFeatureSet(Sample.so, pattern = "^RP[SL]")

head(Sample.so@meta.data, 10)    #display first 10 rows of the metadata table
tail(Sample.so@meta.data, 10)    #display last 10 rows of the metadata table

#Visualize with a violin plot for UMI counts (ncounts), Gene counts (nFeature), percent mt
p1 <- VlnPlot(Sample.so, 
              features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
              ncol = 3, 
              pt.size = 0.1)

p2 <- VlnPlot(Sample.so, 
              features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
              ncol = 3, 
              pt.size = 0)

p1 / p2

ggsave("P1_violin_plot.png", plot = p1)
ggsave("P2_violin_plot.png", plot = p2)
ggsave("p1_p2_violin_plot.png", plot = p1|p2)

quantile(Sample.so@meta.data$percent.mt, seq(0,1,0.1))
quantile(Sample.so@meta.data$percent.mt, seq(0.9,1,0.01))
quantile(Sample.so@meta.data$nCount_RNA, seq(0,1,0.05))
quantile(Sample.so@meta.data$nCount_RNA, seq(0.99,1,0.001))
quantile(Sample.so@meta.data$nFeature_RNA, seq(0,1,0.05))
quantile(Sample.so@meta.data$nFeature_RNA, seq(0,0.05,0.001))
quantile(Sample.so@meta.data$nFeature_RNA, seq(0.9,1,0.001))

#Now we can apply filters on the Seurat object 'so' using the subset function. Based on the quantiles percentages (5-95%).
nCells_before <- nrow(Sample.so@meta.data)

percent.mt_max <- 10
nCount_RNA_min <- 600
nCount_RNA_max <- 13000
nFeature_RNA_min <- 400
#nFeature_RNA_max <- 4000

Sample.so <- subset(Sample.so, percent.mt <= percent.mt_max & 
                      nCount_RNA >= nCount_RNA_min &
                      nCount_RNA <= nCount_RNA_max &
                      nFeature_RNA >= nFeature_RNA_min)

nCells_after <- nrow(Sample.so@meta.data)

message("Number of Cells before Filtering: ", nCells_before)
message("Number of Cells after Filtering: ", nCells_after)
message("Percentage of Cells Filtered: ", round((nCells_before - nCells_after) / nCells_before * 100, 1), "%")


#Carry over the filtered cells and save as an .rds file
saveRDS(Sample.so, file = "Sample.so.rds")

Sample.so
head(Sample.so@meta.data, 10) 

#---------------------------------------------------------------------------------------------------------------

#scTransform
#Data normalization (see the markdown file for more info on SCTransform)
so <- readRDS("Sample.so.rds")
options(future.globals.maxSize = 8000 * 1024^2)
#BiocManager::install('glmGamPoi')
so <- SCTransform(so, verbose = TRUE)


#The function VariableFeaturePlot can be used to visualize the results of the HVGs identification.
hvg_plot <- VariableFeaturePlot(so)


plot(hvg_plot)


top10 <- head(VariableFeatures(so), 10)
hvg_plot <- LabelPoints(plot = hvg_plot, points = top10, repel = TRUE)
hvg.plot<- plot(hvg_plot)

ggsave("Sample_hvg_plot.png", plot = hvg.plot, bg = "white", width = 10, height = 7)


#SCTResults(so, slot = "feature.attributes")
#z <- pbmc_small@assays$SCT@SCTModel.list[[1]]@feature.attributes
#qplot(z$gmean, z$residual_variance) + scale_y_log10() + scale_x_log10()
all.genes <- rownames(so)
#IN008458 <- ScaleData(so, features = all.genes)
#The first step of dimensionality reduction is Principal Component Analysis (PCA).
so <- RunPCA(so, npcs = 100, verbose = FALSE)

elbow_plot<- ElbowPlot(so, ndims = 100)
elbow_plot
ggsave("Sample_elbow_plot.png", plot=elbow_plot, bg = "white",width = 10, height = 7)


#---------------------------------------------------------------------------------------------------------------

#UMAP
npcs <- 20

so <- RunUMAP(so, 
              dims = 1:npcs, 
              verbose = TRUE, 
              seed.use = 42)
dimplot<- DimPlot(so, 
                  reduction = "umap", 
                  pt.size = 0.3)
dimplot

ggsave("Sample_dimplot.png", plot=dimplot, bg = "white",width = 10, height = 7)


#Emphasize on the local structure of the data
so <- RunUMAP(so, 
              dims = 1:npcs, 
              verbose = TRUE, 
              seed.use = 42,
              n.neighbors = 10,
              min.dist = 0.3)

so <- FindNeighbors(so, 
                    k.param = 20)
so <- FindClusters(so,
                   resolution = 0.4,
                   algorithm = 3)
so

#------------------------------------------------------------------------------------------------------------------

#At first you have to annotate! Different script


names<- DimPlot(so.annotated, 
                reduction = "umap", 
                pt.size = 0.3, 
                label = FALSE, 
                group.by = "predicted.celltype.l2")
names

ggsave("umap_names.png", plot= names , bg = "white",width = 10, height = 7)


devtools::install_github('immunogenomics/presto')


#Table with frequency of each cell type - Nr of cells in each cell type
annotated_freq <- table(so.annotated@meta.data$predicted.celltype.l2)
annotated_freq
annotated_freq_df <- as.data.frame(annotated_freq)
write_xlsx(annotated_freq_df, "NrofCells_Bycelltype.xlsx")

#Table with number of cells per cluster
cluster_cells <- table(so.annotated@meta.data$seurat_clusters)
cluster_cells
cluster_cells_df <- as.data.frame(cluster_cells)
colnames(cluster_cells_df) <- c("Seurat_Cluster", "Cell_Number")
write_xlsx(cluster_cells_df, "NrofCells_BySeuratCluster.xlsx")
#--------------------------------------------------------------------------------------------------------------------

#Marker genes of clusters 
so.markers <- FindAllMarkers(so.annotated, 
                             only.pos = TRUE, 
                             min.pct = 0.1, 
                             logfc.threshold = 0.25)

so.markers

write_xlsx(so.markers, "so.annotated.markers.xlsx")


#We first have a look at the results, then we specifically look for the top markers of cluster X
#seurat object to tibble table
so.markers %>% 
  as_tibble()
#ten top genes for whichever cluster needed
so.markers %>% 
  as_tibble() %>% 
  dplyr::filter(cluster == "1")

# tοp 20 Cluster x
top20_cluster6 <- so.markers %>% 
  filter(cluster == 6) %>% 
  head(20)

# Excel save
write_xlsx(top20_cluster6, "Cluster6_Top20_Markers.xlsx")


#--------------------------------------------------------------------------------------------------------------------
#PREDICTED CELL TYPE PLOT FOR REFERENCE
p1 <- DimPlot(so.annotated, 
              reduction = "umap", 
              pt.size = 0.3, 
              label = TRUE, 
              repel = TRUE,
              group.by = "predicted.celltype.l2")
p1

ggsave("predicted.celltypes.l2_umap.png", plot= p1, bg = "white",width = 10, height = 7)


#PLOT WHICHEVER GENE NEEDED

p2 <- FeaturePlot(so, 
                  pt.size = 0.3, 
                  features = c("AICDA"))
p2


p3 <- FeaturePlot(so, 
                  pt.size = 0.3, 
                  features = c("IGHA1"))
p3


#plot gois with predicted cell type together
p1 | p2 | p3
ggsave("p1_p2_p3.png", plot= p1|p2|p3, bg = "white",width = 10, height = 7)

p1 | p3
ggsave("p1_p3.png", plot= p1|p3, bg = "white",width = 10, height = 7)


#plot for 2 genes expression patterns
CD5_CXCR4<- FeaturePlot(so, 
                        features = c("CD5", "CXCR4"), 
                        blend = TRUE, 
                        pt.size = 0.3)
CD5_CXCR4
ggsave("CD5_CXCR.png", plot= CD5_CXCR4, bg = "white",width = 10, height = 7)

#-----------------------------------------------------------------------------------------------------

#Choose a gene of your interest and see in which cluster is expressed
so.markers %>% 
  as_tibble() %>% 
  dplyr::filter(gene == "IGHA1")

#goi should be the signature of the cell population of interest
goi <- c("AICDA","DUSP6","EGR1","IGHA2","IGHA1","IGHG3","IGHG2","IGHG1","IGHM","CD3E", "CD4","CD8B","CCR7","SELL","IL7R","LYAR",
         "CCL5","PRF1","GZMA", "GZMK", "GZMB", "NKG7", "PTPN3", "NRGN",
         "LAG3","ENTPD1","PDCD1","HAVCR2","CTLA4", "TIGIT",
         "CXCR4","CD5" ,"CD24", "CD27", "MIR155HG", "CCND2", "TOP2A", "MKI67",
         "PCNA", "MZB1", "XBP1",
         "FOXP3","IL2RA", "IKZF2", 
         "CD69","CD28", "ITGB1", "KLF2", "KLRB1", "TNFRSF4", "CD40LG",
         "HLA-DPB1","ITGAX",
         "MS4A1","CD79A","CD19",
         "ROR1", "LEF1", "PAX5", "BCL2",
         "BIRC3","BTK", "CD86")

#markers are grouped per cluster
#choose what you want to see in the plot
goi <- so.markers %>% 
  as_tibble() %>% 
  group_by(cluster) %>% 
  slice(1) %>% #Pick the top gene for each cluster
  pull(gene) %>%
  unique()

goi <- so.markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5) %>%   # pick the top 5 by log2FC
  pull(gene) %>%
  unique()

#do the plot, grouped by seurat clusters

p <- DotPlot(
  so.annotated,
  features = goi,
  group.by = "seurat_clusters",
  cols = c("blue", "red"),
  dot.scale = 6
) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(
      angle = 45,       # rotation angle
      hjust = 1,        # horizontal justification
      vjust = 1         # vertical justification
    )
  )
p

# OR do the plot, grouped by predicted cell type

p <- DotPlot(
  so.annotated,
  features = goi,
  group.by = "predicted.celltype.l2",
  cols = c("blue", "red"),
  dot.scale = 6
) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(
      angle = 45,       # rotation angle
      hjust = 1,        # horizontal justification
      vjust = 1         # vertical justification
    )
  )
p


#vlnPlot for gois-groupedby predicted.celltype.l2 - WHAT IS GOI NOW??
VlnPlot(so.annotated, 
        features = goi, 
        group.by = "predicted.celltype.l2")

VlnPlot(so.annotated, features = c("IGHA1"), group.by = "predicted.celltype.l2", pt.size = 0) +
  geom_boxplot()

#OR vlnPlot for gois-groupedby seurat.clusters - WHAT IS GOI NOW??
VlnPlot(so.annotated, 
        features = goi, 
        group.by = "seurat_clusters")

VlnPlot(so.annotated, features = c("DUSP6"), group.by = "seurat_clusters", pt.size = 0) +
  geom_boxplot()


#Save the marker genes
saveRDS(so, file = "Sample_with_markers.so.rds")
saveRDS(so.markers, file = "Sample.markers.rds")
saveRDS(so.annotated, file= "Sample.annotated.rds")



#----------------------------------------------------------------------------
#What cell type labels exist in this column?
unique(so.annotated$predicted.celltype.l2)

#This counts how many cells belong to each current Seurat identity class.
table(Idents(so.annotated))

#number of cells in each predicted cell type
Idents(so.annotated) <- "predicted.celltype.l2"
table(Idents(so.annotated))

#---------------------------------------------------------------------------
# ΒρFind the genes that are more expressed in Tregs
#                            ident.1 = "Treg", 
#                            only.pos = TRUE, 
#                            min.pct = 0.25, 
#                            logfc.threshold = 0.25)

#Top 20 genes
#head(treg.markers, n = 20)
#-----------------------------------------------------------------------------

#Add a patient ID to the seurat object of a patient
so@meta.data$patientID <- "Sample"
so.annotated@meta.data$patientID <- "Sample"
saveRDS(so, file = "Sample.so.rds")
saveRDS(so.annotated, file= "Sample.annotated.rds")

#-------------------------------------------------------------------------------

#Create heatmap based on the top x genes of each cluster
goi <- so.markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10) %>%   # pick the top 5 by log2FC
  pull(gene) %>%
  unique()

#Heatmap fot the DEGs
so.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5


markers_heatmap <- DoHeatmap(so, features = goi)
markers_heatmap

ggsave("markers_heatmap.png", plot= markers_heatmap, bg = "white",width = 12, height = 10)



