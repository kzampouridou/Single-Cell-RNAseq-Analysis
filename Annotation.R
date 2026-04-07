library(Seurat)
library(ggplot2)
library(patchwork)
library(SeuratData)

#Optional load different datasets using SeuratData package
#options(SeuratData.repo.use = "http://seurat.nygenome.org")
#InstallData("pbmc3k")
#pbmc3k <- LoadData('pbmc3k')
#pbmc3k <- UpdateSeuratObject(pbmc3k)
#pbmc3k <- SCTransform(pbmc3k, verbose = FALSE)

#Load the reference dataset from Hao et al. 2021 (downloaded from Zenodo)
reference <- readRDS("C:/pbmc_multimodal_2023.rds")
DimPlot(object = reference, reduction = "wnn.umap", group.by = "celltype.l2", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()

#Find anchors between the query and the reference dataset
anchors <- FindTransferAnchors(
  reference = reference,
  query = so,
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:50
)

#Transfer cell type labels and protein data from the reference to the query
so.annotated <- MapQuery(
  anchorset = anchors,
  query = so,
  reference = reference,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "spca", 
  reduction.model = "wnn.umap"
)

#Vizualise the cells
p2 <-  DimPlot(so.annotated, reduction = "ref.umap", 
             group.by = "predicted.celltype.l2", 
             label = TRUE, 
             label.size = 3 ,
             repel = TRUE) + NoLegend()
p2

ggsave("predicted_cell_types_annotated.png", plot= p2, bg = "white",width = 10, height = 7)

saveRDS(so.annotated, file= "data/sample.annotated.rds")

