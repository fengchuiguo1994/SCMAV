library(Signac)
library(Seurat)
library(stringr)
library(dplyr)
library(ggplot2)
library(clustree)
library(reshape2)
library(DoubletFinder)
library(patchwork)
library(SoupX)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
args=commandArgs(T) # 

args = c("P365.prepareData.rds", "P365")
inputrds = args[1]
x = args[2]
genelistrds = "mouse_cell_cycle_genes.rds" # or "" # https://www.r-bloggers.com/2016/10/converting-mouse-to-human-gene-names-with-biomart-package/
set.seed(1234)

SCG = readRDS(file = inputrds)

# cell cycle
DefaultAssay(SCG) <- "RNA"
SCG <- NormalizeData(object = SCG, verbose=F)
SCG <- ScaleData(object = SCG, features=rownames(SCG))
if(genelistrds == "") {
  s.genes = intersect(rownames(SCG@assays$RNA@counts), cc.genes.updated.2019$s.genes)
  g2m.genes = intersect(rownames(SCG@assays$RNA@counts), cc.genes.updated.2019$g2m.genes)
} else {
  genelist = readRDS(genelistrds) 
  s.genes <- intersect(rownames(SCG@assays$RNA@counts), genelist$s.genes)
  g2m.genes <- intersect(rownames(SCG@assays$RNA@counts), genelist$g2m.genes)
}
cycleGenes <- c(s.genes, g2m.genes)
cycleGenesMat <- GetAssayData(SCG, slot = "counts")[cycleGenes, ]
SCG$cycleGenesCount = colSums(cycleGenesMat)
SCG$cycleGenesCountRatio = SCG$cycleGenesCount/SCG$nCount_RNA

SCG <- CellCycleScoring(SCG, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)
SCG <- RunPCA(object = SCG, features = c(s.genes, g2m.genes), reduction.name="pca.cellcycle", reduction.key = "CCPC_")
SCG <- RunTSNE(SCG, reduction = "pca.cellcycle", dims = 1:15, verbose = F, reduction.name="tsne.cellcycle", reduction.key = "CCTSNE_", check_duplicates = FALSE)
SCG <- RunUMAP(SCG, reduction = "pca.cellcycle", dims = 1:15, verbose = F, reduction.name="umap.cellcycle", reduction.key = "CCUMAP_")
pdf(paste0(x, ".cellcycle.pdf"), width=15)
p1 = DimPlot(SCG,reduction = "pca.cellcycle", group.by = "Phase", pt.size = 1) + ggtitle("cellcyclePCA")
p2 = DimPlot(SCG,reduction = "tsne.cellcycle", group.by = "Phase", pt.size = 1) + ggtitle("cellcycleTSNE")
p3 = DimPlot(SCG,reduction = "umap.cellcycle", group.by = "Phase", pt.size = 1) + ggtitle("cellcycleUMAP")
p4 = FeaturePlot(SCG, features="cycleGenesCountRatio", reduction = "pca.cellcycle", pt.size = 1) + ggtitle("cellcyclePCA")
p5 = FeaturePlot(SCG, features="cycleGenesCountRatio", reduction = "tsne.cellcycle", pt.size = 1) + ggtitle("cellcycleTSNE")
p6 = FeaturePlot(SCG, features="cycleGenesCountRatio", reduction = "umap.cellcycle", pt.size = 1) + ggtitle("cellcycleUMAP")
print(p1+p4)
print(p2+p5)
print(p3+p6)
dev.off()
saveRDS(SCG, file=paste0(x, ".final.cellcycle.rds"))

# RNA cluster
DefaultAssay(SCG) <- "RNA"
SCG <- NormalizeData(object = SCG, verbose=F)
SCG <- FindVariableFeatures(object = SCG, nfeatures = 3000, verbose=F)
SCG <- ScaleData(object = SCG, features=rownames(SCG))
pdf(paste0(x, ".single.RNA.result.pdf"))
SCG <- RunPCA(object = SCG, features = VariableFeatures(SCG), reduction.name="pca.rna", reduction.key = "RNAPC_")
SCG@reductions$pca = SCG@reductions$pca.rna
SCG <- JackStraw(SCG, dims = 50, num.replicate =100)
SCG <- ScoreJackStraw(SCG, dims = 1:50)
p1 = JackStrawPlot(SCG, dims =1:50, ymax = 0.7, xmax = 0.2) + NoLegend()
print(p1)
pct = SCG[["pca"]]@stdev / sum(SCG[["pca"]]@stdev) * 100 # pca 85%
cumu = cumsum(pct)
per85=which(cumu > 85 & pct < 5)[1]
dif01=sort(which(pct[1:length(pct)-1] - pct[2:length(pct)] > 0.1), decreasing = T)[1] + 1
p2 = ElbowPlot(SCG, ndims = 50)$data %>% ggplot() +
    geom_point(aes(x=dims, y=stdev)) +
    geom_vline(xintercept = per85, color="darkred") +
    geom_vline(xintercept = dif01, color="darkgreen") +
    theme_bw() + labs(title="Elbow plot: quantitative approach")
print(p2)
SCG@reductions$pca = NULL
SCG <- FindNeighbors(object = SCG, reduction = "pca.rna", dims = 1:30)
SCG <- FindClusters(object = SCG, resolution = 1)
SCG <- RunTSNE(SCG, reduction = "pca.rna", dims = 1:30, verbose = F, reduction.name="tsne.rna", reduction.key = "RNATSNE_", check_duplicates = FALSE)
SCG <- RunUMAP(SCG, reduction = "pca.rna", dims = 1:30, verbose = F, reduction.name="umap.rna", reduction.key = "RNAUMAP_")
p1 = DimPlot(SCG,reduction = "tsne.rna", group.by = "RNA_snn_res.1", pt.size = 1) + ggtitle("single RNA TSNE")
p2 = DimPlot(SCG,reduction = "umap.rna", group.by = "RNA_snn_res.1", pt.size = 1) + ggtitle("single RNA UMAP")
print(p1)
print(p2)
dev.off()


# ATAC cluster
DefaultAssay(SCG) <- "ATAC"
SCG <- RunTFIDF(SCG, assay="ATAC")
SCG <- FindTopFeatures(SCG, assay="ATAC", min.cutoff = "q5")
SCG <- RunSVD(SCG, assay="ATAC", reduction.key = "ATACLSI_", reduction.name = "lsi.atac")
pdf(paste0(x, ".single.ATAC.result.pdf"))
DepthCor(SCG, reduction = "lsi.atac", n = 15)
SCG@reductions$pca = SCG@reductions$lsi.atac
pct = SCG[["pca"]]@stdev / sum(SCG[["pca"]]@stdev) * 100
cumu = cumsum(pct)
per85=which(cumu > 85 & pct < 5)[1]
dif01=sort(which(pct[1:length(pct)-1] - pct[2:length(pct)] > 0.1), decreasing = T)[1] + 1
p2 = ElbowPlot(SCG, ndims = 50)$data %>% ggplot() +
    geom_point(aes(x=dims, y=stdev)) +
    geom_vline(xintercept = per85, color="darkred") +
    geom_vline(xintercept = dif01, color="darkgreen") +
    theme_bw() + labs(title="Elbow plot: quantitative approach")
print(p2)
SCG@reductions$pca = NULL
SCG <- FindNeighbors(SCG, reduction = 'lsi.atac', dims = 2:20, k.param = 30)
SCG <- FindClusters(SCG, resolution = 1, algorithm = 3)
SCG <- RunTSNE(SCG, reduction = "lsi.atac", dims = 2:20, reduction.name = "tsne.atac", reduction.key = "ATACTSNE_", check_duplicates = FALSE)
SCG <- RunUMAP(SCG, reduction = "lsi.atac", dims = 2:20, reduction.name = "umap.atac", reduction.key = "ATACUMAP_")
p1 = DimPlot(SCG,reduction = "tsne.atac", group.by = "ATAC_snn_res.1", pt.size = 1) + ggtitle("single ATAC TSNE")
p2 = DimPlot(SCG,reduction = "umap.atac", group.by = "ATAC_snn_res.1", pt.size = 1) + ggtitle("single ATAC UMAP")
print(p1)
print(p2)
p1 = DimPlot(SCG,reduction = "tsne.rna", group.by = "ATAC_snn_res.1", pt.size = 1) + ggtitle("single RNA TSNE with ATAC cluster")
p2 = DimPlot(SCG,reduction = "umap.rna", group.by = "ATAC_snn_res.1", pt.size = 1) + ggtitle("single RNA UMAP with ATAC cluster")
print(p1)
print(p2)
p1 = DimPlot(SCG,reduction = "tsne.atac", group.by = "RNA_snn_res.1", pt.size = 1) + ggtitle("single ATAC TSNE with RNA cluster")
p2 = DimPlot(SCG,reduction = "umap.atac", group.by = "RNA_snn_res.1", pt.size = 1) + ggtitle("single ATAC UMAP with RNA cluster")
print(p1)
print(p2)
dev.off()
# saveRDS(SCG, file="tmp.rds")

# PET cluster
DefaultAssay(SCG) <- "PET"
SCG <- FindTopFeatures(SCG, assay="PET", min.cutoff = "q5")
SCG = RunSVD(SCG, assay = "PET", reduction.key = "PETLSI_", reduction.name = "lsi.pet")
pdf(paste0(x, ".single.PET.lsi.result.pdf"))
DepthCor(SCG, reduction = "lsi.pet", n = 15)
SCG@reductions$pca = SCG@reductions$lsi.pet
pct = SCG[["pca"]]@stdev / sum(SCG[["pca"]]@stdev) * 100
cumu = cumsum(pct)
per85=which(cumu > 85 & pct < 5)[1]
dif01=sort(which(pct[1:length(pct)-1] - pct[2:length(pct)] > 0.1), decreasing = T)[1] + 1
p2 = ElbowPlot(SCG, ndims = 50)$data %>% ggplot() +
    geom_point(aes(x=dims, y=stdev)) +
    geom_vline(xintercept = per85, color="darkred") +
    geom_vline(xintercept = dif01, color="darkgreen") +
    theme_bw() + labs(title="Elbow plot: quantitative approach")
print(p2)
SCG@reductions$pca = NULL
SCG <- FindNeighbors(SCG, reduction = 'lsi.pet', dims = 1:20, k.param = 30)
SCG <- FindClusters(SCG, resolution = 1, algorithm = 3)
SCG <- RunTSNE(SCG, reduction = "lsi.pet", dims = 1:20, reduction.name = "tsne.lsi.pet", reduction.key = "PETLSITSNE_", check_duplicates = FALSE)
SCG <- RunUMAP(SCG, reduction = "lsi.pet", dims = 1:20, reduction.name = "umap.lsi.pet", reduction.key = "PETLSIUMAP_")
p1 = DimPlot(SCG,reduction = "tsne.lsi.pet", group.by = "PET_snn_res.1", pt.size = 1) + ggtitle("single PET TSNE")
p2 = DimPlot(SCG,reduction = "umap.lsi.pet", group.by = "PET_snn_res.1", pt.size = 1) + ggtitle("single PET UMAP")
print(p1)
print(p2)
p1 = DimPlot(SCG,reduction = "tsne.rna", group.by = "PET_snn_res.1", pt.size = 1) + ggtitle("single RNA TSNE with PET cluster")
p2 = DimPlot(SCG,reduction = "umap.rna", group.by = "PET_snn_res.1", pt.size = 1) + ggtitle("single RNA UMAP with PET cluster")
print(p1)
print(p2)
p1 = DimPlot(SCG,reduction = "tsne.lsi.pet", group.by = "RNA_snn_res.1", pt.size = 1) + ggtitle("single PET TSNE with RNA cluster")
p2 = DimPlot(SCG,reduction = "umap.lsi.pet", group.by = "RNA_snn_res.1", pt.size = 1) + ggtitle("single PET UMAP with RNA cluster")
print(p1)
print(p2)
dev.off()
saveRDS(SCG, file=paste0(x, ".final.withPETlsi.rds"))

###########################################################################
#         can use pca for PET data, but it's very slow                    #
###########################################################################
# SCG <- FindVariableFeatures(object = SCG, nfeatures = 10000, verbose=F)
# SCG <- ScaleData(object = SCG, features=rownames(SCG))
# SCG <- RunPCA(object = SCG, features = VariableFeatures(SCG), reduction.name="pca.pet", reduction.key = "PETPC_")
# pdf(paste0(x, ".single.PET.pca.result.pdf"))
# SCG@reductions$pca = SCG@reductions$pca.pet
# SCG <- JackStraw(SCG, dims = 50, num.replicate =100)
# SCG <- ScoreJackStraw(SCG, dims = 1:50)
# p1 = JackStrawPlot(SCG, dims =1:50, ymax = 0.7, xmax = 0.2) + NoLegend()
# print(p1)
# pct = SCG[["pca"]]@stdev / sum(SCG[["pca"]]@stdev) * 100 # pca 85%
# cumu = cumsum(pct)
# per85=which(cumu > 85 & pct < 5)[1]
# dif01=sort(which(pct[1:length(pct)-1] - pct[2:length(pct)] > 0.1), decreasing = T)[1] + 1
# p2 = ElbowPlot(SCG, ndims = 50)$data %>% ggplot() +
#     geom_point(aes(x=dims, y=stdev)) +
#     geom_vline(xintercept = per85, color="darkred") +
#     geom_vline(xintercept = dif01, color="darkgreen") +
#     theme_bw() + labs(title="Elbow plot: quantitative approach")
# print(p2)
# SCG@reductions$pca = NULL
# SCG <- FindNeighbors(object = SCG, reduction = "pca.pet", dims = 1:20)
# SCG <- FindClusters(object = SCG, resolution = 1)
# SCG <- RunTSNE(SCG, reduction = "pca.pet", dims = 1:20, verbose = F, reduction.name="tsne.pca.pet", reduction.key = "PETPCATSNE_", check_duplicates = FALSE)
# SCG <- RunUMAP(SCG, reduction = "pca.pet", dims = 1:20, verbose = F, reduction.name="umap.pca.pet", reduction.key = "PETPCAUMAP_")
# p1 = DimPlot(SCG,reduction = "tsne.pca.pet", group.by = "PET_snn_res.1", pt.size = 1) + ggtitle("single PET TSNE")
# p2 = DimPlot(SCG,reduction = "umap.pca.pet", group.by = "PET_snn_res.1", pt.size = 1) + ggtitle("single PET UMAP")
# print(p1)
# print(p2)
# p1 = DimPlot(SCG,reduction = "tsne.rna", group.by = "PET_snn_res.1", pt.size = 1) + ggtitle("single RNA TSNE with PET cluster")
# p2 = DimPlot(SCG,reduction = "umap.rna", group.by = "PET_snn_res.1", pt.size = 1) + ggtitle("single RNA UMAP with PET cluster")
# print(p1)
# print(p2)
# p1 = DimPlot(SCG,reduction = "tsne.pca.pet", group.by = "RNA_snn_res.1", pt.size = 1) + ggtitle("single PET TSNE with RNA cluster")
# p2 = DimPlot(SCG,reduction = "umap.pca.pet", group.by = "RNA_snn_res.1", pt.size = 1) + ggtitle("single PET UMAP with RNA cluster")
# print(p1)
# print(p2)
# dev.off()

# WNN analysis
# SCG <- FindMultiModalNeighbors(SCG, reduction.list = list("pca.rna", "lsi.atac", 'pca.pet'), dims.list = list(1:30, 2:20, 1:20))
SCG <- FindMultiModalNeighbors(SCG, reduction.list = list("pca.rna", "lsi.atac", 'lsi.pet'), dims.list = list(1:30, 2:20, 1:20))
SCG <- FindClusters(SCG, graph.name = "wsnn", algorithm = 3, resolution = 1)
SCG <- RunUMAP(SCG, nn.name = "weighted.nn", reduction.name = "wnn3.umap", reduction.key = "wnn3UMAP_")

pdf(paste0(x, ".multiomics.result.pdf"))
p1 <- DimPlot(SCG, reduction = "wnn3.umap", group.by = "wsnn_res.1", pt.size = 1) + ggtitle("UMAPWNN3")
print(p1)
p1 <- DimPlot(SCG, reduction = "wnn3.umap", group.by = "RNA_snn_res.1", pt.size = 1) + ggtitle("UMAPWNN3 with RNA cluster")
print(p1)
dev.off()

saveRDS(SCG, file=paste0(x, ".final.withPETlsi.joint3.rds"))