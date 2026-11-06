library(Signac)
library(Seurat)
library(stringr)
library(dplyr)
library(ggplot2)
library(clustree)
library(reshape2)
library(DoubletFinder)
library(patchwork)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(RColorBrewer)
set.seed(1234)

x = "P365"

SCG = readRDS(file="P365.final.withPETlsi.joint3.rds")

DefaultAssay(SCG) <- "RNA"
SCG <- RunPCA(object = SCG, features = VariableFeatures(SCG), verbose = F)
SCG <- FindNeighbors(object = SCG, dims = 1:30)
SCG <- FindClusters(object = SCG, resolution = 1)
SCG <- RunUMAP(SCG, dims = 1:30, verbose = F)

SCG = readRDS(file="P365.final.withPETlsi.joint3.rds")
# RNA_snn_res.1
# ATAC_snn_res.1
# PET_snn_res.1
# table(SCG@meta.data[,c("RNA_snn_res.1", "bigClass2")])
DefaultAssay(SCG) <- "RNA"
Idents(SCG) = "RNA_snn_res.1"
RNAmarker = FindAllMarkers(SCG, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(RNAmarker, file=paste0(x, ".RNA.marker.genelist.txt"), sep="\t", quote=F)
SCG1 = FindSubCluster(SCG, cluster="0", graph.name = "RNA_snn", subcluster.name = "RNA_snn_sub")

SCG$celltype2 = SCG$Phase
SCG@meta.data[SCG$RNA_snn_res.1==0,]$celltype2 = "ExN" # Cbln3, Gabra6 # CBGRC
SCG@meta.data[SCG$RNA_snn_res.1==1,]$celltype2 = "ExN" # Gabra6 # CBGRC
SCG@meta.data[SCG$RNA_snn_res.1==2,]$celltype2 = "ExN" # Slc17a6, Oxt, Baiap3
SCG@meta.data[SCG$RNA_snn_res.1==3,]$celltype2 = "Oligo" # Olig1
SCG@meta.data[SCG$RNA_snn_res.1==4,]$celltype2 = "ExN" # Slc17a7, Lamp5
SCG@meta.data[SCG$RNA_snn_res.1==5,]$celltype2 = "Ast" # Slc6a11
SCG@meta.data[SCG$RNA_snn_res.1==6,]$celltype2 = "InN" # Gad2, Syt6, Gad1(low),
SCG@meta.data[SCG$RNA_snn_res.1==7,]$celltype2 = "ExN" # Rprm, Slc17a7
SCG@meta.data[SCG$RNA_snn_res.1==8,]$celltype2 = "ExN" # C1ql2, Nr3c2
SCG@meta.data[SCG$RNA_snn_res.1==9,]$celltype2 = "InN" # Gad1
SCG@meta.data[SCG$RNA_snn_res.1==10,]$celltype2 = "InN" # Gad1, Gad2, Sst
SCG@meta.data[SCG$RNA_snn_res.1==11,]$celltype2 = "OPC" # Pdgfra, Olig1, Sox10
SCG@meta.data[SCG$RNA_snn_res.1==12,]$celltype2 = "InN" # Pvalb, Gad1
SCG@meta.data[SCG$RNA_snn_res.1==13,]$celltype2 = "Microglia" # Cx3cr1

# pdf(paste0(x, ".RNA.gene.exp.pdf"), width=20, height=20)
# genelist = c("Pdgfra", "Olig1", "Aqp4","Hes5","Rorb", "Ctss", "Ptprc", "Apbb1ip", "Cldn5", "Ptgds", "Slc6a13", "Vtn", "Flt1", "Mgst1", "Frzb", "Meg3", "Slc17a7", "Gad2", "Gad1", "Slc17a6", "Cbln3", "Gabra6", "Lamp5", "Cpne7", "Fam107a", "Syt6", "Nr3c2", "C1ql2", "Sst", "Npy", "Vip", "Lhx6", "Pvalb")
# FeaturePlot(SCG, reduction ="tsne.rna", features = genelist[1:4], pt.size=0.1)
# FeaturePlot(SCG, reduction ="tsne.rna", features = genelist[5:8], pt.size=0.1)
# FeaturePlot(SCG, reduction ="tsne.rna", features = genelist[9:12], pt.size=0.1)
# FeaturePlot(SCG, reduction ="tsne.rna", features = genelist[13:16], pt.size=0.1)
# FeaturePlot(SCG, reduction ="tsne.rna", features = genelist[17:20], pt.size=0.1)
# FeaturePlot(SCG, reduction ="tsne.rna", features = genelist[21:24], pt.size=0.1)
# FeaturePlot(SCG, reduction ="tsne.rna", features = genelist[25:28], pt.size=0.1)
# FeaturePlot(SCG, reduction ="tsne.rna", features = genelist[29:33], pt.size=0.1)
# dev.off()

pdf(paste0(x, ".multiomics.result.withanno.3.pdf"), width=15)
p1 = DimPlot(SCG,reduction = "tsne.rna", group.by = "celltype2", pt.size = 0.5) + ggtitle("single RNA TSNE")
p2 = DimPlot(SCG,reduction = "umap.rna", group.by = "celltype2", pt.size = 0.5) + ggtitle("single RNA UMAP")
print(p1+p2)
p1 = DimPlot(SCG,reduction = "tsne.atac", group.by = "celltype2", pt.size = 0.5) + ggtitle("single ATAC TSNE")
p2 = DimPlot(SCG,reduction = "umap.atac", group.by = "celltype2", pt.size = 0.5) + ggtitle("single ATAC UMAP")
print(p1+p2)
p1 = DimPlot(SCG,reduction = "tsne.lsi.pet", group.by = "celltype2", pt.size = 0.5) + ggtitle("single PET TSNE")
p2 = DimPlot(SCG,reduction = "umap.lsi.pet", group.by = "celltype2", pt.size = 0.5) + ggtitle("single PET UMAP")
print(p1+p2)
# p1 = DimPlot(SCG,reduction = "tsne.pca.pet", group.by = "celltype2", pt.size = 0.5) + ggtitle("single PET TSNE")
# p2 = DimPlot(SCG,reduction = "umap.pca.pet", group.by = "celltype2", pt.size = 0.5) + ggtitle("single PET UMAP")
# print(p1+p2)
p1 <- DimPlot(SCG, reduction = "wnn3.umap", group.by = "celltype2", pt.size = 0.5) + ggtitle("UMAPWNN3")
p2 <- DimPlot(SCG, reduction = "wnn3.umap", group.by = "celltype2", pt.size = 0.5) + ggtitle("UMAPWNN3")
print(p1+p2)
dev.off()

# SCGtmp = subset(SCG, subset=(RNA_snn_res.1==0 | RNA_snn_res.1==1 | RNA_snn_res.1==8))
# pdf(paste0(x, ".multiomics.result.withanno.4.pdf"), width=15)
# p1 = DimPlot(SCGtmp,reduction = "tsne.rna", group.by = "RNA_snn_res.1", pt.size = 0.5) + ggtitle("single RNA TSNE")
# p2 = DimPlot(SCGtmp,reduction = "umap.rna", group.by = "RNA_snn_res.1", pt.size = 0.5) + ggtitle("single RNA UMAP")
# print(p1+p2)
# dev.off()
saveRDS(SCG, file="P365.final.withPETlsi.joint3.withcustomanno.rds")
write.table(SCG@meta.data, file="P365.final.withPETlsi.joint3.withcustomanno.metadata.txt", sep="\t", quote=F)

celltypeorder = c(0,1,2,4,7,8,6,9,10,12,5,13,3,11)
topgene = c("Slc17a7", "Slc17a6", "Gabra6", "Oxt", "Rprm", "C1ql2", "Gad2", "Gad1", "Syt6", "Sst", "Pvalb", "Slc6a11", "Cx3cr1", "Olig1", "Mbp", "Pdgfra")
SCG$RNA_snn_res.1 = factor(SCG$RNA_snn_res.1,levels = celltypeorder)

pdf(paste0(x, ".plotmarkgene.3.pdf"), width=15)
DotPlot(SCG, features = topgene, group.by = "RNA_snn_res.1") + scale_colour_gradientn(colours = colorRampPalette(c("white", "#c9c9e6", "#67001F"))(100)) + RotatedAxis() + coord_flip()
DotPlot(SCG, features = topgene, group.by = "RNA_snn_res.1") + scale_colour_gradientn(colours = colorRampPalette(c("white", "#c9c9e6", "red", "red"))(100)) + RotatedAxis() + coord_flip()
dev.off()



SCG = readRDS(file=paste0("P365", ".final.withPETlsi.joint3.withcustomanno.rds"))
DefaultAssay(SCG) <- "ATAC"
bb =  GeneActivity(SCG)
aa = bb
SCG[['GeneACT']] <- CreateAssayObject(counts = aa)
DefaultAssay(SCG) <- "GeneACT"
SCG = NormalizeData(SCG)
SCG <- ScaleData(object = SCG, features=rownames(SCG))
celltypeorder = c(0,1,2,4,7,8,6,9,10,12,5,13,3,11)
topgene = c("Slc17a7", "Slc17a6", "Gabra6", "Oxt", "Rprm", "C1ql2", "Gad2", "Gad1", "Syt6", "Sst", "Pvalb", "Slc6a11", "Cx3cr1", "Olig1", "Mbp", "Pdgfra")
SCG$RNA_snn_res.1 = factor(SCG$RNA_snn_res.1,levels = celltypeorder)

pdf(paste0(x, ".plotmarkgene.inACT.3.pdf"), width=15)
DotPlot(SCG, features = topgene, group.by = "RNA_snn_res.1") + scale_colour_gradientn(colours = colorRampPalette(c("white", "#c9c9e6", "#67001F"))(100)) + RotatedAxis() + coord_flip()
DotPlot(SCG, features = topgene, group.by = "RNA_snn_res.1") + scale_colour_gradientn(colours = colorRampPalette(c("white", "#c9c9e6", "red", "red"))(100)) + RotatedAxis() + coord_flip()
dev.off()
write(t(colnames(SCG@assays$RNA@counts)), file="barcode.list")



celltypeorder = c(0,1,2,4,7,8,6,9,10,12,5,13,3,11)
topgene = c("Slc17a7", "Slc17a6", "Gabra6", "Oxt", "Rprm", "C1ql2", "Gad2", "Gad1", "Syt6", "Sst", "Pvalb", "Slc6a11", "Cx3cr1", "Olig1", "Mbp", "Pdgfra")
SCG$RNA_snn_res.1 = factor(SCG$RNA_snn_res.1,levels = celltypeorder)

pdf(paste0(x, ".plotmarkgene.inACT.3.pdf"), width=15)
promoterbodyex2000fragmentmat = read.table("mouseBrainP365.fragment.gene.promoterbodyex2000.info.bed.mat", check.names = F)
promoterbodyex2000fragmentmat$genename = rownames(promoterbodyex2000fragmentmat)

geneinfo = read.table("~/Tools/cellranger-7.1.0/genome/refdata-gex-mm10-2020-A/genes/mm10.info", check.names = F)
names(geneinfo) = c("chrom", "start", "end", "geneid", "score", "strand", "genename", "genetype")
geneinfo$length = geneinfo$end - geneinfo$start
promoterbodyex2000fragmentmatall = left_join(promoterbodyex2000fragmentmat, geneinfo, by="genename")
rownames(promoterbodyex2000fragmentmatall) = promoterbodyex2000fragmentmatall$genename
promoterbodyex2000fragmentmatall$fc0 <- promoterbodyex2000fragmentmatall$length/5000

pbf = Matrix::Matrix(as.matrix(promoterbodyex2000fragmentmatall[,1:5935]), sparse = T) # promoterbodyfragment
pbffc0 = Matrix::Matrix(pbf / promoterbodyex2000fragmentmatall$fc0, sparse = T)

SCG[['pbf']] <- CreateAssayObject(counts = pbf)
DefaultAssay(SCG) <- "pbf"
SCG = NormalizeData(SCG)
SCG <- ScaleData(object = SCG, features=rownames(SCG))
DotPlot(SCG, features = topgene, group.by = "RNA_snn_res.1") + scale_colour_gradientn(colours = colorRampPalette(c("white", "#c9c9e6", "red"))(100)) + RotatedAxis() + coord_flip() + ggtitle("pbf")

SCG[['pbffc0']] <- CreateAssayObject(counts = pbffc0)
DefaultAssay(SCG) <- "pbffc0"
SCG = NormalizeData(SCG)
SCG <- ScaleData(object = SCG, features=rownames(SCG))
DotPlot(SCG, features = topgene, group.by = "RNA_snn_res.1") + scale_colour_gradientn(colours = colorRampPalette(c("white", "#c9c9e6", "red"))(100)) + RotatedAxis() + coord_flip() + ggtitle("pbffc0")


promoterbodyex2000PETmat = read.table("mouseBrainP365.PET.gene.promoterbodyex2000.info.bed.mat", check.names = F)
promoterbodyex2000PETmat$genename = rownames(promoterbodyex2000PETmat)

geneinfo = read.table("~/Tools/cellranger-7.1.0/genome/refdata-gex-mm10-2020-A/genes/mm10.info", check.names = F)
names(geneinfo) = c("chrom", "start", "end", "geneid", "score", "strand", "genename", "genetype")
geneinfo$length = geneinfo$end - geneinfo$start
promoterbodyex2000PETmatall = left_join(promoterbodyex2000PETmat, geneinfo, by="genename")
rownames(promoterbodyex2000PETmatall) = promoterbodyex2000PETmatall$genename
promoterbodyex2000PETmatall$fc0 <- promoterbodyex2000PETmatall$length/5000

pbp = Matrix::Matrix(as.matrix(promoterbodyex2000PETmatall[,1:5935]), sparse = T) # promoterbodyPET
pbpfc0 = Matrix::Matrix(pbp / promoterbodyex2000PETmatall$fc0, sparse = T)

SCG[['pbp']] <- CreateAssayObject(counts = pbp)
DefaultAssay(SCG) <- "pbp"
SCG = NormalizeData(SCG)
SCG <- ScaleData(object = SCG, features=rownames(SCG))
DotPlot(SCG, features = topgene, group.by = "RNA_snn_res.1") + scale_colour_gradientn(colours = colorRampPalette(c("white", "#c9c9e6", "red"))(100)) + RotatedAxis() + coord_flip() + ggtitle("pbp")

SCG[['pbpfc0']] <- CreateAssayObject(counts = pbpfc0)
DefaultAssay(SCG) <- "pbpfc0"
SCG = NormalizeData(SCG)
SCG <- ScaleData(object = SCG, features=rownames(SCG))
DotPlot(SCG, features = topgene, group.by = "RNA_snn_res.1") + scale_colour_gradientn(colours = colorRampPalette(c("white", "#c9c9e6", "red"))(100)) + RotatedAxis() + coord_flip() + ggtitle("pbpfc0")
dev.off()

saveRDS(SCG, file=paste0("P365", ".final.withPETlsi.joint3.withcustomanno.final.rds"))

SCG = readRDS(file=paste0("P365", ".final.withPETlsi.joint3.withcustomanno.final.rds")) # update 2026.10.25

celltypeorder = c("ExN","InN","Ast","Microglia","Oligo","OPC")
topgene = c("Slc17a7", "Slc17a6", "Gabra6", "Oxt", "Rprm", "C1ql2", "Gad2", "Gad1", "Syt6", "Sst", "Pvalb", "Slc6a11", "Cx3cr1", "Olig1", "Mbp", "Pdgfra")
SCG$celltype2 = factor(SCG$celltype2, levels = celltypeorder)

pdf(paste0(x, ".plotmarkgenemarkpeakmarkloop.final.pdf"), width=15)
DefaultAssay(SCG) <- "RNA"
DotPlot(SCG, features = topgene, group.by = "celltype2") + scale_colour_gradientn(colours = colorRampPalette(c("white", "#c9c9e6", "red"))(100)) + RotatedAxis() + coord_flip() + ggtitle("RNA")

# DefaultAssay(SCG) <- "pbf"
# DotPlot(SCG, features = topgene, group.by = "celltype2") + scale_colour_gradientn(colours = colorRampPalette(c("white", "#c9c9e6", "red"))(100)) + RotatedAxis() + coord_flip() + ggtitle("pbf")
DefaultAssay(SCG) <- "pbf" # update 2026.10.25
DotPlot(SCG, features = topgene, group.by = "celltype2") + scale_colour_gradientn(colours = colorRampPalette(c("white", "#c9c9e6", "red"))(100)) + RotatedAxis() + coord_flip() + ggtitle("pbf") # update 2026.10.25

# DefaultAssay(SCG) <- "pbp"
# DotPlot(SCG, features = topgene, group.by = "celltype2") + scale_colour_gradientn(colours = colorRampPalette(c("white", "#c9c9e6", "red"))(100)) + RotatedAxis() + coord_flip() + ggtitle("pbp")
DefaultAssay(SCG) <- "pbp" # update 2026.10.25
DotPlot(SCG, features = topgene, group.by = "celltype2") + scale_colour_gradientn(colours = colorRampPalette(c("white", "#c9c9e6", "red"))(100)) + RotatedAxis() + coord_flip() + ggtitle("pbp") # update 2026.10.25
dev.off()

DefaultAssay(SCG) <- "RNA"
genelist = VariableFeatures(SCG)
pdf(paste0(x, ".plotmarkgenemarkpeakmarkloop.plotloop.pdf"), width=15)
DefaultAssay(SCG) <- "pbp"
SCG <- FindVariableFeatures(object = SCG, nfeatures = 3000, verbose=F)
SCG <- RunPCA(object = SCG, reduction.name="pca.pbp", reduction.key = "pbpPC_", features = genelist)
SCG <- RunTSNE(SCG, reduction = "pca.pbp", dims = 1:15, verbose = F, reduction.name="tsne.pbp", reduction.key = "pbpTSNE_", check_duplicates = FALSE)
SCG <- RunUMAP(SCG, reduction = "pca.pbp", dims = 1:15, verbose = F, reduction.name="umap.pbp", reduction.key = "pbpUMAP_")
p1 = DimPlot(SCG,reduction = "tsne.pbp", group.by = "celltype2", pt.size = 0.5, cols = c("#b79f00", "#00b938", "#f8766c", "#00bec4", "#609cff", "#f563e2")) + ggtitle("pbp loop TSNE")
p2 = DimPlot(SCG,reduction = "umap.pbp", group.by = "celltype2", pt.size = 0.5, cols = c("#b79f00", "#00b938", "#f8766c", "#00bec4", "#609cff", "#f563e2")) + ggtitle("pbp loop UMAP")
print(p1+p2)

DefaultAssay(SCG) <- "pbf"
SCG <- FindVariableFeatures(object = SCG, nfeatures = 3000, verbose=F)
SCG <- RunPCA(object = SCG, reduction.name="pca.pbf", reduction.key = "pbfPC_")
# SCG <- RunPCA(object = SCG, reduction.name="pca.pbf", reduction.key = "pbfPC_", features = genelist)
SCG <- RunTSNE(SCG, reduction = "pca.pbf", dims = 1:15, verbose = F, reduction.name="tsne.pbf", reduction.key = "pbfTSNE_", check_duplicates = FALSE)
SCG <- RunUMAP(SCG, reduction = "pca.pbf", dims = 1:15, verbose = F, reduction.name="umap.pbf", reduction.key = "pbfUMAP_")
p1 = DimPlot(SCG,reduction = "tsne.pbf", group.by = "celltype2", pt.size = 0.5, cols = c("#b79f00", "#00b938", "#f8766c", "#00bec4", "#609cff", "#f563e2")) + ggtitle("pbf peak TSNE")
p2 = DimPlot(SCG,reduction = "umap.pbf", group.by = "celltype2", pt.size = 0.5, cols = c("#b79f00", "#00b938", "#f8766c", "#00bec4", "#609cff", "#f563e2")) + ggtitle("pbf peak UMAP")
print(p1+p2)

SCG <- FindMultiModalNeighbors(SCG, reduction.list = list("pca.rna", "pca.pbf", 'pca.pbp'), dims.list = list(1:30, 1:20, 1:20), knn.graph.name = "wknn2", snn.graph.name = "wsnn2", weighted.nn.name = "weighted.nn2")
SCG <- RunUMAP(SCG, nn.name = "weighted.nn2", reduction.name = "wnn3.umap2", reduction.key = "wnn3UMAP2_")
p1 <- DimPlot(SCG, reduction = "wnn3.umap2", group.by = "celltype2", pt.size = 0.5, cols = c("#b79f00", "#00b938", "#f8766c", "#00bec4", "#609cff", "#f563e2")) + ggtitle("UMAPWNN32")
print(p1+p1)
dev.off()


saveRDS(SCG, file=paste0("P365", ".final.withPETlsi.joint3.withcustomanno.final.addepi.rds"))
