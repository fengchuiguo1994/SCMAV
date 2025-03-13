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
args=commandArgs(T) # args = c("/data/home/ruanlab/huangxingyu/scChAIR/Pastki/result/SCG0089-15872_SI-NA-D6/10xGEX/outs/raw_feature_bc_matrix","/data/home/ruanlab/huangxingyu/scChAIR/Pastki/result/SCG0089-15872_SI-NA-D6/10xATAC/outs/raw_peak_bc_matrix.h5","/data/home/ruanlab/huangxingyu/scChAIR/Pastki/result/SCG0089-15872_SI-NA-D6/bulkChiatac/atac_fragments.tsv.gz","/data/home/ruanlab/huangxingyu/scChAIR/Pastki/result/SCG0089-15872_SI-NA-D6/SCG0089.1000k.bedpe.matrix.rds","/data/home/ruanlab/huangxingyu/scChAIR/Pastki/result/SCG0089-15872_SI-NA-D6/bulkChiatac/SCG0089-15872_SI-NA-D6.hq.CB.nr.bedpe.stat",200,500,50,"SCG0089","/data/home/ruanlab/huangxingyu/scChAIR/Pastki/result/SCG0089-15872_SI-NA-D6")
# RNA ATAC ATACfragment PET PETinfo RNAcutoff ATACcutoff PETcutoff libraryID outdir

set.seed(1234)

x = args[9]
print(x)
nRNA = as.numeric(args[6])
nATAC = as.numeric(args[7])
nPET = as.numeric(args[8])
print(nRNA)
print(nATAC)
print(nPET)
# Add gene annotations
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
ucsc.levels <- str_replace(string=paste("chr",seqlevels(annotations),sep=""), pattern="chrMT", replacement="chrM")
seqlevels(annotations) <- ucsc.levels
genome(annotations) <- "mm10"
# GRanges(seqinfo(annotations))

RNAs <- Read10X(args[1])
ATACs <- Read10X_h5(args[2])
PETs = readRDS(args[4])
PETinfo = read.table(args[5],header=T)

PETinfoflt = PETinfo[PETinfo$PET_total>=nPET,]

RNA_BC <- colnames(RNAs)
ATAC_BC <- colnames(ATACs)
PET_BC <- rownames(PETinfoflt)
cor_BC <- intersect(intersect(RNA_BC,ATAC_BC),PET_BC)
print(sprintf("There are %d RNA barcodes, there are %d ATAC barcodes, there are %d PET barcodes, and there are %d shared barcodes in %s sample", length(RNA_BC), length(ATAC_BC), nrow(PETs), length(cor_BC), x))

cor_RNAs <- RNAs[,cor_BC]
cor_ATACs <- ATACs[,cor_BC]
cor_PETs <- PETs[,cor_BC]
cor_PETinfo <- PETinfo[cor_BC,]

# RNA
mydat <- CreateSeuratObject(cor_RNAs, project = x, min.cells = 5)
mydat[['percent.mt']] <- PercentageFeatureSet(object = mydat, pattern = "^mt-")
# grep("^mt-", rownames(mydat), value=T)
mydat[['percent.ribo']] <- PercentageFeatureSet(object = mydat, pattern = "^Rp[sl]")
# grep("^Rp[sl]", rownames(mydat), value=T)
mydat[['percent.hb']] <- PercentageFeatureSet(object = mydat, pattern = "^Hb[^(p)]")
# grep("^Hb[^(p)]", rownames(mydat), value=T)
mydat <- mydat[!grepl("^mt-", rownames(mydat)), ]
mydat <- mydat[!grepl("^Rp[sl]", rownames(mydat)), ]

# Add ATAC data
grange.counts <- StringToGRanges(rownames(cor_ATACs), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
cor_ATACs <- cor_ATACs[as.vector(grange.use),]

chrom_assay <- CreateChromatinAssay(counts = cor_ATACs, sep = c(":", "-"), genome = 'mm10', fragments = args[3], min.cells = 5, annotation = annotations )
mydat[["ATAC"]] <- chrom_assay
DefaultAssay(mydat) <- "ATAC"
mydat <- TSSEnrichment(mydat)
mydat <- NucleosomeSignal(mydat)
mydat$blacklist_fraction <- FractionCountsInRegion(object = mydat, assay = 'ATAC', regions = blacklist_mm10)

# Add PET data
mydat[["PET"]] <- CreateAssayObject(cor_PETs)
DefaultAssay(mydat) <- "PET"
mydat@assays$PET@counts = mydat@assays$PET@data # has no raw count, don't use count assay/nCount_PET/nFeature_PET
mydat <- AddMetaData(object = mydat, metadata = cor_PETinfo[rownames(cor_PETinfo) %in% rownames(mydat@meta.data),])


# filter
# sorted_vector <- sort(mydat$nCount_RNA)
# sorted_vector[length(sorted_vector) - 7999]
# sorted_vector <- sort(mydat$nCount_ATAC)
# sorted_vector[length(sorted_vector) - 7999]
# sorted_vector <- sort(mydat$PET_total)
# sorted_vector[length(sorted_vector) - 7999]
# summary(mydat$nCount_RNA)
# summary(mydat$nCount_ATAC)
# summary(mydat$PET_total)
mydatflt <- subset(mydat, subset = nCount_RNA > nRNA & nCount_ATAC > nATAC & PET_total > nPET & percent.mt < 20 & percent.ribo < 50 & nCount_RNA < 20000 & nCount_ATAC < 200000)
cat(sprintf("There are %d cell in %s sample\n", ncol(mydatflt), x))


# DoubletFinder
DefaultAssay(mydatflt) = "RNA"
mydatflt <- NormalizeData(object = mydatflt)
mydatflt <- FindVariableFeatures(object = mydatflt, nfeatures = 3000)
mydatflt <- ScaleData(object = mydatflt, features = rownames(mydatflt))
mydatflt <- RunPCA(mydatflt, features = VariableFeatures(mydatflt))
mydatflt <- FindNeighbors(object = mydatflt, dims = 1:30)
mydatflt <- FindClusters(object = mydatflt, resolution = 0.8)
mydatflt <- RunUMAP(mydatflt, dims = 1:30, verbose = F)

pcselect = 30
sweep = paramSweep_v3(mydatflt, PCs=1:pcselect, sct=F)
sweep.stats <- summarizeSweep(sweep, GT = FALSE)
pdf(paste0(args[10],"/",x,".pk.test.pdf"))
bcmvn <- find.pK(sweep.stats)
dev.off()
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()

# DoubletRate = 0.054 # 7000 cell: doubletsrate ~ 5.4%
DoubletRate = ncol(mydatflt)*8*1e-6
homotypic.prop <- modelHomotypic(mydatflt$RNA_snn_res.0.8) 
nExp_poi <- round(DoubletRate * ncol(mydatflt)) 
nExp_poi.adj <- round(nExp_poi * (1-homotypic.prop))
mydatflt <- doubletFinder_v3(mydatflt, PCs = 1:pcselect, pN = 0.25, pK = pK_bcmvn, nExp = nExp_poi.adj, reuse.pANN=F, sct =F)
# mydatflt <- doubletFinder_v3(mydatflt, pN = 0.25, pK = 0.09, nExp = DoubletRate, PCs = 1:30)

DF.name = colnames(mydatflt@meta.data)[grepl("DF.classification", colnames(mydatflt@meta.data))]
pANN.name = colnames(mydatflt@meta.data)[grepl("pANN", colnames(mydatflt@meta.data))]
mydatflt <- AddMetaData(mydatflt, mydatflt[[pANN.name]], col.name = "RNApANN")
mydatflt <- AddMetaData(mydatflt, mydatflt[[DF.name]], col.name = "RNADF.classifications")
mydatflt[[pANN.name]] <- NULL
mydatflt[[DF.name]] <- NULL

pdf(paste0(args[10],"/",x,".RNA.DoubletFinder.pdf"), width=10)
p4 = DimPlot(mydatflt, reduction ="umap", label = TRUE, label.size=6, group.by = c("RNADF.classifications", "RNA_snn_res.0.8")) + ggtitle(paste(x,"DFRNAumap",sep="."))
print(p4)
dev.off()
saveRDS(mydatflt, file=paste0(args[10],"/",x,".RNA.DoubletFinder.raw.rds"))


# RNA cluster
mydatflt <- subset(mydatflt, subset = RNADF.classifications == "Singlet")
mydatflt@reductions$pca = NULL
mydatflt@reductions$umap = NULL
DefaultAssay(mydatflt) <- "RNA"
mydatflt <- NormalizeData(object = mydatflt, verbose=F)
mydatflt <- FindVariableFeatures(object = mydatflt, nfeatures = 3000, verbose=F)
mydatflt <- ScaleData(object = mydatflt, features=rownames(mydatflt))
mydatflt <- RunPCA(object = mydatflt, features = VariableFeatures(mydatflt))
mydatflt <- RunPCA(object = mydatflt, features = VariableFeatures(mydatflt), reduction.name="pca.rna", reduction.key = "RNAPC_")
# mydatflt <- JackStraw(mydatflt, num.replicate = 100,dims = 30)
# mydatflt <- ScoreJackStraw(mydatflt, dims = 1:30) 
# p1 = JackStrawPlot(mydatflt, dims = 1:30) + ggtitle(x)
# p2 = ElbowPlot(mydatflt, ndims = 30) + ggtitle(x)
# print(p1 + p2)
mydatflt <- FindNeighbors(object = mydatflt, reduction = "pca.rna", dims = 1:30)
mydatflt <- FindClusters(object = mydatflt, resolution = 0.5)
mydatflt <- RunTSNE(mydatflt, reduction = "pca.rna", dims = 1:30, verbose = F, reduction.name="tsne.rna", reduction.key = "RNATSNE_", check_duplicates = FALSE)
mydatflt <- RunUMAP(mydatflt, reduction = "pca.rna", dims = 1:30, verbose = F, reduction.name="umap.rna", reduction.key = "RNAUMAP_")

# ATAC cluster
DefaultAssay(mydatflt) <- "ATAC"
mydatflt <- RunTFIDF(mydatflt, assay="ATAC")
mydatflt <- FindTopFeatures(mydatflt, assay="ATAC", min.cutoff = "q5")
mydatflt <- RunSVD(mydatflt, assay="ATAC", reduction.key = "ATACLSI_", reduction.name = "lsi.atac")
# p5 = ElbowPlot(mydatflt, reduction = "lsi") + ggtitle(x)
# print(p5)
mydatflt <- FindNeighbors(mydatflt, reduction = 'lsi.atac', dims = 2:20, k.param = 30)
mydatflt <- FindClusters(mydatflt, resolution = 0.5, algorithm = 3)
mydatflt <- RunUMAP(mydatflt, reduction = "lsi.atac", dims = 2:20, reduction.name = "umap.atac", reduction.key = "ATACUMAP_")
mydatflt <- RunTSNE(mydatflt, reduction = "lsi.atac", dims = 2:20, reduction.name = "tsne.atac", reduction.key = "ATACTSNE_", check_duplicates = FALSE)

# PET cluster
DefaultAssay(mydatflt) <- "PET"
mydatflt <- FindTopFeatures(mydatflt, assay="PET", min.cutoff = "q5")
mydatflt = RunSVD(mydatflt, assay = "PET", reduction.key = "PETLSI_", reduction.name = "lsi.pet")
mydatflt <- FindNeighbors(mydatflt, reduction = 'lsi.pet', dims = 1:20, k.param = 30)
mydatflt <- FindClusters(mydatflt, resolution = 0.5, algorithm = 3)
mydatflt <- RunUMAP(mydatflt, reduction = "lsi.pet", dims = 1:20, reduction.name = "umap.lsi.pet", reduction.key = "PETLSIUMAP_")
mydatflt <- RunTSNE(mydatflt, reduction = "lsi.pet", dims = 1:20, reduction.name = "tsne.lsi.pet", reduction.key = "PETLSITSNE_", check_duplicates = FALSE)

mydatflt <- FindVariableFeatures(object = mydatflt, nfeatures = 50000, verbose=F)
mydatflt <- ScaleData(object = mydatflt, features=rownames(mydatflt))
mydatflt <- RunPCA(object = mydatflt, features = VariableFeatures(mydatflt), reduction.name="pca.pet", reduction.key = "PETPC_")
mydatflt <- FindNeighbors(object = mydatflt, reduction = "pca.pet", dims = 1:20)
mydatflt <- FindClusters(object = mydatflt, resolution = 0.5)
mydatflt <- RunTSNE(mydatflt, reduction = "pca.pet", dims = 1:20, verbose = F, reduction.name="tsne.pca.pet", reduction.key = "PETPCATSNE_", check_duplicates = FALSE)
mydatflt <- RunUMAP(mydatflt, reduction = "pca.pet", dims = 1:20, verbose = F, reduction.name="umap.pca.pet", reduction.key = "PETPCAUMAP_")

saveRDS(mydatflt, file=paste0(args[10],"/",x,".cluster.rds"))

# data plot
is_decimal <- function(x) {
  is_num <- is.numeric(x)
  is_int <- round(x) == x
  
  if (is_num && !is_int) {
    x <- format(round(x, 3), nsmall = 3)
    return(as.numeric(x))
  } else {
    return(x)
  }
}

pdf(paste0(args[10],"/",x,".raw.count.pdf"), width=15)
p1 = VlnPlot(mydatflt, features = c('nCount_RNA','nFeature_RNA','percent.mt', 'percent.ribo', 'percent.hb','nCount_ATAC', 'nFeature_ATAC', 'TSS.enrichment', 'PET_total', 'PET_inter', 'PET_intragt1000'), pt.size=0, ncol = 4, group.by = 'orig.ident')
print(p1)
p1 = VlnPlot(mydatflt, features = c('nCount_RNA','nFeature_RNA','percent.mt', 'percent.ribo', 'percent.hb','nCount_ATAC', 'nFeature_ATAC', 'TSS.enrichment', 'PET_total', 'PET_inter', 'PET_intragt1000'), pt.size=0, ncol = 4, group.by = 'orig.ident', log = T)
print(p1)

mydat = mydatflt@meta.data
mydat$cell = rownames(mydat)
aa = c('nCount_RNA','nFeature_RNA','percent.mt', 'percent.ribo','percent.hb','nCount_ATAC', 'nFeature_ATAC', 'TSS.enrichment', 'PET_total', 'PET_inter', 'PET_intragt1000')
select = mydat[,c(aa[1], 'cell')]
mt = melt(select, id.vars = "cell")
med = is_decimal(median(mt$value))
mt$value = log10(mt$value)
bb = ggplot(mt, aes(x = variable , y = value)) +
       geom_violin() +
       stat_summary(fun = "median", geom = "text", label = med, position = position_dodge(width = 0.9), vjust = -0.5, size = 4) +
       theme_bw() +
       labs(x = x, y = "log10")

for (i in 2:length(aa)) {
  select = mydat[,c(aa[i], 'cell')]
  mt = melt(select, id.vars = "cell")
  med = is_decimal(median(mt$value))
  mt$value = log10(mt$value)
  cc = ggplot(mt, aes(x = variable , y = value)) +
    geom_violin() +
    stat_summary(fun = "median", geom = "text", label = med, position = position_dodge(width = 0.9), vjust = -0.5, size = 4) +
    theme_bw() +
    labs(x = x, y = "log10")
  bb = bb + cc
}
print(bb)

p3 = FeatureScatter(mydatflt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = 'orig.ident') 
p4 = FeatureScatter(mydatflt, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = 'orig.ident') 
p5 = FeatureScatter(mydatflt, feature1 = "nCount_RNA", feature2 = "nCount_ATAC", group.by = 'orig.ident') 
p6 = FeatureScatter(mydatflt, feature1 = "nCount_RNA", feature2 = "nFeature_ATAC", group.by = 'orig.ident') 
p7 = FeatureScatter(mydatflt, feature1 = "nCount_RNA", feature2 = "PET_total", group.by = 'orig.ident') 
p8 = FeatureScatter(mydatflt, feature1 = "nCount_RNA", feature2 = "PET_inter", group.by = 'orig.ident') 
p9 = FeatureScatter(mydatflt, feature1 = "nCount_RNA", feature2 = "PET_intragt1000", group.by = 'orig.ident') 
print(p3 + p4 + p5 + p6 + p7 + p8 + p9)
dev.off()


pdf(paste0(args[10],"/",x,".cluster.pdf"), width=15)
p1 = DimPlot(mydatflt, reduction ="tsne.rna", label = TRUE, label.size=6, group.by = "RNA_snn_res.0.5") + ggtitle(paste(x,"RNAtsne",sep="."))
p2 = DimPlot(mydatflt, reduction ="umap.rna", label = TRUE, label.size=6, group.by = "RNA_snn_res.0.5") + ggtitle(paste(x,"RNAumap",sep="."))
print(p1+p2)
p1 = DimPlot(mydatflt, reduction ="tsne.atac", label = TRUE, label.size=6, group.by = "ATAC_snn_res.0.5") + ggtitle(paste(x,"ATACtsne",sep="."))
p2 = DimPlot(mydatflt, reduction ="umap.atac", label = TRUE, label.size=6, group.by = "ATAC_snn_res.0.5") + ggtitle(paste(x,"ATACumap",sep="."))
print(p1+p2)
p1 = DimPlot(mydatflt, reduction ="tsne.lsi.pet", label = TRUE, label.size=6, group.by = "PET_snn_res.0.5") + ggtitle(paste(x,"PETlsitsne",sep="."))
p2 = DimPlot(mydatflt, reduction ="umap.lsi.pet", label = TRUE, label.size=6, group.by = "PET_snn_res.0.5") + ggtitle(paste(x,"PETlsiumap",sep="."))
print(p1+p2)
p1 = DimPlot(mydatflt, reduction ="tsne.pca.pet", label = TRUE, label.size=6, group.by = "PET_snn_res.0.5") + ggtitle(paste(x,"PETpcatsne",sep="."))
p2 = DimPlot(mydatflt, reduction ="umap.pca.pet", label = TRUE, label.size=6, group.by = "PET_snn_res.0.5") + ggtitle(paste(x,"PETpcaumap",sep="."))
print(p1+p2)
dev.off()



# WNN analysis
mydatflt <- FindMultiModalNeighbors(mydatflt, reduction.list = list("pca.rna", "lsi.atac", 'pca.pet'), dims.list = list(1:30, 2:20, 1:20))
mydatflt <- FindClusters(mydatflt, graph.name = "wsnn", algorithm = 3, resolution = 0.5)
mydatflt <- RunUMAP(mydatflt, nn.name = "weighted.nn", reduction.name = "wnn3.umap", reduction.key = "wnn3UMAP_")

pdf(paste0(args[10],"/",x,".wnn.cluster.pdf"))
p9 <- DimPlot(mydatflt, reduction = "wnn3.umap", label = TRUE, label.size = 6) + ggtitle("UMAPWNN3")
print(p9)
dev.off()

saveRDS(mydatflt, file=paste0(args[10],"/",x,".final.joint3.rds"))
