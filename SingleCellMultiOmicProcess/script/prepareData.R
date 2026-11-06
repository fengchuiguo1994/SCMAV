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
args=commandArgs(T) # 
# args=c("RNA/raw_matrix", "ATAC/raw_matrix", "DNAalign/mouseBrainP365_fragments.tsv.gz", "DNAalign/mouseBrainP365.scPET.alldata.bedpe.mat.svd.rds", "DNAalign/mouseBrainP365.scPET.alldata.bedpe.stat", 500, 3000, 500, "P365")
# args=c("RNA/raw_matrix", "ATAC/raw_matrix", "fragment/mouseBrainP365_fragments.tsv.gz", "PET/mouseBrainP365.scPET.alldata.bedpe.mat.svd.rds", "PET/mouseBrainP365.scPET.alldata.bedpe.stat", 400, 4000, 300, "P365")

set.seed(1234)
x = args[9]
print(x)
gexpath = args[1]
atacpath = args[2]
fragmentpath = args[3]
PETpath = args[4]
PETinfopath = args[5]
nRNA = as.numeric(args[6])
nATAC = as.numeric(args[7])
nPET = as.numeric(args[8])
print(nRNA)
print(nATAC)
print(nPET)

mt=20
ribo=50
maxnRNA=50000
maxnATAC=100000
DoubletRate = 0.046 # 7000 cell: doubletsrate ~ 5.4%; 6000 cell: doubletsrate ~ 4.6%
# DoubletRate = ncol(SCG)*8*1e-6 

# Add gene annotations
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
ucsc.levels <- str_replace(string=paste("chr",seqlevels(annotations),sep=""), pattern="chrMT", replacement="chrM")
seqlevels(annotations) <- ucsc.levels
genome(annotations) <- "mm10"

cutoff = function(x) {
  x = sort(x, decreasing = T)
  cut = x[1]
  index = 1
  for (i in seq(1, length(x)-1)) {
    if (x[i] > x[i+1]*5) {
      cut = x[i]
      index = i
      break
    }
  }
  print(c(x[index], x[index+1]))
}

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






RNAs <- Read10X(gexpath)
ATACs <- Read10X(atacpath)
PETs = readRDS(PETpath)
PETinfo = read.table(PETinfopath, header=T)

RNA_BC <- colnames(RNAs)
ATAC_BC <- colnames(ATACs)
PET_BC <- colnames(PETs$rawMat)
cor_BC <- intersect(intersect(RNA_BC,ATAC_BC),PET_BC)
print(sprintf("There are %d RNA barcodes, there are %d ATAC barcodes, there are %d PET barcodes, and there are %d shared barcodes in %s sample", length(RNA_BC), length(ATAC_BC), nrow(PETs$rawMat), length(cor_BC), x))
cor_RNAs <- RNAs[,cor_BC]
cor_ATACs <- ATACs[,cor_BC]
rawMat <- PETs$rawMat[,cor_BC]
SVDMat <- PETs$SVDMat[,cor_BC]
cor_PETinfo <- PETinfo[cor_BC,]

# RNA
mydat <- CreateSeuratObject(cor_RNAs, project = x, min.cells = 0, assay = "RNA")
mydat[['percent.mt']] <- PercentageFeatureSet(object = mydat, pattern = "^mt-")
mydat[['percent.ribo']] <- PercentageFeatureSet(object = mydat, pattern = "^Rp[sl]")
mydat[['percent.hb']] <- PercentageFeatureSet(object = mydat, pattern = "^Hb[^(p)]")
# grep("^mt-", rownames(mydat), value=T)
# grep("^Rp[sl]", rownames(mydat), value=T)
# grep("^Hb[^(p)]", rownames(mydat), value=T)
# mydat <- mydat[!grepl("Gm42418", rownames(mydat)), ]
# mydat <- mydat[!grepl("^mt-", rownames(mydat)), ]
# mydat <- mydat[!grepl("^Rp[sl]", rownames(mydat)), ]

# Add ATAC data
grange.counts <- StringToGRanges(rownames(cor_ATACs), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
cor_ATACs <- cor_ATACs[as.vector(grange.use),]
chrom_assay <- CreateChromatinAssay(counts = cor_ATACs, sep = c(":", "-"), genome = 'mm10', fragments = fragmentpath, min.cells = 0, annotation = annotations )
mydat[["ATAC"]] <- chrom_assay

# Add PET data
mydat[["PET"]] <- CreateAssayObject(counts = rawMat)
# table(rownames(mydat@assays$PET@counts) == gsub("_", "-", rownames(SVDMat)))
# table(colnames(mydat@assays$PET@counts) == colnames(SVDMat))
rownames(SVDMat) = rownames(mydat@assays$PET@counts)
mydat@assays$PET@data = SVDMat
mydat <- AddMetaData(object = mydat, metadata = cor_PETinfo[rownames(cor_PETinfo) %in% rownames(mydat@meta.data),])

# filter low quality cell  
# filter # cell_number = 1.2 * length(unique(mydat$orig.ident)) * 6000
# sorted_vector <- sort(mydat$nCount_RNA)
# sorted_vector[length(sorted_vector) - cell_number][1]
# sorted_vector <- sort(mydat$nCount_ATAC)
# sorted_vector[length(sorted_vector) - cell_number][1]
# sorted_vector <- sort(mydat$PET_total)
# sorted_vector[length(sorted_vector) - cell_number][1]
# summary(mydat$nCount_RNA)
# summary(mydat$nCount_ATAC)
# summary(mydat$PET_total)
# cutoff(mydat$nCount_RNA)
# cutoff(mydat$nCount_ATAC)
# cutoff(mydat$PET_total)


SCG <- subset(mydat, subset = nCount_RNA >= nRNA & nCount_ATAC >= nATAC & PET_total >= nPET & percent.mt < mt & percent.ribo < ribo & nCount_RNA < maxnRNA & nCount_ATAC < maxnATAC)
cat(sprintf("There are %d cell in %s sample\n", ncol(SCG), x))
DefaultAssay(SCG) <- "ATAC"
SCG <- TSSEnrichment(SCG)
SCG <- NucleosomeSignal(SCG)
SCG$blacklist_fraction <- FractionCountsInRegion(object = SCG, assay = 'ATAC', regions = blacklist_mm10)

DefaultAssay(mydat) <- "RNA"
# Identify abnormally high expressed genes, can find the gene Gm42418, need to be removed
C <- SCG@assays$RNA@counts
C <- Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100
most_expressed <- order(apply(C, 1, median), decreasing = T)[20:1]
pdf(paste0(x, ".raw.topexpgene.pdf"), width=15)
boxplot(as.matrix(t(as.data.frame(C[most_expressed, ]))), cex = 0.1, las = 1, xlab = "% total count per cell", col = (scales::hue_pal())(20)[20:1], horizontal = TRUE, main=x)
dev.off()

# RNA cluster
DefaultAssay(SCG) = "RNA"
SCG <- NormalizeData(object = SCG, verbose=F)
SCG <- FindVariableFeatures(object = SCG, nfeatures = 3000, verbose=F)
SCG <- ScaleData(object = SCG, features=rownames(SCG))
SCG <- RunPCA(SCG, features = VariableFeatures(SCG))
SCG <- FindNeighbors(object = SCG, dims = 1:30)
SCG <- FindClusters(object = SCG, resolution = 1)
SCG <- RunUMAP(SCG, dims = 1:30, verbose = F)

# DoubletFinder
pcselect = 30
sweep = paramSweep_v3(SCG, PCs=1:pcselect, sct = F)
sweep.stats <- summarizeSweep(sweep, GT = FALSE)
pdf(paste0(x,".pk.test.pdf"))
bcmvn <- find.pK(sweep.stats)
dev.off()
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()

homotypic.prop <- modelHomotypic(SCG$RNA_snn_res.1) 
nExp_poi <- round(DoubletRate * ncol(SCG)) 
nExp_poi.adj <- round(nExp_poi * (1-homotypic.prop))
SCG <- doubletFinder_v3(SCG, PCs = 1:pcselect, pN = 0.25, pK = pK_bcmvn, nExp = nExp_poi.adj, reuse.pANN=F, sct = F)
# SCG <- doubletFinder_v3(SCG, pN = 0.25, pK = 0.09, nExp = DoubletRate, PCs = 1:30)

DF.name = colnames(SCG@meta.data)[grepl("DF.classification", colnames(SCG@meta.data))]
pANN.name = colnames(SCG@meta.data)[grepl("pANN", colnames(SCG@meta.data))]
SCG <- AddMetaData(SCG, SCG[[pANN.name]], col.name = "RNApANN")
SCG <- AddMetaData(SCG, SCG[[DF.name]], col.name = "RNADF")
SCG[[pANN.name]] <- NULL
SCG[[DF.name]] <- NULL

pdf(paste0(x,".RNA.DoubletFinder.pdf"), width=10)
p4 = DimPlot(SCG, reduction ="umap", label = TRUE, label.size=6, group.by = c("RNADF", "RNA_snn_res.1")) + ggtitle(paste(x,"DFRNAumap",sep="."))
print(p4)
dev.off()
saveRDS(SCG, file=paste0(x,".DoubletFinder.raw.rds"))
save(SCG, file = paste0(x,".DoubletFinder.raw.rdata"))
SCG <- subset(SCG, subset = RNADF == "Singlet")

# plot stat
pdf(paste0(x,".count.pdf"), width=15)
p1 = VlnPlot(SCG, features = c('nCount_RNA','nFeature_RNA','percent.mt', 'percent.ribo', 'percent.hb','nCount_ATAC', 'nFeature_ATAC', 'TSS.enrichment', 'PET_total', 'PET_inter', 'PET_intragt1000'), pt.size=0, ncol = 4, group.by = 'orig.ident')
print(p1)
p1 = VlnPlot(SCG, features = c('nCount_RNA','nFeature_RNA','percent.mt', 'percent.ribo', 'percent.hb','nCount_ATAC', 'nFeature_ATAC', 'TSS.enrichment', 'PET_total', 'PET_inter', 'PET_intragt1000'), pt.size=0, ncol = 4, group.by = 'orig.ident', log = T)
print(p1)

mydat = SCG@meta.data
mydat$cell = rownames(mydat)
aa = c('nCount_RNA','nFeature_RNA','percent.mt', 'percent.ribo','percent.hb','nCount_ATAC', 'nFeature_ATAC', 'TSS.enrichment', 'PET_total', 'PET_inter', 'PET_intragt1000')
select = mydat[,c(aa[1], 'cell')]
mt = melt(select, id.vars = "cell")
med = is_decimal(median(mt$value))
# mt$value = log10(mt$value)
bb = ggplot(mt, aes(x = variable , y = value)) +
      geom_violin() +
      geom_boxplot(width = 0.1, outlier.size = 1, show.legend = FALSE, fatten = 2, position = position_nudge(x = 0.0), outlier.shape = NA) +
      stat_summary(fun = "median", geom = "text", label = med, position = position_dodge(width = 0.9), vjust = -0.5, size = 4) +
      scale_y_log10() +
      theme_bw() +
      theme(panel.grid=element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())

for (i in 2:length(aa)) {
  select = mydat[,c(aa[i], 'cell')]
  mt = melt(select, id.vars = "cell")
  med = is_decimal(median(mt$value))
  # mt$value = log10(mt$value)
  cc = ggplot(mt, aes(x = variable , y = value)) +
    geom_violin() +
    geom_boxplot(width = 0.1, outlier.size = 1, show.legend = FALSE, fatten = 2, position = position_nudge(x = 0.0), outlier.shape = NA) +
    stat_summary(fun = "median", geom = "text", label = med, position = position_dodge(width = 0.9), vjust = -0.5, size = 4) +
    scale_y_log10() +
    # stat_summary(fun = "median", geom = "crossbar", colour = "black") +
    theme_bw() +
    theme(panel.grid=element_blank(), , axis.title.x = element_blank(), axis.title.y = element_blank())
  bb = bb + cc
}
print(bb)

p3 = FeatureScatter(SCG, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = 'orig.ident') 
p4 = FeatureScatter(SCG, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = 'orig.ident') 
p5 = FeatureScatter(SCG, feature1 = "nCount_RNA", feature2 = "nCount_ATAC", group.by = 'orig.ident') 
p6 = FeatureScatter(SCG, feature1 = "nCount_RNA", feature2 = "nFeature_ATAC", group.by = 'orig.ident') 
p7 = FeatureScatter(SCG, feature1 = "nCount_RNA", feature2 = "PET_total", group.by = 'orig.ident') 
p8 = FeatureScatter(SCG, feature1 = "nCount_RNA", feature2 = "PET_inter", group.by = 'orig.ident') 
p9 = FeatureScatter(SCG, feature1 = "nCount_RNA", feature2 = "PET_intragt1000", group.by = 'orig.ident') 
print(p3 + p4 + p5 + p6 + p7 + p8 + p9)
dev.off()

saveRDS(SCG, file = paste0(x,".DoubletFinder.raw.rds"))
write.table(SCG@meta.data, file = paste0(x,".DoubletFinder.raw.metacell.txt"), quote=F, sep="\t")
DefaultAssay(SCG) = "RNA"
SCG@reductions$pca = NULL
SCG@reductions$umap = NULL
saveRDS(SCG, file = paste0(x,".prepareData.rds"))

# SCG@assays$RNA@counts[1:20,1:3]
# SCG@assays$RNA@data[1:20,1:3]
# SCG@assays$RNA@scale.data[1:20,1:3]

# SCG@assays$ATAC@counts[1:20,1:3]
# SCG@assays$ATAC@data[1:20,1:3]
# SCG@assays$ATAC@scale.data[1:20,1:3]

# SCG@assays$PET@counts[1:20,1:3]
# SCG@assays$PET@data[1:20,1:3]
# SCG@assays$PET@scale.data[1:20,1:3]
