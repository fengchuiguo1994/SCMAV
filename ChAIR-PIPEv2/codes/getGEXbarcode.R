library(dplyr)
library(tidyr)
library(Matrix)

args <- commandArgs(trailingOnly = TRUE);

dir <- args[1]
gf = gzfile(paste0(dir,"/barcodes.tsv.gz"), 'rt')
barcodes = read.table(gf)
close(gf)
names(barcodes) = c("barcode")
barcodes$V2 = 1:nrow(barcodes)
  
gf = gzfile(paste0(dir,"/features.tsv.gz"), 'rt')
features = read.table(gf)
close(gf)
features = features[c(1)]
names(features) = c("gene")
features$V1 = 1:nrow(features)

gf = gzfile(paste0(dir,"/matrix.mtx.gz"), 'rt')
mydat = read.table(gf, skip=3)
close(gf)
mat = left_join(left_join(mydat,features,by="V1"), barcodes, by = "V2")

sparse_matrix <- sparseMatrix(
  i = match(mat$gene, unique(mat$gene)),
  j = match(mat$barcode, unique(mat$barcode)),
  x = mat$V3,
  dimnames = list(unique(mat$gene), unique(mat$barcode))
)

RNA = as.data.frame(colSums(sparse_matrix))
names(RNA) = ("UMI")
RNA$barcode = rownames(RNA)

write.table(RNA, file = args[2], sep="\t", quote=F)