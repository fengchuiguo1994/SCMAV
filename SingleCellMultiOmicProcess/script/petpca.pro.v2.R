# library(BandNorm)
# library(Signac)
library(dplyr)
library(data.table)
library(Matrix)
args=commandArgs(T)
set.seed(1234)

# ref: https://github.com/sshen82/BandNorm/blob/main/R/bandnorm.R
# bandnorm is not installed well in my docker image

bandnormxyh = function(path = NULL, hic_df = NULL, save = TRUE, save_path = NULL) {
  if (is.null(path) && is.null(hic_df)){
    stop("Please specify one of path or hic_df.")
  }
  if (!is.null(path) && !is.null(hic_df)){
    warning("Specified both path and hic_df. Will use path file as default for saving memory.")
  }
  if (save) {
    if (is.null(save_path)){
      stop("path for saving the data is NULL!")
    }
    if (!dir.exists(save_path)){
      warning("path for saving the data doesn't exist, will create one according to the directory.")
      dir.create(save_path, recursive = TRUE)
    }
  }
  # Get path and name for all the cells in this path.
  if (!is.null(path)){
    if (!dir.exists(path)){
      stop("path for data doesn't exist!")
    }
    paths = list.files(path, recursive = TRUE, full.names = TRUE)
    names = basename(list.files(path, recursive = TRUE))
    # The input format of the cell should be [chr1, bin1, chr2, bin2, count].
    load_cell = function(i) {
      return(fread(paths[i], select = c(1, 2, 4, 5)) %>% dplyr::rename(chrom = V1, binA = V2, binB = V4, count = V5) %>%
               mutate(diag = abs(binB - binA), cell = names[i]))
    }
    hic_df = rbindlist(lapply(1:length(paths), load_cell))
  }
  # Calculate band depth and the mean of band depth for bandnorm.
  band_info <- hic_df %>% group_by(chrom, diag, cell) %>% summarise(band_depth = sum(count))
  alpha_j <- band_info %>% group_by(chrom, diag) %>% summarise(depth = mean(band_depth))

  hic_df <- hic_df %>% left_join(alpha_j, by = c("chrom", "diag")) %>% left_join(band_info, by = c("chrom", "diag", "cell")) %>% mutate(BandNorm = count/band_depth * depth) %>% select(-c(band_depth, depth))
  if (save) {
    save_path = file.path(save_path, unique(names))
    for (i in 1:length(unique(names))) {
      write.table(hic_df[cell == unique(names)[i], c("chrom", "binA", "chrom", "binB", "BandNorm")], file = save_path[i], col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
    }
  }
  return(hic_df)
}

create_embeddingxyh = function (path = NULL, hic_df = NULL, chrs = paste0("chr", c(1:19, "X")), dim_pca = 50, do_harmony = FALSE, batch = NULL) 
{
    # hic_df = hic_df[hic_df$diag != 0, ] # hxy
    setDT(hic_df)
    # names = unique(hic_df$cell)
    # n = length(names)
    hash.cell <- data.frame(cell = unique(hic_df$cell), cellnum = seq(length(unique(hic_df$cell))))
    hic_df = left_join(hic_df, hash.cell, by="cell")
    hic_df = hic_df %>% mutate(featureIndex = paste(chrom, binA, binB, sep = "_")) %>% select(-c(binA, binB, diag))
    hash.featureIndex <- data.frame(featureIndex = unique(hic_df$featureIndex), featureIndexnum = seq(length(unique(hic_df$featureIndex))))
    hic_df = left_join(hic_df, hash.featureIndex, by="featureIndex")
    # hic_df = sparseMatrix(i = hic_df$cellnum, j = hic_df$featureIndexnum, x = hic_df$BandNorm, dims = c(max(hic_df$cellnum), max(hic_dfCHR$featureIndexnum)), index1 = TRUE)
    rawInput = sparseMatrix(i = hic_df$cellnum, j = hic_df$featureIndexnum, x = hic_df$count)
    SVDInput = sparseMatrix(i = hic_df$cellnum, j = hic_df$featureIndexnum, x = hic_df$BandNorm)
    rownames(rawInput) = hash.cell$cell
    colnames(rawInput) = hash.featureIndex$featureIndex
    rownames(SVDInput) = hash.cell$cell
    colnames(SVDInput) = hash.featureIndex$featureIndex
    matInput = list()
    matInput$rawMat = t(rawInput)
    matInput$SVDMat = t(SVDInput)
    return(matInput)
}


PETdata = read.table(args[1], header = T, sep = "\t", check.names = F, stringsAsFactors = F) # DNAalign/mouseBrainP365.scPET.alldata.bedpe.mat
PETdata$resolution = as.integer(PETdata$resolution)
PETdata$pos1 = as.integer(PETdata$pos1)
PETdata$pos2 = as.integer(PETdata$pos2)
PETdata$count = as.integer(PETdata$count)
PETdata = PETdata[PETdata$chrom1 == PETdata$chrom2,] # https://github.com/sshen82/BandNorm/issues/4
# table(PETdata$pos1 > PETdata$pos2)

PETdata$diag = PETdata$pos2 - PETdata$pos1
batch = unique(PETdata[, c("cell_id", "library")]) %>% as.matrix
PETdata = PETdata[, c("chrom1", "pos1", "pos2", "count", "diag", "cell_id")]
colnames(PETdata) = c("chrom", "binA", "binB", "count", "diag", "cell")

# bandnorm_res = bandnorm(hic_df = PETdata, save = F)
bandnorm_res = bandnormxyh(hic_df = PETdata, save = F)
output = create_embeddingxyh(hic_df = bandnorm_res, do_harmony = F, batch = batch, chrs = sort(unique(PETdata$chrom)))
saveRDS(output, file = args[2])






# library(BandNorm)
# PETdata = read.table(args[1], header = T, sep = "\t", check.names = F, stringsAsFactors = F) # DNAalign/mouseBrainP365.scPET.alldata.bedpe.mat
# PETdata$resolution = as.integer(PETdata$resolution)
# PETdata$pos1 = as.integer(PETdata$pos1)
# PETdata$pos2 = as.integer(PETdata$pos2)
# PETdata$count = as.integer(PETdata$count)
# PETdata = PETdata[PETdata$chrom1 == PETdata$chrom2,] # https://github.com/sshen82/BandNorm/issues/4
# # table(PETdata$pos1 > PETdata$pos2)

# PETdata$diag = PETdata$pos2 - PETdata$pos1
# batch = unique(PETdata[, c("cell_id", "library")]) %>% as.matrix
# names(batch) = c("cell_name", "batch")
# PETdata = PETdata[, c("chrom1", "pos1", "pos2", "count", "diag", "cell_id")]
# colnames(PETdata) = c("chrom", "binA", "binB", "count", "diag", "cell")

# bandnorm_result = bandnorm(hic_df = PETdata, save = FALSE)
# write.table(bandnorm_result, file="raw.1.txt", quote=F, sep="\t")
# embedding = create_embedding(hic_df = bandnorm_result, do_harmony = F, batch = batch, chrs = sort(unique(PETdata$chrom)))