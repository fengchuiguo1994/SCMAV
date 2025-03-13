library(BandNorm)
# library(Signac)
args=commandArgs(T)
set.seed(1234)

create_embeddingxyh = function (path = NULL, hic_df = NULL, chrs = paste0("chr", c(1:19, "X")), dim_pca = 50, do_harmony = FALSE, batch = NULL) 
{
    n = length(unique(hic_df$cell))
    setDT(hic_df)
    # hic_df = hic_df[hic_df$diag != 0, ] # hxy
    names = unique(hic_df$cell)
    # unique(str_match(names,"^(SCG\\w+)_")[,2])
    n = length(names)
    hic_df$cell = as.numeric(factor(hic_df$cell, level = names))
    hic_df = hic_df %>% mutate(featureIndex = paste(binA, 
        binB, sep = "_")) %>% select(-c(binA, binB, diag))
    output = list()
    # return(hic_df)
    for (c in 1:length(chrs)) {
        hic_dfCHR = hic_df[chrom == chrs[c], ]
        hic_dfCHR$featureIndex = as.numeric(factor(hic_dfCHR$featureIndex))
        hic_dfCHR = sparseMatrix(i = hic_dfCHR$cell, j = hic_dfCHR$featureIndex, 
            x = hic_dfCHR$BandNorm, dims = c(max(hic_df$cell), 
                max(hic_dfCHR$featureIndex)), index1 = TRUE)
        output[[chrs[c]]] = hic_dfCHR
    }

    # return(output)
    
    SVDInput = c()
    for (i in 1:length(output)) {
        SVDInput = cbind(SVDInput, output[[i]])
    }
    rownames(SVDInput) = names
    return(t(SVDInput))
}


PETdata = read.table(args[1], header = T, sep = "\t", check.names = F, stringsAsFactors = F)
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

bandnorm_res = bandnorm(hic_df = PETdata, save = F)
SVDInput = create_embeddingxyh(hic_df = bandnorm_res, do_harmony = F, batch = batch, chrs = sort(unique(PETdata$chrom)))
rownames(SVDInput) = paste0("cbin",seq(1,nrow(SVDInput)))
# SVD = RunSVD(SVDInput,assay="ATAC",reduction.key="PETLSI_",reduction.name = "petlsi")

saveRDS(SVDInput, file = args[2])
