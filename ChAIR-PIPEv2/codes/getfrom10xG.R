args <- commandArgs(trailingOnly = TRUE)

mydat = read.csv(args[1])

message('RNA_raw_reads\t', mydat[1,]$Number.of.Reads)
message('RNA_barcode\t', mydat[1,]$Valid.Barcodes)
message('RNA_mapped_ratio\t', mydat[1,]$Reads.Mapped.to.Genome)
message('RNA_overlap_gene\t', mydat[1,]$Total.Genes.Detected)
message('RNA_Estimated_cells\t', mydat[1,]$Estimated.Number.of.Cells)
message('RNA_UMI\t', mydat[1,]$Median.UMI.Counts.per.Cell)
message('RNA_nGeneFeature\t', mydat[1,]$Median.Genes.per.Cell)
