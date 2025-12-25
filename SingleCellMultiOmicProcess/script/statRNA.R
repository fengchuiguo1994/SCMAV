args <- commandArgs(trailingOnly = TRUE)

if (args[2] == "10X"){
    mydat = read.csv(args[1], stringsAsFactors=F)
    cat('RNA_raw_reads\t', mydat[1,]$Number.of.Reads, "\n", sep="")
    cat('RNA_barcode_reads_ratio\t', mydat[1,]$Valid.Barcodes, "\n", sep="")
    cat('RNA_mapped_ratio\t', mydat[1,]$Reads.Mapped.to.Genome, "\n", sep="")
    cat('RNA_detect_gene\t', mydat[1,]$Total.Genes.Detected, "\n", sep="")
    cat('RNA_Estimated_cells\t', mydat[1,]$Estimated.Number.of.Cells, "\n", sep="")
    cat('RNA_UMI\t', mydat[1,]$Median.UMI.Counts.per.Cell, "\n", sep="")
    cat('RNA_nGeneFeature\t', mydat[1,]$Median.Genes.per.Cell, "\n", sep="")
} else if (args[2] == "STARsolo") {
    mydat = read.csv(args[1], header=F, stringsAsFactors=F)
    cat('RNA_raw_reads\t', prettyNum(mydat[mydat$V1=="Number of Reads","V2"], big.mark = ",", scientific = FALSE), "\n", sep="")
    cat('RNA_barcode_reads_ratio\t', sprintf("%.2f%%", mydat[mydat$V1=="Reads With Valid Barcodes","V2"] * 100), "\n", sep="")
    cat('RNA_mapped_ratio\t', sprintf("%.2f%%", mydat[mydat$V1=="Reads Mapped to Genome: Unique+Multiple","V2"] * 100), "\n", sep="")
    cat('RNA_detect_gene\t', prettyNum(mydat[mydat$V1=="Total Gene Detected","V2"], big.mark = ",", scientific = FALSE), "\n", sep="")
    cat('RNA_Estimated_cells\t', prettyNum(mydat[mydat$V1=="Estimated Number of Cells","V2"], big.mark = ",", scientific = FALSE), "\n", sep="")
    cat('RNA_UMI\t', prettyNum(mydat[mydat$V1=="Median UMI per Cell","V2"], big.mark = ",", scientific = FALSE), "\n", sep="")
    cat('RNA_nGeneFeature\t', prettyNum(mydat[mydat$V1=="Median Gene per Cell","V2"], big.mark = ",", scientific = FALSE), "\n", sep="")
}