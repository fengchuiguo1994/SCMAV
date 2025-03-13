library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(Matrix)
library(plotly)
library(plot3D)
args <- commandArgs(trailingOnly = TRUE)

RNA = read.table(args[1], header=T)
ATAC = read.table(args[2], header=T)
PET = read.table(args[3], header=T)
PET = PET[, c("nPET",	"barcode")]

pdf(args[4])

bb = as.data.frame(table(RNA$UMI))
bb$Var1 = as.numeric(as.character(bb$Var1))
# full_numbers <- data.frame(Var1 = 1:max(bb$Var1))
# cc <- full_numbers %>% left_join(bb, by = "Var1") %>% replace_na(list(Freq = 0))
# cc = cc[order(-cc$Var1), ]
cc = bb[order(-bb$Var1), ]
Cumulative <- cc %>% arrange(desc(Var1)) %>% mutate(cumulative_count = cumsum(Freq))
Cumulative$Var1_log = log10(Cumulative$Var1)
Cumulative$cumulative_count_log = log10(Cumulative$cumulative_count)
p1 = ggplot(data = Cumulative, mapping = aes(cumulative_count_log, Var1_log)) +
  # geom_line(color="black") +
  geom_step(linetype = "solid", color = "black") + 
  labs(title="scRNA Barcode Rank Plot") +
  ylab('UMI count') +
  xlab("Barcodes") +
  scale_x_continuous(breaks = seq(0, 8, by = 1), labels = c("1", "10", "100", "1k", "10k", "100k","1M","10M","100M")) +
  scale_y_continuous(breaks = seq(0, 8, by = 1), labels = c("1", "10", "100", "1k", "10k", "100k","1M","10M","100M")) +
  theme_bw() +
  theme(panel.grid=element_blank())
print(p1)

p <- plot_ly(Cumulative, x = ~cumulative_count_log, y = ~Var1_log, mode = 'lines', type = 'scatter', text = ~paste("barcode: ", Var1, "UMI: ", cumulative_count), hoverinfo = "text", line = list(shape = 'hv')) %>%
  layout(title = "scRNA Barcode Rank Plot",
         xaxis = list(title = "Barcodes", automargin = TRUE, tickvals = log10(c(1,10,100,1000,10000,100000,1000000,10000000,100000000)), ticktext=c("1", "10", "100", "1k", "10k", "100k","1M","10M","100M")),
         yaxis = list(title = "UMI count", automargin = TRUE, tickvals = log10(c(1,10,100,1000,10000,100000,1000000,10000000,100000000)),ticktext=c("1", "10", "100", "1k", "10k", "100k","1M","10M","100M"))
  )
htmlwidgets::saveWidget(p, paste0(args[4], ".scRNA.html"))


bb = as.data.frame(table(ATAC$nFragment))
bb$Var1 = as.numeric(as.character(bb$Var1))
# full_numbers <- data.frame(Var1 = 1:max(bb$Var1))
# cc <- full_numbers %>% left_join(bb, by = "Var1") %>% replace_na(list(Freq = 0))
# cc = cc[order(-cc$Var1), ]
cc = bb[order(-bb$Var1), ]
Cumulative <- cc %>% arrange(desc(Var1)) %>% mutate(cumulative_count = cumsum(Freq))
Cumulative$Var1_log = log10(Cumulative$Var1)
Cumulative$cumulative_count_log = log10(Cumulative$cumulative_count)
p1 = ggplot(data = Cumulative, mapping = aes(cumulative_count_log, Var1_log)) +
  # geom_line(color="black") +
  geom_step(linetype = "solid", color = "black") + 
  labs(title="scATAC Barcode Rank Plot") +
  ylab('nFragment count') +
  xlab("Barcodes") +
  scale_x_continuous(breaks = seq(0, 8, by = 1), labels = c("1", "10", "100", "1k", "10k", "100k","1M","10M","100M")) +
  scale_y_continuous(breaks = seq(0, 8, by = 1), labels = c("1", "10", "100", "1k", "10k", "100k","1M","10M","100M")) +
  theme_bw() +
  theme(panel.grid=element_blank())
print(p1)

p <- plot_ly(Cumulative, x = ~cumulative_count_log, y = ~Var1_log, mode = 'lines', type = 'scatter', text = ~paste("barcode: ", Var1, "nFragment: ", cumulative_count), hoverinfo = "text", line = list(shape = 'hv')) %>%
  layout(title = "scATAC Barcode Rank Plot",
         xaxis = list(title = "Barcodes", automargin = TRUE, tickvals = log10(c(1,10,100,1000,10000,100000,1000000,10000000,100000000)), ticktext=c("1", "10", "100", "1k", "10k", "100k","1M","10M","100M")),
         yaxis = list(title = "nFragment count", automargin = TRUE, tickvals = log10(c(1,10,100,1000,10000,100000,1000000,10000000,100000000)),ticktext=c("1", "10", "100", "1k", "10k", "100k","1M","10M","100M"))
  )
htmlwidgets::saveWidget(p, paste0(args[4], ".scATAC.html"))


bb = as.data.frame(table(PET$nPET))
bb$Var1 = as.numeric(as.character(bb$Var1))
# full_numbers <- data.frame(Var1 = 1:max(bb$Var1))
# cc <- full_numbers %>% left_join(bb, by = "Var1") %>% replace_na(list(Freq = 0))
# cc = cc[order(-cc$Var1), ]
cc = bb[order(-bb$Var1), ]
Cumulative <- cc %>% arrange(desc(Var1)) %>% mutate(cumulative_count = cumsum(Freq))
Cumulative$Var1_log = log10(Cumulative$Var1)
Cumulative$cumulative_count_log = log10(Cumulative$cumulative_count)
p1 = ggplot(data = Cumulative, mapping = aes(cumulative_count_log, Var1_log)) +
  # geom_line(color="black") +
  geom_step(linetype = "solid", color = "black") + 
  labs(title="scPET Barcode Rank Plot") +
  ylab('nPET count') +
  xlab("Barcodes") +
  scale_x_continuous(breaks = seq(0, 8, by = 1), labels = c("1", "10", "100", "1k", "10k", "100k","1M","10M","100M")) +
  scale_y_continuous(breaks = seq(0, 8, by = 1), labels = c("1", "10", "100", "1k", "10k", "100k","1M","10M","100M")) +
  theme_bw() +
  theme(panel.grid=element_blank())
print(p1)

p <- plot_ly(Cumulative, x = ~cumulative_count_log, y = ~Var1_log, mode = 'lines', type = 'scatter', text = ~paste("barcode: ", Var1, "nPET: ", cumulative_count), hoverinfo = "text", line = list(shape = 'hv')) %>%
  layout(title = "scPET Barcode Rank Plot",
         xaxis = list(title = "Barcodes", automargin = TRUE, tickvals = log10(c(1,10,100,1000,10000,100000,1000000,10000000,100000000)), ticktext=c("1", "10", "100", "1k", "10k", "100k","1M","10M","100M")),
         yaxis = list(title = "nPET count", automargin = TRUE, tickvals = log10(c(1,10,100,1000,10000,100000,1000000,10000000,100000000)),ticktext=c("1", "10", "100", "1k", "10k", "100k","1M","10M","100M"))
  )
htmlwidgets::saveWidget(p, paste0(args[4], ".scPET.html"))



outer = merge(merge(RNA,ATAC,all = TRUE),PET, all = TRUE)
outer[is.na(outer)] = 0
inner = merge(merge(RNA,ATAC),PET)

with(outer, scatter3D(UMI, nFragment, nPET, col="#E64B35FF", theta = 60, phi = 20,
                        xlab = "UMI", ylab = "nFragment", zlab = "nPET",main = paste0("There are: ", nrow(outer), " barcodes"), cex = 1,
                        bty = "f",ticktype = "detailed",d = 3, clab = c("type","class"),
                        adj = 0.5,font = 2,type = "s",pch = 16,box = TRUE,bg="black",alpha=0.4))

outercg = outer[outer$UMI > 200 & outer$nFragment > 1000 & outer$nPET > 50, ]
if (nrow(outercg) > 0) {
  if (nrow(outercg[outercg$UMI > 5000,]) > 0) {
    outercg[outercg$UMI > 5000,]$UMI = 5000
  }
  if (nrow(outercg[outercg$nFragment > 50000,]) > 0) {
    outercg[outercg$nFragment > 50000,]$nFragment = 50000
  }
  if (nrow(outercg[outercg$nPET > 500,]) > 0) {
    outercg[outercg$nPET > 500,]$UMI = 500
  }
  with(outercg, scatter3D(UMI, nFragment, nPET, col="#E64B35FF", theta = 60, phi = 20,
                    xlab = "UMI", ylab = "nFragment", zlab = "nPET",main = paste0("There are: ", nrow(outercg), " barcodes with UMI>200, nFragment>1000, nPET>50"), cex = 1,
                    bty = "f",ticktype = "detailed",d = 3, clab = c("type","class"),
                    adj = 0.5,font = 2,type = "s",pch = 16,box = TRUE,bg="black",alpha=0.4))
} else {
  p <- ggplot() +
  annotate("text", x = 0, y = 0, label = "There are 0 barcode with UMI>200, nFragment>1000, nPET>50", size = 5, color = "blue") +
  xlim(-1, 1) + ylim(-1, 1) +
  theme_void()
print(p)
}
dev.off()
