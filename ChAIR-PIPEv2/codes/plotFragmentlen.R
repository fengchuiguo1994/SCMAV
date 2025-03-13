library(ggplot2)
library(patchwork)
args=commandArgs(T)
mydat = read.table(args[1])

nf = sum(mydat$V2)
cutoff = mean(mydat[mydat$V1 %in% 140:150 ,]$V2)
rm = mydat[mydat$V1 %in% 140:150 & mydat$V2 > 2*cutoff,]$V1
mydat = mydat[!mydat$V1 %in% rm,]
xmax = max(mydat$V1)
nf2 = sum(mydat$V2)

p1 = ggplot(data = mydat, mapping = aes(V1,V2)) +
  geom_point(color="#f8b62b") +
  labs(title=paste0("Fragment length distribution (n=",nf,") (filter n=",nf2,")")) +
  ylab('Number') +
  xlab("Fragment length") +
  scale_x_continuous(breaks = seq(0, xmax, by = 200)) +
  theme_bw() +
  theme(panel.grid=element_blank())

pdf(args[2])
print(p1)
dev.off()
