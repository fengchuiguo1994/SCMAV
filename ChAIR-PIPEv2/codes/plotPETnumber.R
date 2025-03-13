library(ggplot2)
library(patchwork)
args=commandArgs(T)
mydat=read.table(args[1],header=F,sep="\t",stringsAsFactors=F)

result = mydat[grep("pets_number",mydat$V1),]
result$V2 = as.numeric(gsub(",", "", result$V2))
intra = result[1:10,]
inter = result[11:20,]

names(intra) = c("PETcount","number")
intra$PETcount = factor(intra$PETcount,levels=unique(intra$PETcount))
intra$number = as.numeric(intra$number)
names(inter) = c("PETcount","number")
inter$PETcount = factor(inter$PETcount,levels=unique(inter$PETcount))
inter$number = as.numeric(inter$number)

p1 = ggplot(data=intra,mapping=aes(PETcount,number)) +
  geom_bar(stat="identity",width = 0.5,position='dodge',colour="black",size=0.25,fill="#f8b62b",alpha=1) +
  geom_text(stat='identity',mapping=aes(label=number),color="black",size=4,position=position_dodge(0.5),vjust=-0.1) +
  ylim(c(0,max(intra$number*1.02))) +
  labs(title="Intra-interaction") +
  ylab('Number') +
  xlab("PET count per interaction") +
  theme_bw() +
  theme(panel.grid=element_blank(),axis.text.x = element_text(angle = 45, hjust = 1))
p2 = ggplot(data=inter,mapping=aes(PETcount,number)) +
  geom_bar(stat="identity",width = 0.5,position='dodge',colour="black",size=0.25,fill="#f8b62b",alpha=1) +
  geom_text(stat='identity',mapping=aes(label=number),color="black",size=4,position=position_dodge(0.5),vjust=-0.1) +
  ylim(c(0,max(inter$number*1.02))) +
  labs(title="Inter-interaction") +
  ylab('Number') +
  xlab("PET count per interaction") +
  theme_bw() +
  theme(panel.grid=element_blank(),axis.text.x = element_text(angle = 45, hjust = 1))

pdf(args[2], width=12)
print(p1 + p2)
dev.off()
