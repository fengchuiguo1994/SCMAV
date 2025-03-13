library(plotrix)
args=commandArgs(T) 
r=0.3

mydat = read.table(args[1])
total = as.numeric(gsub(",", "", mydat[mydat$V1 == "DNA_raw_reads",]$V2))
TwoTags = as.numeric(gsub(",", "", mydat[mydat$V1 == "DNA_raw_2tag_with_linker",]$V2))
OneTag = as.numeric(gsub(",", "", mydat[mydat$V1 == "DNA_raw_1tag_with_linker",]$V2))
NoLinker = as.numeric(gsub(",", "", mydat[mydat$V1 == "DNA_raw_without_linker",]$V2))
other = total - NoLinker - OneTag - TwoTags

TwoTags_UU = as.numeric(gsub(",", "", mydat[mydat$V1 == "PET_UU_all",]$V2))
# TwoTags_uu_hq = as.numeric(gsub(",", "", mydat[mydat$V1 == "PET_UU_mapped_q30",]$V2))
TwoTags_uu_hq = as.numeric(gsub(",", "", mydat[grep("PET_UU_mapped_q\\d+$", mydat$V1), ]$V2))
TwoTags_UU_low = TwoTags_UU - TwoTags_uu_hq
# TwoTags_uu_cb = as.numeric(gsub(",", "", mydat[mydat$V1 == "PET_UU_mapped_q30_validCB",]$V2))
TwoTags_uu_cb = as.numeric(gsub(",", "", mydat[grep("PET_UU_mapped_q\\d+_validCB$", mydat$V1),]$V2))
TwoTags_uu_nocb = TwoTags_uu_hq - TwoTags_uu_cb
# TwoTags_uu_nr = as.numeric(gsub(",", "", mydat[mydat$V1 == "PET_UU_mapped_q30_validCB_dedup",]$V2))
TwoTags_uu_nr = as.numeric(gsub(",", "", mydat[grep("PET_UU_mapped_q\\d+_validCB_dedup$", mydat$V1),]$V2))
TwoTags_uu_nonr = TwoTags_uu_cb - TwoTags_uu_nr
TwoTags_uu_nopet = TwoTags_UU_low + TwoTags_uu_nocb + TwoTags_uu_nonr
TwoTags_uu_trans = as.numeric(gsub(",", "", mydat[mydat$V1 == "PET_trans",]$V2))
TwoTags_uu_cis = as.numeric(gsub(",", "", mydat[mydat$V1 == "PET_cis",]$V2))
TwoTags_uu_cisgt1k = as.numeric(gsub(",", "", mydat[mydat$V1 == "PET_cis>1kb",]$V2))
TwoTags_uu_cisle1k = TwoTags_uu_cis - TwoTags_uu_cisgt1k
TwoTags_UxxU = as.numeric(gsub(",", "", mydat[mydat$V1 == "PET_UxxU_all",]$V2))
# TwoTags_uxxu_hq = as.numeric(gsub(",", "", mydat[mydat$V1 == "PET_UxxU_mapped_q30",]$V2))
TwoTags_uxxu_hq = as.numeric(gsub(",", "", mydat[grep("PET_UxxU_mapped_q\\d+$", mydat$V1),]$V2))
TwoTags_uxxu_low = TwoTags_UxxU - TwoTags_uxxu_hq
# TwoTags_uxxu_cb = as.numeric(gsub(",", "", mydat[mydat$V1 == "PET_UxxU_mapped_q30_validCB",]$V2))
TwoTags_uxxu_cb = as.numeric(gsub(",", "", mydat[grep("PET_UxxU_mapped_q\\d+_validCB$", mydat$V1),]$V2))
TwoTags_uxxu_nocb = TwoTags_uxxu_hq - TwoTags_uxxu_cb
# TwoTags_uxxu_nr = as.numeric(gsub(",", "", mydat[mydat$V1 == "PET_UxxU_mapped_q30_validCB_dedup",]$V2))
TwoTags_uxxu_nr = as.numeric(gsub(",", "", mydat[grep("PET_UxxU_mapped_q\\d+_validCB_dedup$", mydat$V1),]$V2))
TwoTags_uxxu_nonr = TwoTags_uxxu_cb - TwoTags_uxxu_nr
TwoTags_other = TwoTags - TwoTags_UU - TwoTags_UxxU
TwoTags_UxxU_nopet = TwoTags_uxxu_low + TwoTags_uxxu_nocb + TwoTags_uxxu_nonr + TwoTags_uxxu_nr + TwoTags_other


OneTag_UxxU = as.numeric(gsub(",", "", mydat[mydat$V1 == "single_UxxU_all",]$V2))
# OneTag_uxxu_hq = as.numeric(gsub(",", "", mydat[mydat$V1 == "single_UxxU_mapped_q30",]$V2))
OneTag_uxxu_hq = as.numeric(gsub(",", "", mydat[grep("single_UxxU_mapped_q\\d+$", mydat$V1),]$V2))
OneTag_uxxu_low = OneTag_UxxU - OneTag_uxxu_hq
# OneTag_uxxu_cb = as.numeric(gsub(",", "", mydat[mydat$V1 == "single_UxxU_mapped_q30_validCB",]$V2))
OneTag_uxxu_cb = as.numeric(gsub(",", "", mydat[grep("single_UxxU_mapped_q\\d+_validCB$", mydat$V1),]$V2))
OneTag_uxxu_nocb = OneTag_uxxu_hq - OneTag_uxxu_cb
# OneTag_uxxu_nr = as.numeric(gsub(",", "", mydat[mydat$V1 == "single_UxxU_mapped_q30_validCB_dedup",]$V2))
OneTag_uxxu_nr = as.numeric(gsub(",", "", mydat[grep("single_UxxU_mapped_q\\d+_validCB_dedup$", mydat$V1),]$V2))
OneTag_uxxu_nonr = OneTag_uxxu_cb - OneTag_uxxu_nr
OneTag_other = OneTag - OneTag_UxxU
OneTag_UxxU_nopet = OneTag

NoLinker_UU = as.numeric(gsub(",", "", mydat[mydat$V1 == "NoLinker_UU_all",]$V2))
# NoLinker_uu_hq = as.numeric(gsub(",", "", mydat[mydat$V1 == "NoLinker_UU_mapped_q30",]$V2))
NoLinker_uu_hq = as.numeric(gsub(",", "", mydat[grep("NoLinker_UU_mapped_q\\d+$", mydat$V1),]$V2))
NoLinker_UU_low = NoLinker_UU - NoLinker_uu_hq
# NoLinker_uu_cb = as.numeric(gsub(",", "", mydat[mydat$V1 == "NoLinker_UU_mapped_q30_validCB",]$V2))
NoLinker_uu_cb = as.numeric(gsub(",", "", mydat[grep("NoLinker_UU_mapped_q\\d+_validCB$", mydat$V1),]$V2))
NoLinker_uu_nocb = NoLinker_uu_hq - NoLinker_uu_cb
# NoLinker_uu_nr = as.numeric(gsub(",", "", mydat[mydat$V1 == "NoLinker_UU_mapped_q30_validCB_dedup",]$V2))
NoLinker_uu_nr = as.numeric(gsub(",", "", mydat[grep("NoLinker_UU_mapped_q\\d+_validCB_dedup$", mydat$V1),]$V2))
NoLinker_uu_nonr = NoLinker_uu_cb - NoLinker_uu_nr
NoLinker_uu_trans = as.numeric(gsub(",", "", mydat[mydat$V1 == "NoLinker_trans",]$V2))
NoLinker_uu_cisgt1k = as.numeric(gsub(",", "", mydat[mydat$V1 == "NoLinker_cis>1kb",]$V2))
NoLinker_uu_nopet = NoLinker_UU_low + NoLinker_uu_nocb + NoLinker_uu_nonr + (NoLinker_uu_nr - NoLinker_uu_trans- NoLinker_uu_cisgt1k)
NoLinker_UxxU = as.numeric(gsub(",", "", mydat[mydat$V1 == "NoLinker_UxxU_all",]$V2))
# NoLinker_uxxu_hq = as.numeric(gsub(",", "", mydat[mydat$V1 == "NoLinker_UxxU_mapped_q30",]$V2))
NoLinker_uxxu_hq = as.numeric(gsub(",", "", mydat[grep("NoLinker_UxxU_mapped_q\\d+$", mydat$V1),]$V2))
NoLinker_uxxu_low = NoLinker_UxxU - NoLinker_uxxu_hq
# NoLinker_uxxu_cb = as.numeric(gsub(",", "", mydat[mydat$V1 == "NoLinker_UxxU_mapped_q30_validCB",]$V2))
NoLinker_uxxu_cb = as.numeric(gsub(",", "", mydat[grep("NoLinker_UxxU_mapped_q\\d+_validCB$", mydat$V1),]$V2))
NoLinker_uxxu_nocb = NoLinker_uxxu_hq - NoLinker_uxxu_cb
# NoLinker_uxxu_nr = as.numeric(gsub(",", "", mydat[mydat$V1 == "NoLinker_UxxU_mapped_q30_validCB_dedup",]$V2))
NoLinker_uxxu_nr = as.numeric(gsub(",", "", mydat[grep("NoLinker_UxxU_mapped_q\\d+_validCB_dedup$", mydat$V1),]$V2))
NoLinker_uxxu_nonr = NoLinker_uxxu_cb - NoLinker_uxxu_nr
NoLinker_other = NoLinker - NoLinker_UU - NoLinker_UxxU
NoLinker_UxxU_nopet = NoLinker_UxxU + NoLinker_other

# perl -e '$aa="TwoTags_UU_low,TwoTags_uu_nocb,TwoTags_uu_nonr,TwoTags_uu_nr,TwoTags_uxxu_low,TwoTags_uxxu_nocb,TwoTags_uxxu_nonr,TwoTags_uxxu_nr,TwoTags_other,OneTag_uxxu_low,OneTag_uxxu_nocb,OneTag_uxxu_nonr,OneTag_uxxu_nr,OneTag_other,NoLinker_UU_low,NoLinker_uu_nocb,NoLinker_uu_nonr,NoLinker_uu_nr,NoLinker_uxxu_low,NoLinker_uxxu_nocb,NoLinker_uxxu_nonr,NoLinker_uxxu_nr,NoLinker_other"; @bb=split(/,/,$aa); for $cc (@bb){print "round($cc/total*100,2),"}'
# perl -e '$aa="TwoTags_UU_low,TwoTags_uu_nocb,TwoTags_uu_nonr,TwoTags_uu_nr,TwoTags_uxxu_low,TwoTags_uxxu_nocb,TwoTags_uxxu_nonr,TwoTags_uxxu_nr,TwoTags_other,OneTag_uxxu_low,OneTag_uxxu_nocb,OneTag_uxxu_nonr,OneTag_uxxu_nr,OneTag_other,NoLinker_UU_low,NoLinker_uu_nocb,NoLinker_uu_nonr,NoLinker_uu_nr,NoLinker_uxxu_low,NoLinker_uxxu_nocb,NoLinker_uxxu_nonr,NoLinker_uxxu_nr,NoLinker_other"; @bb=split(/,/,$aa); for $cc (@bb){print "$cc+"}'

pdf(args[2],width=10)
cols <- c(TwoTags="#00a0e9", OneTag="#da1c6f", NoLinker="#f8b62b",Other='white',
          lowQ="#4DBBD5FF", noCB="#E64B35FF", PCR="#8491B4FF", HQ="#E18727FF",TwoTagsOther='white',OneTagOther='white',NoLinkerOther='white',
          transPET="#ED0000FF", `cisPET<1kb`="#7876B1FF", `cisPET>=1kb`="#20854EFF", NoPET="white")
pie(1, radius=r, init.angle=90, col="white", border=NA, labels='')

nn = floating.pie(0,0, c(TwoTags_uu_nopet,TwoTags_uu_trans,TwoTags_uu_cisle1k,TwoTags_uu_cisgt1k,TwoTags_UxxU_nopet,
                         OneTag_UxxU_nopet,
                         NoLinker_uu_nopet,NoLinker_uu_trans,NoLinker_uu_cisgt1k,NoLinker_UxxU_nopet,
                         other),radius=3*r,startpos=pi/2,
                  col=cols[c("NoPET", "transPET",'cisPET<1kb',"cisPET>=1kb", "NoPET",  
                             "NoPET",
                             "NoPET", "transPET","cisPET>=1kb", "NoPET",  
                             "NoPET")],border=NA)
# id = paste(c(round(TwoTags_uu_nopet/total*100,2),round(TwoTags_uu_trans/total*100,2),round(TwoTags_uu_cisle1k/total*100,2),round(TwoTags_uu_cisgt1k/total*100,2),round(TwoTags_UxxU_nopet/total*100,2),round(OneTag_UxxU_nopet/total*100,2),round(NoLinker_uu_nopet/total*100,2),round(NoLinker_uu_trans/total*100,2),round(NoLinker_uu_cisgt1k/total*100,2),round(NoLinker_UxxU_nopet/total*100,2),round(other/total*100,2)),"%",sep="")
id = paste(c(round(TwoTags_uu_nopet/total*100,2),round(TwoTags_uu_trans/total*100,2),round(TwoTags_uu_cisle1k/total*100,2),round(TwoTags_uu_cisgt1k/total*100,2),round(TwoTags_UxxU_nopet/total*100,2),
             round(OneTag_UxxU_nopet/total*100,2),
             round(NoLinker_uu_nopet/total*100,2),round(NoLinker_uu_trans/total*100,2),round(NoLinker_uu_cisgt1k/total*100,2),round(NoLinker_UxxU_nopet/total*100,2),
             round(other/total*100,2)),"%",sep="")
id = paste(id,c("NoPET","transPET", "cisPET<1kb", "cisPET>=1kb","NoPET","NoPET", "NoPET", "transPET","cisPET>=1kb","NoPET","NoPET"),sep="\n")
pie.labels(0,0,nn,id,radius=2*r+0.08,minangle=0,explode=0.25,border=NA, cex = 0.5)

floating.pie(0,0, c(total),radius=2.5*r,startpos=pi/2,col=cols[c("Other")],border=NA)

nn = floating.pie(0,0, c(TwoTags_UU_low,TwoTags_uu_nocb,TwoTags_uu_nonr,TwoTags_uu_nr,
                         TwoTags_uxxu_low,TwoTags_uxxu_nocb,TwoTags_uxxu_nonr,TwoTags_uxxu_nr,TwoTags_other,
                         OneTag_uxxu_low,OneTag_uxxu_nocb,OneTag_uxxu_nonr,OneTag_uxxu_nr,OneTag_other,
                         NoLinker_UU_low,NoLinker_uu_nocb,NoLinker_uu_nonr,NoLinker_uu_nr,
                         NoLinker_uxxu_low,NoLinker_uxxu_nocb,NoLinker_uxxu_nonr,NoLinker_uxxu_nr,NoLinker_other,other),radius=2*r,startpos=pi/2,
                  col=cols[c("lowQ", "noCB",'PCR',"HQ", 
                             "lowQ", "noCB",'PCR',"HQ", 'TwoTagsOther',
                             "lowQ", "noCB",'PCR',"HQ", "OneTagOther",
                             "lowQ", "noCB",'PCR',"HQ",
                             "lowQ", "noCB",'PCR',"HQ", "NoLinkerOther", "Other")],border=NA)
# id = paste(c(round(TwoTags_UU_low/total*100,2),round(TwoTags_uu_nocb/total*100,2),round(TwoTags_uu_nonr/total*100,2),round(TwoTags_uu_nr/total*100,2),round(TwoTags_uxxu_low/total*100,2),round(TwoTags_uxxu_nocb/total*100,2),round(TwoTags_uxxu_nonr/total*100,2),round(TwoTags_uxxu_nr/total*100,2),round(TwoTags_other/total*100,2),round(OneTag_uxxu_low/total*100,2),round(OneTag_uxxu_nocb/total*100,2),round(OneTag_uxxu_nonr/total*100,2),round(OneTag_uxxu_nr/total*100,2),round(OneTag_other/total*100,2),round(NoLinker_UU_low/total*100,2),round(NoLinker_uu_nocb/total*100,2),round(NoLinker_uu_nonr/total*100,2),round(NoLinker_uu_nr/total*100,2),round(NoLinker_uxxu_low/total*100,2),round(NoLinker_uxxu_nocb/total*100,2),round(NoLinker_uxxu_nonr/total*100,2),round(NoLinker_uxxu_nr/total*100,2),round(NoLinker_other/total*100,2),round(other/total*100,2)),"%",sep="")
id = paste(c(round(TwoTags_UU_low/total*100,2),round(TwoTags_uu_nocb/total*100,2),round(TwoTags_uu_nonr/total*100,2),round(TwoTags_uu_nr/total*100,2),
             round(TwoTags_uxxu_low/total*100,2),round(TwoTags_uxxu_nocb/total*100,2),round(TwoTags_uxxu_nonr/total*100,2),round(TwoTags_uxxu_nr/total*100,2),
             round(TwoTags_other/total*100,2),round(OneTag_uxxu_low/total*100,2),round(OneTag_uxxu_nocb/total*100,2),round(OneTag_uxxu_nonr/total*100,2),
             round(OneTag_uxxu_nr/total*100,2),round(OneTag_other/total*100,2),round(NoLinker_UU_low/total*100,2),round(NoLinker_uu_nocb/total*100,2),
             round(NoLinker_uu_nonr/total*100,2),round(NoLinker_uu_nr/total*100,2),round(NoLinker_uxxu_low/total*100,2),round(NoLinker_uxxu_nocb/total*100,2),
             round(NoLinker_uxxu_nonr/total*100,2),round(NoLinker_uxxu_nr/total*100,2),round(NoLinker_other/total*100,2),round(other/total*100,2)),"%",sep="")
id = paste(id,c("UU_lowQ","UU_noCB", "UU_PCR", "UU_HQ","UxxU_lowQ","UxxU_noCB", "UxxU_PCR", "UxxU_HQ","TwoTagsOther","UxxU_lowQ","UxxU_noCB", "UxxU_PCR", "UxxU_HQ","OneTagOther","UU_lowQ","UU_noCB", "UU_PCR", "UU_HQ","UxxU_lowQ","UxxU_noCB", "UxxU_PCR", "UxxU_HQ","NoLinkerOther"),sep="\n")
pie.labels(0,0,nn,id,radius=1*r+0.08,minangle=0,explode=0.25,border=NA, cex = 0.5)

floating.pie(0,0, c(total),radius=1.5*r,startpos=pi/2,col=cols[c("Other")],border=NA)

nn = floating.pie(0, 0, c(TwoTags, OneTag, NoLinker, other),radius=1*r,startpos=pi/2,col=cols[c("TwoTags", "OneTag", "NoLinker", "Other")],border=NA)
id = paste(c(round(TwoTags/total*100,2),round(OneTag/total*100,2),round(NoLinker/total*100,2),round(other/total*100,2)),"%",sep="")
id = paste(id,c("TwoTags","OneTag", "NoLinker", "Other"),sep="\n")
pie.labels(0,0,nn,id,radius=1*r-0.05,minangle=0,explode=0.05,border=NA, cex = 0.5)
legend(3.5*r, 3.5*r, legend=sub("_", " ", names(cols)[-1]),fill=cols[-1], bty="n", cex=1.2)




pie(1, radius=r, init.angle=90, col="white", border=NA, labels='')

nn = floating.pie(0,0, c(TwoTags_uu_nopet,TwoTags_uu_trans,TwoTags_uu_cisle1k,TwoTags_uu_cisgt1k,TwoTags_UxxU_nopet,
                         OneTag_UxxU_nopet,
                         NoLinker_uu_nopet,NoLinker_uu_trans,NoLinker_uu_cisgt1k,NoLinker_UxxU_nopet,
                         other),radius=3*r,startpos=pi/2,
                  col=cols[c("NoPET", "transPET",'cisPET<1kb',"cisPET>=1kb", "NoPET",  
                             "NoPET",
                             "NoPET", "transPET","cisPET>=1kb", "NoPET",  
                             "NoPET")],border=NA)
# id = paste(c(round(TwoTags_uu_nopet/total*100,2),round(TwoTags_uu_trans/total*100,2),round(TwoTags_uu_cisle1k/total*100,2),round(TwoTags_uu_cisgt1k/total*100,2),round(TwoTags_UxxU_nopet/total*100,2),round(OneTag_UxxU_nopet/total*100,2),round(NoLinker_uu_nopet/total*100,2),round(NoLinker_uu_trans/total*100,2),round(NoLinker_uu_cisgt1k/total*100,2),round(NoLinker_UxxU_nopet/total*100,2),round(other/total*100,2)),"%",sep="")
id = paste(c(round(TwoTags_uu_nopet/total*100,2),round(TwoTags_uu_trans/total*100,2),round(TwoTags_uu_cisle1k/total*100,2),round(TwoTags_uu_cisgt1k/total*100,2),round(TwoTags_UxxU_nopet/total*100,2),
             round(OneTag_UxxU_nopet/total*100,2),
             round(NoLinker_uu_nopet/total*100,2),round(NoLinker_uu_trans/total*100,2),round(NoLinker_uu_cisgt1k/total*100,2),round(NoLinker_UxxU_nopet/total*100,2),
             round(other/total*100,2)),"%",sep="")
# id = paste(id,c("NoPET","transPET", "cisPET<1kb", "cisPET>=1kb","NoPET","NoPET", "NoPET", "transPET","cisPET>=1kb","NoPET","NoPET"),sep="\n")
id[c(1,5,6,7,10,11)] = ""
pie.labels(0,0,nn,id,radius=2*r+0.08,minangle=0,explode=0.25,border=NA, cex = 0.5)

floating.pie(0,0, c(total),radius=2.5*r,startpos=pi/2,col=cols[c("Other")],border=NA)

nn = floating.pie(0,0, c(TwoTags_UU_low,TwoTags_uu_nocb,TwoTags_uu_nonr,TwoTags_uu_nr,
                         TwoTags_uxxu_low,TwoTags_uxxu_nocb,TwoTags_uxxu_nonr,TwoTags_uxxu_nr,TwoTags_other,
                         OneTag_uxxu_low,OneTag_uxxu_nocb,OneTag_uxxu_nonr,OneTag_uxxu_nr,OneTag_other,
                         NoLinker_UU_low,NoLinker_uu_nocb,NoLinker_uu_nonr,NoLinker_uu_nr,
                         NoLinker_uxxu_low,NoLinker_uxxu_nocb,NoLinker_uxxu_nonr,NoLinker_uxxu_nr,NoLinker_other,other),radius=2*r,startpos=pi/2,
                  col=cols[c("lowQ", "noCB",'PCR',"HQ", 
                             "lowQ", "noCB",'PCR',"HQ", 'TwoTagsOther',
                             "lowQ", "noCB",'PCR',"HQ", "OneTagOther",
                             "lowQ", "noCB",'PCR',"HQ",
                             "lowQ", "noCB",'PCR',"HQ", "NoLinkerOther", "Other")],border=NA)
# id = paste(c(round(TwoTags_UU_low/total*100,2),round(TwoTags_uu_nocb/total*100,2),round(TwoTags_uu_nonr/total*100,2),round(TwoTags_uu_nr/total*100,2),round(TwoTags_uxxu_low/total*100,2),round(TwoTags_uxxu_nocb/total*100,2),round(TwoTags_uxxu_nonr/total*100,2),round(TwoTags_uxxu_nr/total*100,2),round(TwoTags_other/total*100,2),round(OneTag_uxxu_low/total*100,2),round(OneTag_uxxu_nocb/total*100,2),round(OneTag_uxxu_nonr/total*100,2),round(OneTag_uxxu_nr/total*100,2),round(OneTag_other/total*100,2),round(NoLinker_UU_low/total*100,2),round(NoLinker_uu_nocb/total*100,2),round(NoLinker_uu_nonr/total*100,2),round(NoLinker_uu_nr/total*100,2),round(NoLinker_uxxu_low/total*100,2),round(NoLinker_uxxu_nocb/total*100,2),round(NoLinker_uxxu_nonr/total*100,2),round(NoLinker_uxxu_nr/total*100,2),round(NoLinker_other/total*100,2),round(other/total*100,2)),"%",sep="")
id = paste(c(round(TwoTags_UU_low/total*100,2),round(TwoTags_uu_nocb/total*100,2),round(TwoTags_uu_nonr/total*100,2),round(TwoTags_uu_nr/total*100,2),
             round(TwoTags_uxxu_low/total*100,2),round(TwoTags_uxxu_nocb/total*100,2),round(TwoTags_uxxu_nonr/total*100,2),round(TwoTags_uxxu_nr/total*100,2),
             round(TwoTags_other/total*100,2),round(OneTag_uxxu_low/total*100,2),round(OneTag_uxxu_nocb/total*100,2),round(OneTag_uxxu_nonr/total*100,2),
             round(OneTag_uxxu_nr/total*100,2),round(OneTag_other/total*100,2),round(NoLinker_UU_low/total*100,2),round(NoLinker_uu_nocb/total*100,2),
             round(NoLinker_uu_nonr/total*100,2),round(NoLinker_uu_nr/total*100,2),round(NoLinker_uxxu_low/total*100,2),round(NoLinker_uxxu_nocb/total*100,2),
             round(NoLinker_uxxu_nonr/total*100,2),round(NoLinker_uxxu_nr/total*100,2),round(NoLinker_other/total*100,2),round(other/total*100,2)),"%",sep="")
id = paste(id,c("UU_lowQ","UU_noCB", "UU_PCR", "UU_HQ","UxxU_lowQ","UxxU_noCB", "UxxU_PCR", "UxxU_HQ","TwoTagsOther","UxxU_lowQ","UxxU_noCB", "UxxU_PCR", "UxxU_HQ","OneTagOther","UU_lowQ","UU_noCB", "UU_PCR", "UU_HQ","UxxU_lowQ","UxxU_noCB", "UxxU_PCR", "UxxU_HQ","NoLinkerOther"),sep="\n")
id[c(1:3,5:7,9,10:12,14,15:17,19:21,23)] = ""
pie.labels(0,0,nn,id,radius=1*r+0.08,minangle=0,explode=0.25,border=NA, cex = 0.5)

floating.pie(0,0, c(total),radius=1.5*r,startpos=pi/2,col=cols[c("Other")],border=NA)

nn = floating.pie(0, 0, c(TwoTags, OneTag, NoLinker, other),radius=1*r,startpos=pi/2,col=cols[c("TwoTags", "OneTag", "NoLinker", "Other")],border=NA)
id = paste(c(round(TwoTags/total*100,2),round(OneTag/total*100,2),round(NoLinker/total*100,2),round(other/total*100,2)),"%",sep="")
id = paste(id,c("TwoTags","OneTag", "NoLinker", "Other"),sep="\n")
pie.labels(0,0,nn,id,radius=1*r-0.05,minangle=0,explode=0.05,border=NA, cex = 0.5)
legend(3.5*r, 3.5*r, legend=sub("_", " ", names(cols)[-1]),fill=cols[-1], bty="n", cex=1.2)
dev.off()
