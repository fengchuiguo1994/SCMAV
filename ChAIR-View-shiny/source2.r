library(DBI)
library(RMySQL)
library(RColorBrewer)
library(plotrix)
library(pheatmap)
library(strawr)
library(dplyr)
library(GenomicRanges)

options(scipen = 200)


run2 <- function(region="chr10:3600000-6600000",db="hg38",tbl="MM10.RNA.COVER",petCount=1,adjustBinNumber=T,maxCellNumber=800,mode="all",downsample=-1,gene="",shareBarcode="Don't share barcode",anno="",showLabels=T,trackInfo=c(),parfig = c(0,1,0,1),plotnew = F,hic_part=1/4,filterDis=8000){
  
  ############################
  # input parameter
  ############################
  
  region = gsub("-|:"," ",region)
  regionarr = strsplit(trimws(region)," |\t")
  chrom = regionarr[[1]][1]
  start = as.numeric(regionarr[[1]][2])
  end = as.numeric(regionarr[[1]][3])

  read_list <- readData(chrom,start,end,tbl,db,downsample,filterDis=filterDis)
  result <- read_list[[1]]
  format <- read_list[[2]]
  nrecords <- read_list[[3]]
  
  direction = F
  directions = F
  color_class=c()
  transcriptinfo=""
  
  ##############################
  # single cell barcode order
  ##############################
  #print(proc.time())
  
  if(nrow(result)!=0){
    
    if(shareBarcode == "Share Barcode And Sort By PET"&&(format=="scPET"||format=="scRNA"||format=="scATAC")){
      con <- dbConnect(MySQL(), username = 'root', password = 'kkltmax8', dbname = db,host = 'localhost')
      group <- dbGetQuery(con,paste0("SELECT `group` FROM header WHERE sample = '",tbl,"'"))
      if(!is.na(group)){
        tbl2 <- dbGetQuery(con,paste0("SELECT sample FROM header WHERE `group` = '",group,"' AND format = 'scPET'"))
        if(!is.na(tbl2)){
          read_list2 <- readData(chrom,start,end,tbl2,db,downsample)
          result2 <- read_list2[[1]]
          format2 <- read_list2[[2]]
          ### limit max cell number
          if(format2 =="scRNA"||format2 == "scATAC"){
            order <- data.frame(barcode=names(table(result2$name)),barcodecount=as.vector(table(result2$name)))
            order <- order[order(-order$barcodecount),]
            if(nrow(order)>maxCellNumber) {
              order <- order[1:maxCellNumber,]
              result2 <- result2[result2$name %in% order$barcode,]
            }
          }else if(format2 == "scPET"){
            order <- data.frame(barcode=names(table(result2$barcode)),barcodecount=as.vector(table(result2$barcode)))
            order <- order[order(-order$barcodecount),]
            if(nrow(order)>maxCellNumber) {
              order <- order[1:maxCellNumber,]
              result2 <- result2[result2$barcode %in% order$barcode,]
            }
          }

          if(mode=="TargetOneGene"&&format2=="scPET"&&gene!="") direction = T
          if(mode=="TargetGenes"&&format2=="scPET") directions = T
          
          if(format2=="scPET"){
            if(direction){
              transcriptinfo <- readAnno(chrom=chrom,start=start,end=end,anno=anno,db=db)
              result2 <- filterData(result=result2,format=format,mode=mode,gene=gene,transcriptinfo=transcriptinfo)
              geneinfo <- list(transcriptinfo[gene,]$starts,transcriptinfo[gene,]$stops,transcriptinfo[gene,]$strand)
              tmp_list <- sortedByDirection(result=result2,geneinfo = geneinfo)
              orderBarcode <- tmp_list[[1]]$orderBarcode
              color_class <- tmp_list[[1]]$color_class
              result2 <- tmp_list[[2]]
              if(format=="scPET") result <- result2

            }else if(directions){
              transcriptinfo <- readAnno(chrom=chrom,start=start,end=end,anno=anno,db=db)
              result2 <- filterData(result=result2,format=format,mode=mode,gene=gene,transcriptinfo=transcriptinfo)
              geneinfo <- list(transcriptinfo[gene,]$starts,transcriptinfo[gene,]$stops,transcriptinfo[gene,]$strand)
              tmp_dataframe <- clusterPET(result=result2,start=start,end=end)
              orderBarcode <- tmp_dataframe$orderBarcode
              color_class <- tmp_dataframe$color_class
            }else{
              tmp_dataframe <- clusterPET(result=result2,start=start,end=end)
              orderBarcode <- tmp_dataframe$orderBarcode
              color_class <- tmp_dataframe$color_class
            }
          }        
        }
      }
      dbDisconnect(con)
    }else if((shareBarcode == "Share Barcode And Sort By ATAC")&&(format=="scPET"||format=="scRNA"||format=="scATAC")){
      con <- dbConnect(MySQL(), username = 'root', password = 'kkltmax8', dbname = db,host = 'localhost')
      group <- dbGetQuery(con,paste0("SELECT `group` FROM header WHERE sample = '",tbl,"'"))
      if(!is.na(group)){
        tbl2 <- dbGetQuery(con,paste0("SELECT sample FROM header WHERE `group` = '",group,"' AND format = 'scATAC'"))
        if(!is.na(tbl2)){
          read_list2 <- readData(chrom,start,end,tbl2,db,downsample)
          result2 <- read_list2[[1]]
          format2 <- read_list2[[2]]
          ### limit max cell number
          if(format2 =="scRNA"||format2 == "scATAC"){
            order <- data.frame(barcode=names(table(result2$name)),barcodecount=as.vector(table(result2$name)))
            order <- order[order(-order$barcodecount),]
              if(nrow(order)>maxCellNumber) {
                order <- order[1:maxCellNumber,]
                result2 <- result2[result2$name %in% order$barcode,]
              }
          }else if(format2 == "scPET"){
            order <- data.frame(barcode=names(table(result2$barcode)),barcodecount=as.vector(table(result2$barcode)))
            order <- order[order(-order$barcodecount),]
              if(nrow(order)>maxCellNumber) {
                order <- order[1:maxCellNumber,]
                result2 <- result2[result2$barcode %in% order$barcode,]
              }
          }
          if(format2=="scATAC"){
            orderBarcode <- sortedByNumber(result2)
          }    
        }
      }
      dbDisconnect(con)
    }else if((shareBarcode == "Share Barcode And Sort By RNA")&&(format=="scPET"||format=="scRNA"||format=="scATAC")){
      con <- dbConnect(MySQL(), username = 'root', password = 'kkltmax8', dbname = db,host = 'localhost')
      group <- dbGetQuery(con,paste0("SELECT `group` FROM header WHERE sample = '",tbl,"'"))
      if(!is.na(group)){
        tbl2 <- dbGetQuery(con,paste0("SELECT sample FROM header WHERE `group` = '",group,"' AND format = 'scRNA'"))
        if(!is.na(tbl2)){
          read_list2 <- readData(chrom,start,end,tbl2,db,downsample)
          result2 <- read_list2[[1]]
          format2 <- read_list2[[2]]
          ### limit max cell number
          if(format2 =="scRNA"||format2 == "scATAC"){
            order <- data.frame(barcode=names(table(result2$name)),barcodecount=as.vector(table(result2$name)))
            order <- order[order(-order$barcodecount),]
              if(nrow(order)>maxCellNumber) {
                order <- order[1:maxCellNumber,]
                result2 <- result2[result2$name %in% order$barcode,]
              }
          }else if(format2 == "scPET"){
            order <- data.frame(barcode=names(table(result2$barcode)),barcodecount=as.vector(table(result2$barcode)))
            order <- order[order(-order$barcodecount),]
              if(nrow(order)>maxCellNumber) {
                order <- order[1:maxCellNumber,]
                result2 <- result2[result2$barcode %in% order$barcode,]
              }
          }
          if(format2=="scRNA"){
            orderBarcode <- sortedByNumber(result2)
          }    
        }
      }
      dbDisconnect(con)
    }
    
    if(shareBarcode == "Don't Share Barcode"){
        ### limit max cell number
      if(format =="scRNA"||format == "scATAC"){
        order <- data.frame(barcode=names(table(result$name)),barcodecount=as.vector(table(result$name)))
        order <- order[order(-order$barcodecount),]
          if(nrow(order)>maxCellNumber) {
            order <- order[1:maxCellNumber,]
            result <- result[result$name %in% order$barcode,]
          }        
      }else if(format == "scPET"){
        order <- data.frame(barcode=names(table(result$barcode)),barcodecount=as.vector(table(result$barcode)))
        order <- order[order(-order$barcodecount),]
          if(nrow(order)>maxCellNumber) {
            order <- order[1:maxCellNumber,]
            result <- result[result$barcode %in% order$barcode,]
          }
      }else if(format == "singleton"){
        order <- data.frame(barcode=names(table(result$name)),barcodecount=as.vector(table(result$name)))
        order <- order[order(-order$barcodecount),]
        
        order <- order[order$barcodecount>1,]
        result <- result[result$name %in% order$barcode,]
        if(nrow(order)>maxCellNumber) {
          order <- order[1:maxCellNumber,]
          result <- result[result$name %in% order$barcode,]
        }        
      }
      
      if(mode=="TargetOneGene"&&format=="scPET"&&gene!="") direction = T
      if(mode=="TargetGenes"&&format=="scPET") {directions = T}

      if(format=="scPET"){
        if(direction){
          transcriptinfo <- readAnno(chrom=chrom,start=start,end=end,anno=anno,db=db)
          result <- filterData(result=result,format=format,mode=mode,gene=gene,transcriptinfo=transcriptinfo)
          geneinfo <- list(transcriptinfo[gene,]$starts,transcriptinfo[gene,]$stops,transcriptinfo[gene,]$strand)
          tmp_list <- sortedByDirection(result=result,geneinfo = geneinfo)
          orderBarcode <- tmp_list[[1]]$orderBarcode
          color_class <- tmp_list[[1]]$color_class
          result <- tmp_list[[2]]
        }else if(directions){
          transcriptinfo <- readAnno(chrom=chrom,start=start,end=end,anno=anno,db=db)
          result <- filterData(result=result,format=format,mode=mode,gene=gene,transcriptinfo=transcriptinfo)
          geneinfo <- list(transcriptinfo[gene,]$starts,transcriptinfo[gene,]$stops,transcriptinfo[gene,]$strand)
          tmp_dataframe <- clusterPET(result=result,start=start,end=end)
          orderBarcode <- tmp_dataframe$orderBarcode
          color_class <- tmp_dataframe$color_class
        }else{
          tmp_dataframe <- clusterPET(result=result,start=start,end=end)
          orderBarcode <- tmp_dataframe$orderBarcode
          color_class <- tmp_dataframe$color_class
          print(orderBarcode)
        }
      } else if(format=="scRNA" || format=="scATAC"){
        orderBarcode <- sortedByNumber(result)
      } else if(format=="singleton"){
        tmp_dataframe <- clusterSingleton(result=result,start=start,end=end)
        orderBarcode <- tmp_dataframe$orderBarcode
      }
    }
  }else{
    orderBarcode <- c()
  }
  
  
  ##############
  # Plot
  ##############
  #print(proc.time())
  par(cex=0.5, mai=c(0,0.22,0,0.05))

  if(plotnew){
    par(fig=parfig,new=T)
  }else{
    par(fig=parfig)  
  }
  
  min_vec = trackInfo$min_vec
  max_vec = trackInfo$max_vec
  color1_vec = trackInfo$color1_vec
  color2_vec = trackInfo$color2_vec
  norm_vec <- trackInfo$norm_vec
  res_vec <- trackInfo$res_vec
  
  if(format=="bedgraph"||format=="RNAbedgraph"||format=="ATACbedgraph"){
    if(nrow(result)!=0){
      if(adjustBinNumber) result <- adjustBin(result,format,chrom,start,end)
      plot(c(0,0),type='n',ylab="",xlim=c(start,end),ylim=c(min_vec,max_vec),xaxt="n",xaxs="i",yaxs="i")
      lapply(1:nrow(result), function(row){
        xleft = result$start[row]
        ybottom = 0
        xright = result$end[row]
        ytop = result$value[row]
        #color = "blue"
        # color = rgb(21/255,29/255,252/255)
        color = color1_vec
        if(xleft < start) xleft = start
        if(xright > end) xright = end
        rect(xleft, ybottom, xright, ytop,col=color,border=color)
      })
    }else{
      plot(c(0,0),type='n',xlim=c(start,end),ylim=c(0,1),xaxt="n",xaxs="i",yaxs="i")  
    }
    if(showLabels)  mtext(paste0(tbl), side = 3,cex=0.8,line=-1.5,adj=0)
  } else if(format == "ssbedgraph"){
    if(nrow(result)!=0){
      if(adjustBinNumber) result <- adjustBin(result,format,chrom,start,end)
      plot(c(0,0),type='n',ylab="",xlim=c(start,end),ylim=c(min_vec,max_vec),xaxt="n",xaxs="i",yaxs="i")
      lapply(1:nrow(result), function(row){
        xleft = result$start[row]
        ybottom = 0
        xright = result$end[row]
        ytop1 = result$value1[row]
        ytop2 = result$value2[row]
        color1 = color1_vec
        color2 = color2_vec
        if(xleft < start) xleft = start
        if(xright > end) xright = end
        rect(xleft, ybottom, xright, ytop1,col=color1,border=color1)
        rect(xleft, ytop2, xright, ybottom,col=color2,border=color2)
      })
    }else{
      plot(c(0,0),type='n',xlim=c(start,end),ylim=c(-1,1),xaxt="n",xaxs="i",yaxs="i")   
    }
    if(showLabels)  mtext(paste0(tbl), side = 3,cex=0.8,line=-1.5,adj=0)
  } else if(format == "chromhmm"||format == "bed6"){
    plot(c(0,0),type='n',ylab="",xlim=c(start,end), ylim=c(0,1),xaxt="n",xaxs="i",yaxs="i")
    if(nrow(result) != 0){
      lapply(1:nrow(result), function(row){
        if(result$color[row]!="-"){
          color_chromhmm = colorts(result$color[row])
        }else{
          color_chromhmm = "red"
        }
        rect(result$start[row], 0, result$end[row], 1,col=color_chromhmm,border=color_chromhmm)       
      })
    }
    if(showLabels)  mtext(paste0(tbl), side = 3,cex=0.8,line=-1.5,adj=0)
  } else if(format == "loop"){
    if (nrow(result) != 0){
      result <- result[result$value>=as.numeric(petCount),]
      plot(c(0,0),type='n',ylab="",xlim=c(start,end),ylim=c(min_vec,max_vec),xaxt="n",xaxs="i",yaxs="i")
      # result$value = log10(result$value + 1)
      result$value = result$value
      lapply(1:nrow(result), function(row){
        x1 = (result$start1[row] + result$end1[row])/2
        x2 = (result$start2[row] + result$end2[row])/2
        vert_y  = result$value[row]
        vert_x= (x1 + x2) / 2
        a = vert_y / ((vert_x - x1) * (vert_x - x2))
        color = color1_vec
        if(x1 < start) xleft = start
        if(x2 > end) xright = end
        curve(a * (x - x1) * (x - x2), from=x1, to=x2, add=TRUE, col=color,lwd = 1.5)
      })
    }else{
      plot(c(0,0),type='n',xlim=c(start,end),ylim=c(0,1),xaxt="n",xaxs="i",yaxs="i")
    }
    if(showLabels)  mtext(paste0(tbl), side = 3,cex=0.8,line=-1.5,adj=0)
  } else if(format == "anno"){
    #########################
    ### process plotrow 
    #########################
    numberofgeneinfo = length(names(table(result$geneName)))
    namesgeneinfo = names(table(result$geneName))
    starts = c()
    stops  = c()
    sizes  = c()
    strands = c()
    for (i in (1:numberofgeneinfo)) {
      subgeneinfo  = result[which(result$geneName == namesgeneinfo[i]),]
      starts = c(starts,min(subgeneinfo$start,subgeneinfo$end))
      stops  = c(stops, max(subgeneinfo$start,subgeneinfo$end))
      sizes  = c(sizes, stops[i] - starts[i])
      strands = c(strands,subgeneinfo$strand[1])
    }
    
    transcriptinfo = data.frame(names=namesgeneinfo,starts=starts,stops=stops,sizes=sizes,strand=strands)
    transcriptinfo = transcriptinfo[order(transcriptinfo$names),]
    typeshow = "squished"
    if (typeshow == "collapsed") {
      transcriptinfo$plotrow = 1
    } else if(typeshow == "squished") {
      if (nrow(transcriptinfo) > 0) {
        maxn = 1
        transcriptinfo$plotrow = 0
        for (i in (1:nrow(transcriptinfo))){
          overlapflag = T
          for (j in (1:maxn)) {
            transcriptinfotmp = transcriptinfo[transcriptinfo$plotrow==j &
                                                 transcriptinfo$starts[i] < transcriptinfo$stops & transcriptinfo$stops[i] > transcriptinfo$starts,]
            if (nrow(transcriptinfotmp) == 0) {
              transcriptinfo$plotrow[i] = j
              overlapflag = F
              break
            }
          }
          # update this
          if (overlapflag) {
            maxn = maxn + 1
            transcriptinfo$plotrow[i] = maxn
          }
        }
        # if (overlapflag) {
        #   maxn = maxn + 1
        #   transcriptinfo$plotrow[i] = maxn
        # }
      }
    } else if (typeshow == "flatten") {
      transcriptinfo$plotrow = seq(1:nrow(transcriptinfo))
    }
    
    #########################
    ### start plot 
    #########################

    plot(c(1,1),type='n',ylab="",xlim=c(start,end), ylim=c(-0.5,max(transcriptinfo$plotrow)+1),xaxt="n",xaxs="i",yaxs="i")
    if (nrow(transcriptinfo) != 0){
      lapply(1:nrow(transcriptinfo), function(row){
        geneName = transcriptinfo$names[row]
        subgeneinfo  = result[which(result$geneName == geneName),]
        if (transcriptinfo$strand[row]=="+"){
          #color = "green"
          color = rgb(23/255,128/255,15/255)
        } else {
          color = rgb(21/255,29/255,252/255)
          #color = "blue"
        }
        rect(subgeneinfo$start, transcriptinfo$plotrow[row] + 0, subgeneinfo$end, transcriptinfo$plotrow[row] + 0 + 0.5 ,col=color,border=color)
        segments(transcriptinfo$starts[row], transcriptinfo$plotrow[row] + 0.25 , transcriptinfo$stops[row], transcriptinfo$plotrow[row] + 0.25 ,col=color)
      })
      transcriptinfo$gene <- transcriptinfo$names
     genes <- unique(transcriptinfo$gene)
      lapply(1:length(genes), function(i){
        if (transcriptinfo[transcriptinfo$gene %in% genes[i],]$strand[1]=="+"){
          color = rgb(23/255,128/255,15/255)
        } else {
          color = rgb(21/255,29/255,252/255)
        }
        textrow <- min(transcriptinfo[transcriptinfo$gene %in% genes[i],]$plotrow)
        maxstop <- max(transcriptinfo[transcriptinfo$gene %in% genes[i],]$stops)
        minstart <- min(transcriptinfo[transcriptinfo$gene %in% genes[i],]$starts)
        text((maxstop+minstart)/2,textrow-0.25,labels=genes[i],col=color,cex=1.5)
      })
    }
  } else if(format == "scRNA"){
    cellNumber = length(orderBarcode)
    adjust_pos <-  0
    if(max(cellNumber) >= 150){
      point_ratio = 2
    }else if(max(cellNumber) >=50){
      point_ratio <- 150 / max(cellNumber) * 2 
    }else{
      point_ratio = 6
    }

    plot(c(1,1),type='n',ylab="",xlim=c(start,end), ylim=c(0,cellNumber+1),xaxt="n",xaxs="i")
    if(length(orderBarcode)==0){
      upvalue = 1
    }else{
      tmpdata = result[result$name %in% orderBarcode,]
      upvalue = floor(quantile(tmpdata$score, 0.99))    
    }

    color = colorRampPalette(brewer.pal(9,"Blues"))(upvalue+1)
    
    if(length(orderBarcode)!=0){
      lapply(1:length(orderBarcode), function(i){
        tmpdata = result[result$name==orderBarcode[i],]
        if (nrow(tmpdata) != 0){
          lapply(1:nrow(tmpdata), function(row){
            if (tmpdata$end[row] - tmpdata$start[row] < 150) {
              xleft = tmpdata$start[row]
              ybottom = 1*i + 0.25 + adjust_pos
              xright = tmpdata$end[row]
              ytop = 1*i + 0.75 + adjust_pos
              if(tmpdata$score[row]+4 > upvalue+1){
                color_count=upvalue+1
              }
              else{
                color_count=tmpdata$score[row]+4
              }
              points((xleft+xright)/2,(ybottom+ytop)/2,lwd=point_ratio * 0.15 / 0.4,cex=0.45,col=rgb(col2rgb(color[color_count])[1]/255,col2rgb(color[color_count])[2]/255,col2rgb(color[color_count])[3]/255,alpha = 0.8),bg=rgb(col2rgb(color[color_count])[1]/255,col2rgb(color[color_count])[2]/255,col2rgb(color[color_count])[3]/255,alpha = 0.8),pch=19)
            } else {
              xleft = tmpdata$start[row]
              ybottom = 1*i + 0.25 + adjust_pos
              xright = tmpdata$end[row]
              ytop = 1*i + 0.75 + adjust_pos
              if(tmpdata$score[row]+4 > upvalue+1){
                color_count=upvalue+1
              }
              else{
                color_count=tmpdata$score[row]+4
              }
              points((xleft+xleft+80)/2,(ybottom+ytop)/2,lwd=point_ratio * 0.15 / 0.4,cex=0.45,col=rgb(col2rgb(color[color_count])[1]/255,col2rgb(color[color_count])[2]/255,col2rgb(color[color_count])[3]/255,alpha = 0.8),bg=rgb(col2rgb(color[color_count])[1]/255,col2rgb(color[color_count])[2]/255,col2rgb(color[color_count])[3]/255,alpha = 0.8),pch=19)
              points((xright+xright+80)/2,(ybottom+ytop)/2,lwd=point_ratio * 0.15 / 0.4,cex=0.45,col=rgb(col2rgb(color[color_count])[1]/255,col2rgb(color[color_count])[2]/255,col2rgb(color[color_count])[3]/255,alpha = 0.8),bg=rgb(col2rgb(color[color_count])[1]/255,col2rgb(color[color_count])[2]/255,col2rgb(color[color_count])[3]/255,alpha = 0.8),pch=19)
            }
          })
        }
      })
    }
    if(showLabels)  mtext(paste0(tbl,"-",cellNumber,"/",nrecords," cells"), side = 3,cex=0.8,line=-1.5,adj=0)
    legend("topright", legend = c(0,"","","",upvalue), fill = brewer.pal(9,"Blues")[c(1,3,5,7,9)], bty = "n",border = brewer.pal(9,"Blues")[c(1,3,5,7,9)],y.intersp = 0.2,cex = 1.8)
  } else if(format == "scATAC"){
    cellNumber = length(orderBarcode)
    adjust_pos <-  0
    if(max(cellNumber) >= 150){
      point_ratio = 2
    }else if(max(cellNumber) >=50){
      point_ratio <- 150 / max(cellNumber) * 2 
    }else{
      point_ratio = 6
    }
    
    plot(c(1,1),type='n',ylab="",xlim=c(start,end), ylim=c(0,cellNumber+1),xaxt="n",xaxs="i")
    if(length(orderBarcode)!=0){
      lapply(1:length(orderBarcode), function(i){
        tmpdata = result[result$name==orderBarcode[i],]
        if (nrow(tmpdata) != 0){
          lapply(1:nrow(tmpdata), function(row){
            xleft = tmpdata$start[row]
            ybottom = 1*i - 0 + adjust_pos
            xright = tmpdata$end[row]
            ytop = 1*i + 0.05 + adjust_pos
            color = "black"
            if (tmpdata$score[row] == 1){
              color = "brown"
            }
            points((xleft+xright)/2,(ybottom+ytop)/2,lwd=point_ratio * 0.15 / 0.4,cex=0.45,col=color[tmpdata$score[row]],bg=color[tmpdata$score[row]],pch=19)
          })
        }
      })
    }
    if(showLabels)  mtext(paste0(tbl,"-",cellNumber,"/",nrecords," cells"), side = 3,cex=0.8,line=-1.5,adj=0)
    legend("topright", legend = c(1,2),pch = c(16,16), col = c("brown","black"), bty = "n",border = c("brown","black"),y.intersp = 0.7,cex = 2)
  } else if(format == "scPET"){
    cellNumber = length(orderBarcode)
    adjust_pos <-  0
    if(max(cellNumber) >= 150){
      point_ratio = 2
    }else if(max(cellNumber) >=50){
      point_ratio <- 150 / max(cellNumber) * 2 
    }else{
      point_ratio = 6
    }
    
    CC=c("#FF0000","#1E1EFF","#007A00","#800080","#FFBD46","#38D1D8","#89410E")

    plot(c(1,1),type='n',ylab="",xlim=c(start,end), ylim=c(0,cellNumber+1),xaxt="n",xaxs="i")
    if(length(orderBarcode)!=0){
      lapply(1:length(orderBarcode), function(i){
        tmpdata = result[result$barcode==orderBarcode[i],]
        if (nrow(tmpdata) != 0){
          lapply(1:nrow(tmpdata), function(row){
            x1 = (tmpdata$start1[row] + tmpdata$end1[row])/2
            x2 = (tmpdata$start2[row] + tmpdata$end2[row])/2
            height = 0.7
            vert_y  = height
            vert_x= (x1 + x2) / 2
            a = vert_y / ((vert_x - x1) * (vert_x - x2))
            if(direction){
              color_list <- c("#25CC69","#2F7E4F","#8080804A")
              color_list_curve <- c("#25CC69","#25CC69","#8080804A")
              ypos <- 1*i + 0.2 + adjust_pos
              curve(ypos + 0 * a * (x - x1) * (x - x2), from=x1, to=x2, add=TRUE, col=color_list_curve[tmpdata$anchor1[row]],lwd=0.25)
              draw.circle(x1,ypos,0.45,lwd = point_ratio,col=color_list[tmpdata$anchor1[row]],border=color_list[tmpdata$anchor1[row]])
              draw.circle(x2,ypos,0.45,lwd = point_ratio,col=color_list[tmpdata$anchor2[row]],border=color_list[tmpdata$anchor2[row]])
            }else{
              if(length(color_class)==0) {
                color = "#FF0000"
              }else{
                color = CC[color_class[i] %% length(CC) + 1]
              }
              ypos <- 1*i + 0.2 + adjust_pos
              curve(ypos + 0 * a * (x - x1) * (x - x2), from=x1, to=x2, add=TRUE, col=color,lwd=0.25)
              points(x1,ypos,lwd=point_ratio,cex=0.45,col=color,bg=color,pch=19)
              points(x2,ypos,lwd=point_ratio,cex=0.45,col=color,bg=color,pch=19)
            }
          })
          
        }
      })
    }
    if(showLabels)  mtext(paste0(tbl,"- ",cellNumber,"/",nrecords," cells; ",nrow(result)," PETs"), side = 3,cex=0.8,line=-1.5,adj=0)
  }else if(format == "hic"){
    file = result$filePath
    norm <- norm_vec
    res = res_vec
    if(is.na(res)) res = ress2[length(ress2)]
    chrom <- gsub("^chr","",chrom)

    file_hic <- straw(norm, file, paste0(chrom,":",start,":",end), paste0(chrom,":",start,":",end), "BP", res)

    file_hic$x2 <- (file_hic$x - start)/res + 1
    file_hic$y2 <- (file_hic$y - start)/res + 1
    file_hic$x3 <- sqrt(2)*file_hic$x2 + sqrt(2)/2*(file_hic$y2-file_hic$x2)
    file_hic$y3 <- sqrt(2)/2*(file_hic$y2-file_hic$x2)
    file_hic$modcounts <- ceiling(file_hic$counts)

    color = colorRampPalette(c(color1_vec,color2_vec),interpolate="linear")(max_vec)
    
    plot(c(min(file_hic$x3), max(file_hic$x3)), c(0,max(file_hic$y3)*hic_part), ylab="",type = "n",xaxt="n",xaxs="i",yaxs="i")
    lapply(1:nrow(file_hic), function(i){
      x = file_hic$x3[i]
      y = file_hic$y3[i]
      xvec=c(x-sqrt(2)/2,x,x+sqrt(2)/2,x)
      yvec=c(y,y+sqrt(2)/2,y,y-sqrt(2)/2)
      if(max_vec<file_hic$modcounts[i]){
        polygon(xvec,yvec,col=color[max_vec],border = color[max_vec]) 
      }else{
        polygon(xvec,yvec,col=color[file_hic$modcounts[i]],border = color[file_hic$modcounts[i]])   
      }
    })
    if(showLabels)  mtext(paste0(tbl,"- ","res: ",res," / ","norm: ",norm), side = 3,cex=0.8,line=-1.5,adj=0)
    legend("topright", legend = c(max_vec,"","","",0), fill = colorRampPalette(c(color1_vec,color2_vec),interpolate="linear")(5)[c(5,4,3,2,1)], bty = "n",border = colorRampPalette(c(color1_vec,color2_vec),interpolate="linear")(5)[c(5,4,3,2,1)],y.intersp = 0.2,cex = 1.8)
  }else if(format == "singleton"){
    cellNumber = length(orderBarcode)
    adjust_pos <-  0
    if(max(cellNumber) >= 150){
      point_ratio = 2
    }else if(max(cellNumber) >=50){
      point_ratio <- 150 / max(cellNumber) * 2 
    }else{
      point_ratio = 6
    }
    CC=c("#FF0000","#1E1EFF","#007A00","#800080","#FFBD46","#38D1D8","#89410E")
    ci = 0
    
    plot(c(1,1),type='n',ylab="",xlim=c(start,end), ylim=c(0,cellNumber+1),xaxt="n",xaxs="i")
    if(length(orderBarcode)!=0){
      lapply(1:length(orderBarcode), function(i){
        tmpdata = result[result$name==orderBarcode[i],]
        if (nrow(tmpdata) != 0){
          color = CC[i%%7+1]
          ci = ci + 1
          lapply(1:nrow(tmpdata), function(row){
            xleft = tmpdata$start[row]
            ybottom = 1*i - 0 + adjust_pos
            xright = tmpdata$end[row]
            ytop = 1*i + 0.05 + adjust_pos
            points(xleft,(ybottom+ytop)/2,lwd=point_ratio * 0.15 / 0.4,cex=0.45,col=color,bg=color,pch=19)
            points(xright,(ybottom+ytop)/2,lwd=point_ratio * 0.15 / 0.4,cex=0.45,col=color,bg=color,pch=19)
          })
          x1 = min(tmpdata$start)
          x2 = max(tmpdata$end)
          ypos <- 1*i + 0.025 + adjust_pos
          height = 0.7
          vert_y  = height
          vert_x= (x1 + x2) / 2
          a = vert_y / ((vert_x - x1) * (vert_x - x2))
          
          curve(ypos + 0 * a * (x - x1) * (x - x2), from=x1, to=x2, add=TRUE, col=color,lwd=0.25)
        }
      })
    }
    if(showLabels)  mtext(paste0(tbl,"-",cellNumber,"/",nrecords," cells"), side = 3,cex=0.8,line=-1.5,adj=0)

  }
  
  #print(proc.time())
  #print("finish")
}

sortedByNumber = function(result){
  order <- data.frame(barcode=names(table(result$name)),barcodecount=as.vector(table(result$name)))
  if(nrow(order)!=0){
    order <- order[order(order$barcodecount),]  
  }

  return(order$barcode)
}

sortedByDirection <- function(result,geneinfo,showAllPET=T){
  genestart = geneinfo[[1]]
  geneend = geneinfo[[2]]
  genestrand = geneinfo[[3]]
  if(genestrand == "-"){
    cut_pos = geneend
  }else{
    cut_pos = genestart
  }
  cut_start = cut_pos - 5000
  cut_end = cut_pos + 5000
  
  if(nrow(result)!=0){
    
    for(i in 1:nrow(result)){
      
      if((result$start1[i] > cut_start && result$end1[i] < cut_end) || (result$start2[i] > cut_start && result$end2[i] < cut_end)){
        dis1 <- (cut_start+cut_end) - (result$start1[i] + result$end1[i])
        dis2 <- (result$start2[i] + result$end2[i]) - (cut_start+cut_end)
        if(dis1 > dis2){
          result$cor[i] <- result$start1[i] + result$end1[i]
          result$anchor1[i] <- 2
          result$anchor2[i] <- 1
        }else {
          result$cor[i] <- result$start2[i] + result$end2[i]
          result$anchor1[i] <- 1
          result$anchor2[i] <- 2
        }
        if(result$start1[i] > cut_start && result$end1[i] < cut_end){
          result$anchor1[i] <- 1
        }
        if(result$start2[i] > cut_start && result$end2[i] < cut_end){
          result$anchor2[i] <- 1
        }    
      }else{
        result$cor[i] <- 0
        result$anchor1[i] <- 3
        result$anchor2[i] <- 3
      }
    }
    
    tmp <- result
    result <- result[(result$start1 > cut_start & result$end1 < cut_end)|(result$start2 > cut_start & result$end2 < cut_end),]
    result <- result[order(result$cor),]
    
    tmp1 <- result[(result$cor) < (cut_start+cut_end),][order(result$cor),]$barcode
    tmp1 <- unique(tmp1)
    tmp2 <- result[(result$cor) > (cut_start+cut_end),][order(-result$cor),]$barcode
    tmp2 <- unique(tmp2)
    tmp2 <- rev(tmp2)
    
    final_barcode <- c(na.omit(tmp1),na.omit(tmp2))
    out_order <- data.frame(orderBarcode=unique(final_barcode),color_class=rep(6,length(unique(final_barcode))))  
    
    if(showAllPET){
      return(list(out_order,tmp))
    }else{
      return(list(out_order,result))
    }
  }else{
    return(list(out_order<-data.frame(),result))
  }
}

clusterPET = function(result, start, end, nbin=50, nclass = 20,direction = F){
  begin = start
  last = end
  bin = (last - begin)/nbin
  binid = paste("bin", seq(1,nbin),sep="")
  tmpdat <- data.frame(matrix(0, nrow = length(unique(result$barcode)), ncol = length(binid), dimnames = list(unique(result$barcode), binid)))
  
  if(length(unique(result$barcode))>=2){
    if(length(unique(result$barcode)) <= nclass ){
      nclass = length(unique(result$barcode))
    }
    
    if(direction){
      for (i in (1:nrow(result))){
        s = floor(((result$start1[i] + result$end1[i])/2-begin)/bin) + 1
        e = floor(((result$start2[i] + result$end2[i])/2-begin)/bin) + 1
        if(s < 1){
          s = 1
        } 
        if(e > nbin){
          e = nbin
        } 
        tmpdat[result$barcode[i],paste("bin",s:e,sep="")] = 1
      }
    }else{
      for (i in (1:nrow(result))){
        s = floor(((result$start1[i] + result$end1[i])/2-begin)/bin) + 1
        e = floor(((result$start2[i] + result$end2[i])/2-begin)/bin) + 1
        if(s < 1){
          s = 1
        } 
        if(e > nbin){
          e = nbin
        } 
        tmpdat[result$barcode[i],paste("bin",e,sep="")] = 1
        tmpdat[result$barcode[i],paste("bin",s,sep="")] = 1
      }
    }
    
    
    PET_order <- data.frame(barcode=names(table(result$barcode)),barcodecount=as.vector(table(result$barcode)))
    PET_order <- PET_order[order(-PET_order$barcodecount),]
    PET_order_barcode <- PET_order$barcode
    
    tmpdat <- tmpdat[PET_order_barcode, ]
    # PH<-pheatmap(tmpdat, legend = FALSE, cluster_rows=TRUE, cluster_cols=FALSE, cutree_rows = nclass,
    #              border_color = NA,
    #              treeheight_row = 50,
    #              treeheight_col = 0,
    #              annotation_legend = FALSE,
    #              annotation_names_row = FALSE,
    #              show_rownames=FALSE,
    #              show_colnames = FALSE,
    #              silent = TRUE)
    # res_clu <- cutree(PH$tree_row, k = nclass)
    
    hc <- hclust(dist(tmpdat))
    res_clu <- cutree(hc, k = nclass)

    tmp = as.data.frame(res_clu)
    tmp$barcode = rownames(tmp)
    tmp = tmp[order(tmp$res_clu),]
    
    ### count each barcode length
    barcode_length_of_each_class <- data.frame(class = c(1:nclass),length = rep(0,nclass),num = rep(0,nclass))
    
    Total_PET = result
    
    ## mean
    # for(i in 1:(nrow(tmp))){
    #   barcode_length_of_each_class$length[tmp$res_clu[i]] = barcode_length_of_each_class$length[tmp$res_clu[i]] + (sum(Total_PET[Total_PET$barcode %in% tmp$barcode[i],]$end1)+sum(Total_PET[Total_PET$barcode %in% tmp$barcode[i],]$end2))/2 - (sum(Total_PET[Total_PET$barcode %in% tmp$barcode[i],]$start1)+sum(Total_PET[Total_PET$barcode %in% tmp$barcode[i],]$start2))/2
    #   barcode_length_of_each_class$num[tmp$res_clu[i]] = barcode_length_of_each_class$num[tmp$res_clu[i]] + 1
    # }
    # barcode_length_of_each_class$length_aver = barcode_length_of_each_class$length / barcode_length_of_each_class$num
    # barcode_length_of_each_class <- barcode_length_of_each_class[order(barcode_length_of_each_class$length_aver),]
    
    ## max
    for(i in 1:(nrow(tmp))){
      record = max(Total_PET[Total_PET$barcode %in% tmp$barcode[i],]$end2) - min(Total_PET[Total_PET$barcode %in% tmp$barcode[i],]$start1)
      if(barcode_length_of_each_class$length[tmp$res_clu[i]] < record) {
        barcode_length_of_each_class$length[tmp$res_clu[i]] = record
      }
    }
    barcode_length_of_each_class <- barcode_length_of_each_class[order(barcode_length_of_each_class$length),]
    
    ### smallest is the first
    barcode_length_of_each_class$length_order <- c(1:nclass)
    barcode_length_of_each_class <- barcode_length_of_each_class[order(barcode_length_of_each_class$class),]  
    for(i in 1:nrow(tmp)){
      tmp$length_order[i] <- barcode_length_of_each_class$length_order[tmp$res_clu[i]]
    }
    tmp <- tmp[order(tmp$length_order),]
    
    out_order <- data.frame(orderBarcode=tmp$barcode,color_class=tmp$length_order)
  }else{
    out_order <- data.frame(orderBarcode=result$barcode,color_class=c(1))
  }
  
  return(out_order)
}

clusterSingleton = function(result, start, end, nbin=50, nclass = 20){
  begin = start
  last = end
  bin = (last - begin)/nbin
  binid = paste("bin", seq(1,nbin),sep="")
  
  tmpdat <- data.frame(matrix(0, nrow = length(unique(result$name)), ncol = length(binid), dimnames = list(unique(result$name), binid)))

  if(length(unique(result$name))>=2){
    if(length(unique(result$name)) <= nclass ){
      nclass = length(unique(result$name))
    }
    
    orderBarcode <- unique(result$name)

    # lapply(1:length(orderBarcode), function(i){
    #   tmpdata = result[result$name==orderBarcode[i],]
    #   # lapply(1:nrow(tmpdata), function(row){
    #   #   s = floor(((tmpdata$start[row] + tmpdata$end[row])/2-begin)/bin) + 1
    #   #   if(s < 1){
    #   #     s = 1
    #   #   }else if(s > nbin){
    #   #     s = nbin
    #   #   }
    #   #   tmpdat[tmpdata$name[row],paste("bin",s,sep="")] = 1
    #   # })
    #   
    #   for(row in 1:nrow(tmpdata)){
    #     s = floor(((tmpdata$start[row] + tmpdata$end[row])/2-begin)/bin) + 1
    #     if(s < 1){
    #       s = 1
    #     }else if(s > nbin){
    #       s = nbin
    #     }
    #     tmpdat[tmpdata$name[row],paste("bin",s,sep="")] = 1
    #   }
    #   
    # })
    for(i in 1:length(orderBarcode)){
      tmpdata = result[result$name==orderBarcode[i],]
      for(row in 1:nrow(tmpdata)){
        s = floor(((tmpdata$start[row] + tmpdata$end[row])/2-begin)/bin) + 1
        if(s < 1){
          s = 1
        }else if(s > nbin){
          s = nbin
        }
        tmpdat[tmpdata$name[row],paste("bin",s,sep="")] = 1
      }
    }


  }
  
  write.csv(tmpdat,"tmpdat",quote=F)
  PET_order <- data.frame(barcode=names(table(result$name)),barcodecount=as.vector(table(result$name)))
  PET_order <- PET_order[order(-PET_order$barcodecount),]
  PET_order_barcode <- PET_order$barcode
  
  tmpdat <- tmpdat[PET_order_barcode, ]
  
  hc <- hclust(dist(tmpdat))
  res_clu <- cutree(hc, k = nclass)
  
  tmp = as.data.frame(res_clu)
  tmp$barcode = rownames(tmp)
  tmp = tmp[order(tmp$res_clu),]

  ### count each barcode length
  # barcode_length_of_each_class <- data.frame(class = c(1:nclass),length = rep(0,nclass),num = rep(0,nclass))
  # 
  # Total_PET = result
  # 
  # ## max
  # for(i in 1:(nrow(tmp))){
  #   record = max(Total_PET[Total_PET$name %in% tmp$barcode[i],]$end) - min(Total_PET[Total_PET$name %in% tmp$barcode[i],]$start)
  #   if(barcode_length_of_each_class$length[tmp$res_clu[i]] < record) {
  #     barcode_length_of_each_class$length[tmp$res_clu[i]] = record
  #   }
  # }
  # barcode_length_of_each_class <- barcode_length_of_each_class[order(barcode_length_of_each_class$length),]
  # 
  # ### smallest is the first
  # barcode_length_of_each_class$length_order <- c(1:nclass)
  # barcode_length_of_each_class <- barcode_length_of_each_class[order(barcode_length_of_each_class$class),]  
  # for(i in 1:nrow(tmp)){
  #   tmp$length_order[i] <- barcode_length_of_each_class$length_order[tmp$res_clu[i]]
  # }
  # tmp <- tmp[order(tmp$length_order),]
  # 
  # out_order <- data.frame(orderBarcode=tmp$barcode,color_class=tmp$length_order)
  out_order <- data.frame(orderBarcode=tmp$barcode,color_class=tmp$res_clu)
  
  return(out_order)
}

colorts <- function(color_str) {
  # 判断是否为颜色字符串
  # black
  color_names <- colors()
  if (tolower(color_str) %in% tolower(color_names)) {
    # library(gplots)
    # color_str = col2hex(color_str)
    return(color_str)
  }
  
  # 判断是否为 16 进制颜色代码
  # #234567
  hex_pattern <- "^#([A-Fa-f0-9]{6}|[A-Fa-f0-9]{3})$"
  if (grepl(hex_pattern, color_str)) {
    return(color_str)
  }
  # 234567
  hex_pattern <- "^([A-Fa-f0-9]{6}|[A-Fa-f0-9]{3})$"
  if (grepl(hex_pattern, color_str)) {
    return(paste0("#",color_str))
  }
  
  # 判断是否为 RGB 字符串
  # rgb(1,2,3) rgb(125,242,32)
  rgb_pattern <- "^rgb\\(\\s*([0-9]{1,3})\\s*,\\s*([0-9]{1,3})\\s*,\\s*([0-9]{1,3})\\s*\\)$"
  if (grepl(rgb_pattern, color_str, perl = TRUE)) {
    # color_nums <- as.numeric(gsub("rgb\\(|\\)", "", color_str, perl=TRUE, fixed=FALSE))
    color_nums <- as.numeric(strsplit(gsub("rgb\\(|\\)", "", color_str), ",")[[1]])
    return (rgb(color_nums[1], color_nums[2], color_nums[3], maxColorValue = 255))
  }
  # 125,242,32
  rgb_pattern <- "^\\s*([0-9]{1,3})\\s*,\\s*([0-9]{1,3})\\s*,\\s*([0-9]{1,3})\\s*$"
  if (grepl(rgb_pattern, color_str, perl = TRUE)) {
    color_nums <- as.numeric(strsplit(color_str, ",")[[1]])
    return (rgb(color_nums[1], color_nums[2], color_nums[3], maxColorValue = 255))
  }
}

getTrackInfo<- function(region="chr10:3600000-6600000",db="./data/test2/ChAIR.db",tbls=c(),adjustBinNumber=T){

  region = gsub("-|:"," ",region)
  regionarr = strsplit(trimws(region)," |\t")
  chrom = regionarr[[1]][1]
  start = as.numeric(regionarr[[1]][2])
  end = as.numeric(regionarr[[1]][3])
  
  con <- dbConnect(MySQL(), username = 'root', password = 'kkltmax8', dbname = db,host = 'localhost')
  min_vec <- c()
  max_vec <- c()
  format_vec <- c()
  color1_vec <- c()
  color2_vec <- c()
  tbl_vec <- c()
  norm_vec <- c()
  res_vec <- c()
  if(length(tbls) >=1){
    for(i in 1:length(tbls)){
      tbl <- tbls[i]
      
      
      format <- dbGetQuery(con, paste0("SELECT format FROM header WHERE sample = '",tbl,"'"))
      if(format=="loop"||format=="scPET"){
        result <- dbGetQuery(con,paste0("SELECT * FROM `",tbl,"` WHERE chromosome1 = '",chrom,"' AND end1 >= ",start," AND start1 <= ",end," AND ","chromosome2 = '",chrom,"' AND end2 >= ",start," AND start2 <= ",end))
        print(result)
      }else if(format=="hic"){
        result <- dbGetQuery(con,paste0("SELECT * FROM `",tbl,"`"))
      }else if(format=="bedgraph"||format=="RNAbedgraph"||format=="ATACbedgraph"){
        result <- dbGetQuery(con,paste0("SELECT * FROM `",tbl,"` WHERE chromosome = '",chrom,"' AND end >= ",start," AND start <= ",end))
        if(adjustBinNumber) result <- adjustBin(result,format,chrom,start,end)
      }else{
        result <- dbGetQuery(con,paste0("SELECT * FROM `",tbl,"` WHERE chromosome = '",chrom,"' AND end >= ",start," AND start <= ",end))
      }
      if(format=="bedgraph"||format=="RNAbedgraph"||format=="ATACbedgraph"){
        min = 0
        max = max(result$value)*1.05
        color1 = rgb(21/255,29/255,252/255)
        color2 = NA
        norm = NA
        res = NA
      }else if(format=="ssbedgraph"){
        min = min(result$value2)*1.05
        max = max(result$value1)*1.05
        color1 = "green"
        color2 = rgb(21/255,29/255,252/255)
        norm = NA
        res = NA
      }else if(format=="loop"){
        min = 0
        # max = log10(max(result$value)+1)*1.05
        max = max(result$value)*1.05
        color1 = "#FF0000"
        color2 = NA
        norm = NA
        res = NA
      }else if(format=="hic"){
        file <- result$filePath
        ress = readHicBpResolutions(file)
        ress2 = ress[order(-ress)]
        res = ress2[ress2<(end-start)/500][1]
        if(is.na(res)) res = ress2[length(ress2)]
        # chrom <- gsub("^chr","",chrom)
        # 
        # indexTab <- readHicChroms(file)
        # chrom <- indexTab[indexTab$name == chrom,"index"]
        
        norm = "NONE"
        file_hic <- straw(norm, file, paste0(gsub("^chr","",chrom),":",start,":",end), paste0(gsub("^chr","",chrom),":",start,":",end), "BP", res)
        
        
        min = 0
        max = quantile(file_hic$count,0.93)
        color1 = "#FFFFFF"
        color2 = "#FF0000"
        
      }else{
        min = NA
        max = NA
        color1 = NA
        color2 = NA
        norm = NA
        res = NA
      }
      
      tbl_vec <- c(tbl_vec,tbl)
      min_vec <- c(min_vec,min)
      max_vec <- c(max_vec,max)
      color1_vec <- c(color1_vec,color1)
      color2_vec <- c(color2_vec,color2)
      norm_vec <- c(norm_vec,norm)
      res_vec <- c(res_vec,res)
      ### R3.6 will treat "format" as data.frame so force it to transform into vector type
      format_vec <- c(format_vec,format[1,1]) 
    }  
  }
  
  dbDisconnect(con)
  
  dat <- data.frame(tbl_vec=tbl_vec,format_vec=format_vec,min_vec=min_vec,max_vec=max_vec,color1_vec=color1_vec,color2_vec=color2_vec,norm_vec=norm_vec,res_vec=res_vec)
  return(dat)
}

adjustBin <- function(result,format,chrom,start,end){
  
  step = as.integer((end-start)/2000)
  bin <- data.frame(chr=chrom,start=seq(start,end-step,step))
  bin$end <- bin$start+step
  bin$queryHits <- seq(1,nrow(bin),by=1)
  bin.gr <- makeGRangesFromDataFrame(bin)
  
  if(format=="bedgraph"||format=="RNAbedgraph"||format=="ATACbedgraph"){
    names(result) <- c("ID","chr","start","end","value","color","start_bin")
    result$subjectHits <- seq(1,nrow(result),by=1)
    result.gr <- makeGRangesFromDataFrame(result)
    x <- findOverlaps(bin.gr,result.gr)
    x.data <- as.data.frame(x)
    x.merge.data <- merge(x.data,result,by="subjectHits")
    x.final.data <- x.merge.data %>% group_by(queryHits) %>%  summarise(average = sum(value * (end-start))/step)
    x.bin.final.data <- merge(bin,x.final.data,by="queryHits")
    result_adjustBin <- x.bin.final.data[,c("chr","start","end","average")]
    names(result_adjustBin) <- c("chromosome","start","end","value")
  }else if(format=="ssbedgraph"){
    names(result) <- c("ID","chr","start","end","value1","value2","color1","color2","start_bin")
    result$subjectHits <- seq(1,nrow(result),by=1)
    result.gr <- makeGRangesFromDataFrame(result)
    x <- findOverlaps(bin.gr,result.gr)
    x.data <- as.data.frame(x)
    x.merge.data <- merge(x.data,result,by="subjectHits")
    x.final1.data <- x.merge.data %>% group_by(queryHits) %>%  summarise(average1 = sum(value1 * (end-start))/step)
    x.final2.data <- x.merge.data %>% group_by(queryHits) %>%  summarise(average2 = sum(value2 * (end-start))/step)
    x.bin.final.data <- merge(bin,merge(x.final1.data,x.final2.data,by="queryHits"),by="queryHits")
    result_adjustBin <- x.bin.final.data[,c("chr","start","end","average1","average2")]
    names(result_adjustBin) <- c("chromosome","start","end","value1","value2")
  }
  #print((result_adjustBin))
  #print(nrow(result))
  #print(nrow(result_adjustBin))
  return(result_adjustBin)
}

readData <- function(chrom,start,end,tbl,db,downsample,filterDis=8000){
  con <- dbConnect(MySQL(), username = 'root', password = 'kkltmax8', dbname = db,host = 'localhost')
  format <- dbGetQuery(con, paste0("SELECT format FROM header WHERE sample = '",tbl,"'"))
  
  if(format == "scPET"||format == "scATAC"||format == "scRNA"||format == "singleton"){
    nrecords <- dbGetQuery(con, paste0("SELECT nrecords FROM header WHERE sample = '",tbl,".barcodeindex'"))
  }else{
    nrecords <- 0
  }
  
  print(paste0(
    "SELECT * FROM `",tbl,"` WHERE chromosome1 = '",chrom,"' AND end1 >= ",start," AND start1 <= ",end," AND ","chromosome2 = '",chrom,"' AND end2 >= ",start," AND start2 <= ",end," AND barcode IN (",
    "SELECT tt.* FROM",
    "(SELECT barcode FROM `",tbl,".barcodeindex` LIMIT ",downsample,") as tt)"
  ))
  if(downsample > 0 && (format == "scPET"||format == "scATAC"||format == "scRNA")){
    nrecords <- downsample
    if(format == "scPET"){
      result <- dbGetQuery(con,paste0(
        "SELECT * FROM `",tbl,"` WHERE chromosome1 = '",chrom,"' AND end1 >= ",start," AND start1 <= ",end," AND ","chromosome2 = '",chrom,"' AND end2 >= ",start," AND start2 <= ",end," AND barcode IN (",
        "SELECT tt.* FROM",
        "(SELECT barcode FROM `",tbl,".barcodeindex` LIMIT ",downsample,") as tt)"
      ))
    }else{
      print(paste0(
        "SELECT * FROM `",tbl,"` WHERE chromosome = '",chrom,"' AND end >= ",start," AND start <= ",end," AND name IN (",
        "SELECT tt.* FROM",
        "(SELECT barcode FROM `",tbl,".barcodeindex` LIMIT ",downsample,") as tt)"
      ))
      result <- dbGetQuery(con,paste0(
        "SELECT * FROM `",tbl,"` WHERE chromosome = '",chrom,"' AND end >= ",start," AND start <= ",end," AND name IN (",
        "SELECT tt.* FROM",
        "(SELECT barcode FROM `",tbl,".barcodeindex` LIMIT ",downsample,") as tt)"
      ))
    }
  }else{
    if(format=="loop"||format=="scPET"){
      result <- dbGetQuery(con,paste0("SELECT * FROM `",tbl,"` WHERE chromosome1 = '",chrom,"' AND end1 >= ",start," AND start1 <= ",end," AND ","chromosome2 = '",chrom,"' AND end2 >= ",start," AND start2 <= ",end))
    }else if(format=="hic"){
      result <- dbGetQuery(con,paste0("SELECT * FROM `",tbl,"`"))
    }else{
      result <- dbGetQuery(con,paste0("SELECT * FROM `",tbl,"` WHERE chromosome = '",chrom,"' AND end >= ",start," AND start <= ",end))
    }    
  }
  dbDisconnect(con)
  
  return(list(result,format,nrecords))
}

readAnno <- function(chrom,start,end,anno,db){
  if(!is.null(anno)){
    con <- dbConnect(MySQL(), username = 'root', password = 'kkltmax8', dbname = db,host = 'localhost')
    #print(anno)
    resultanno <- dbGetQuery(con,paste0("SELECT * FROM `",anno,"` WHERE chromosome = '",chrom,"' AND end >= ",start," AND start <= ",end))
    dbDisconnect(con)
    
    numberofgeneinfo = length(names(table(resultanno$geneName)))
    namesgeneinfo = names(table(resultanno$geneName))
    starts = c()
    stops  = c()
    sizes  = c()
    strands = c()
    for (i in (1:numberofgeneinfo)) {
      subgeneinfo  = resultanno[which(resultanno$geneName == namesgeneinfo[i]),]
      starts = c(starts,min(subgeneinfo$start,subgeneinfo$end))
      stops  = c(stops, max(subgeneinfo$start,subgeneinfo$end))
      sizes  = c(sizes, stops[i] - starts[i])
      strands = c(strands,subgeneinfo$strand[1])
    }
    transcriptinfo = data.frame(names=namesgeneinfo,starts=starts,stops=stops,sizes=sizes,strand=strands)
    transcriptinfo = transcriptinfo[order(transcriptinfo$names),]
    rownames(transcriptinfo) <- transcriptinfo$names
    return(transcriptinfo)    
  }else{
    return(NULL)
  }

}

filterData <- function(result,format,mode,transcriptinfo,gene=""){

    genelist_loci <- c()
    genelist_name <- c()
    for(i in 1:nrow(transcriptinfo)){
      if(transcriptinfo$strand[i]=="+"){
        genelist_loci[i] <- transcriptinfo$starts[i]
      }else if(transcriptinfo$strand[i]=="-"){
        genelist_loci[i] <- transcriptinfo$stops[i]
      }
      genelist_name[i] <- transcriptinfo$names[i]
    }
    genelist <- data.frame(name=genelist_name,loci=genelist_loci)
    genelist$start <- genelist$loci - 5000
    genelist$end <- genelist$loci + 5000
    rownames(genelist) <- genelist$name
    if(mode=="TargetOneGene"){
      genelist <- genelist[gene,]  
    }

    record <- c()
    for(i in nrow(result):1){
      for(j in 1:nrow(genelist)){
        if((result$start1[i]>genelist$start[j]&&result$end1[i]<genelist$end[j])||(result$start2[i]>genelist$start[j]&&result$end2[i]<genelist$end[j])){
          record <- c(record,i)
        }
      }
    }
    selectBarcode <- result[record,]$barcode
    result <- result[result$barcode %in% selectBarcode,]
    return(result)


}

getFormatAndHeight <- function(chrom,start,end,db="hg38",tbl="MM10.RNA.COVER",hic_height=0){
  con <- dbConnect(MySQL(), username = 'root', password = 'kkltmax8', dbname = db,host = 'localhost')
  format <- dbGetQuery(con, paste0("SELECT format FROM header WHERE sample = '",tbl,"'"))
  dbDisconnect(con)
  if(format == "bed6"||format=="chromhmm"){
    height = "22px"
  }else if(format == "scPET"||format == "singleton"){
    height = "500px"
  }else if(format == "scRNA"||format == "scATAC"){
    height = "200px" 
  }else if(format == "anno"){
    #height = "80px"
    
    tmp_list <- readData(chrom=chrom,start=start,end=end,tbl=tbl,db=db,downsample=-1)    
    result <- tmp_list[[1]]

    #########################
    ### process plotrow  
    #########################
    numberofgeneinfo = length(names(table(result$geneName)))
    namesgeneinfo = names(table(result$geneName))
    starts = c()
    stops  = c()
    sizes  = c()
    strands = c()
    for (i in (1:numberofgeneinfo)) {
      subgeneinfo  = result[which(result$geneName == namesgeneinfo[i]),]
      starts = c(starts,min(subgeneinfo$start,subgeneinfo$end))
      stops  = c(stops, max(subgeneinfo$start,subgeneinfo$end))
      sizes  = c(sizes, stops[i] - starts[i])
      strands = c(strands,subgeneinfo$strand[1])
    }
    transcriptinfo = data.frame(names=namesgeneinfo,starts=starts,stops=stops,sizes=sizes,strand=strands)
    transcriptinfo = transcriptinfo[order(transcriptinfo$names),]
    typeshow = "squished"
    if (typeshow == "collapsed") {
      transcriptinfo$plotrow = 1
    } else if(typeshow == "squished") {
      if (nrow(transcriptinfo) > 0) {
        maxn = 1
        transcriptinfo$plotrow = 0
        for (i in (1:nrow(transcriptinfo))){
          overlapflag = T
          for (j in (1:maxn)) {
            transcriptinfotmp = transcriptinfo[transcriptinfo$plotrow==j &
                                                 transcriptinfo$starts[i] < transcriptinfo$stops & transcriptinfo$stops[i] > transcriptinfo$starts,]
            if (nrow(transcriptinfotmp) == 0) {
              transcriptinfo$plotrow[i] = j
              overlapflag = F
              break
            }
          }
          if (overlapflag) {
            maxn = maxn + 1
            transcriptinfo$plotrow[i] = maxn
          }
        }
      }
    } else if (typeshow == "flatten") {
      transcriptinfo$plotrow = seq(1:nrow(transcriptinfo))
    }
    height = paste0(max(transcriptinfo$plotrow)*40,"px")
  }else if(format == "hic"){
    height = paste0(hic_height,"px")

  }else{
    height = "80px"
  }
  
  return(height)
}

# run2(region="chr1:3600000-4900000",db="hg38",tbl="k562_S.k562_cellcycle_S.scPET",anno="hg38.simple.anno",mode="All",gene="AJAP1")
