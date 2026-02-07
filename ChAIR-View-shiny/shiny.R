pkgs=c("shiny","shinydashboard","shinyWidgets","shinyalert","sortable","scales","htmltools","dashboardthemes","plotly","DBI","RMySQL","RColorBrewer","plotrix","pheatmap","dplyr","GenomicRanges")
for(pkg in pkgs){
  if(require(pkg, character.only = TRUE)){
    print(paste("Package", pkg, "is installed and loaded correctly", sep = ""))
  } else{
    install.packages(pkg, repos="http://cran.rstudio.com/", dep = TRUE)
    if(require(pkg, character.only = TRUE)){
      print(paste("Package", pkg, "is installed and loaded correctly", sep = ""))
    } else{
      stop(paste("Couldn't install package", pkg, sep = " "));
    }
  }  
}

library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(shinyalert)
library(sortable)
library(scales)
library(htmltools)
library(dashboardthemes)
library(plotly)
library(patchwork)
library(shinyjs)
options(bitmapType='cairo') 　　

setwd("/home/rstudio/scChIATAC")
tmppath = "."

source("source2.r")

checkboxInput2 <- function (inputId, label, value = FALSE, width = NULL, height = NULL,font_size = NULL) {
  value <- restoreInput(id = inputId, default = value)
  inputTag <- tags$input(id = inputId, type = "checkbox", class = "shiny-input-checkbox")
  if (!is.null(value) && value) 
    inputTag$attribs$checked <- "checked"
  div(class = "form-group shiny-input-container", style = css(width = validateCssUnit(width),height = validateCssUnit(height), margin = paste0(c(0, 0, 0, 0), "px"), font_size = validateCssUnit(font_size)), 
      div(class = "checkbox", tags$label(inputTag, tags$span(label))))
}

checkboxInput3 <- function (inputId, label, value = FALSE, width = NULL, height = NULL,font_size = NULL) {
  value <- restoreInput(id = inputId, default = value)
  inputTag <- tags$input(id = inputId, type = "checkbox", class = "shiny-input-checkbox")
  if (!is.null(value) && value) 
    inputTag$attribs$checked <- "checked"
  div(class = "form-group shiny-input-container", style = css(width = validateCssUnit(width),height = validateCssUnit(height), margin = paste0(c(0, 0, 0, 0), "px"), font_size = validateCssUnit(font_size)), 
      div(class = "checkbox"))
}
#style="overflow-x: hiden;overflow-y:auto;"

plotOutput2 <- function (outputId, width = "100%", height = "400px", click = NULL, 
                         dblclick = NULL, hover = NULL, brush = NULL, inline = FALSE, 
                         fill = !inline){
  res <- imageOutput(outputId, width, height, click, dblclick, 
                     hover, brush, inline, fill)
  res$attribs$class <- "shiny-plot-output"
  res
}

con <- dbConnect(MySQL(), username = 'root', password = 'kkltmax8', dbname = "assembly",host = 'localhost')
assembly <- dbGetQuery(con,"SELECT * from assembly")
genomes <- assembly$assembly
dbDisconnect(con)

customTheme <- shinyDashboardThemeDIY(
  
  ### general
  appFontFamily = "Arial"
  ,appFontColor = "rgb(0,0,0)"
  ,primaryFontColor = "rgb(0,0,0)"
  ,infoFontColor = "rgb(0,0,0)"
  ,successFontColor = "rgb(0,0,0)"
  ,warningFontColor = "rgb(0,0,0)"
  ,dangerFontColor = "rgb(0,0,0)"
  ,bodyBackColor = "rgb(248,248,248)"
  
  ### header
  # ,logoBackColor = "rgb(23,103,124)"
  ,logoBackColor = "rgb(68,142,187)"
  
  ,headerButtonBackColor = "rgb(238,238,238)"
  ,headerButtonIconColor = "rgb(75,75,75)"
  ,headerButtonBackColorHover = "rgb(210,210,210)"
  ,headerButtonIconColorHover = "rgb(0,0,0)"
  
  ,headerBackColor = "rgb(238,238,238)"
  ,headerBoxShadowColor = "#aaaaaa"
  ,headerBoxShadowSize = "2px 2px 2px"
  
  ### sidebar
  ,sidebarBackColor = cssGradientThreeColors(
    direction = "down"
    # ,colorStart = "rgb(20,97,117)"
    # ,colorMiddle = "rgb(56,161,187)"
    # ,colorEnd = "rgb(3,22,56)"
    # ,colorStart = "rgb(68,142,187)"
    # ,colorMiddle = "rgb(68,142,187)"
    # ,colorEnd = "rgb(68,142,187)"
    ,colorStart = "rgb(255,255,255)"
    ,colorMiddle = "rgb(255,255,255)"
    ,colorEnd = "rgb(255,255,255)"
    ,colorStartPos = 0
    ,colorMiddlePos = 50
    ,colorEndPos = 100
  )
  ,sidebarPadding = 0
  
  ,sidebarMenuBackColor = "transparent"
  ,sidebarMenuPadding = 0
  ,sidebarMenuBorderRadius = 0
  
  ,sidebarShadowRadius = "3px 5px 5px"
  ,sidebarShadowColor = "#aaaaaa"
  
  ,sidebarUserTextColor = "rgb(255,255,255)"
  
  ,sidebarSearchBackColor = "rgb(55,72,80)"
  ,sidebarSearchIconColor = "rgb(153,153,153)"
  ,sidebarSearchBorderColor = "rgb(55,72,80)"
  
  ,sidebarTabTextColor = "rgb(0,0,0)"
  ,sidebarTabTextSize = 13
  ,sidebarTabBorderStyle = "none none solid none"
  ,sidebarTabBorderColor = "rgb(35,106,135)"
  ,sidebarTabBorderWidth = 1
  
  ,sidebarTabBackColorSelected = cssGradientThreeColors(
    direction = "right"
    # ,colorStart = "rgba(44,222,235,1)"
    # ,colorMiddle = "rgba(44,222,235,1)"
    # ,colorEnd = "rgba(0,255,213,1)"
    ,colorStart = "rgb(68,142,187)"
    ,colorMiddle = "rgb(68,142,187)"
    ,colorEnd = "rgb(68,142,187)"
    ,colorStartPos = 0
    ,colorMiddlePos = 30
    ,colorEndPos = 100
  )
  ,sidebarTabTextColorSelected = "rgb(0,0,0)"
  ,sidebarTabRadiusSelected = "0px 20px 20px 0px"
  
  ,sidebarTabBackColorHover = cssGradientThreeColors(
    direction = "right"
    # ,colorStart = "rgba(44,222,235,1)"
    # ,colorMiddle = "rgba(44,222,235,1)"
    # ,colorEnd = "rgba(0,255,213,1)"
    ,colorStart = "rgb(68,142,187)"
    ,colorMiddle = "rgb(68,142,187)"
    ,colorEnd = "rgb(68,142,187)"
    ,colorStartPos = 0
    ,colorMiddlePos = 30
    ,colorEndPos = 100
  )
  ,sidebarTabTextColorHover = "rgb(50,50,50)"
  ,sidebarTabBorderStyleHover = "none none solid none"
  ,sidebarTabBorderColorHover = "rgb(75,126,151)"
  ,sidebarTabBorderWidthHover = 1
  ,sidebarTabRadiusHover = "0px 20px 20px 0px"
  
  ### boxes
  ,boxBackColor = "rgb(255,255,255)"
  ,boxBorderRadius = 5
  ,boxShadowSize = "0px 1px 1px"
  ,boxShadowColor = "rgba(0,0,0,.1)"
  ,boxTitleSize = 16
  ,boxDefaultColor = "rgb(210,214,220)"
  ,boxPrimaryColor = "rgba(44,222,235,1)"
  ,boxInfoColor = "rgb(210,214,220)"
  ,boxSuccessColor = "rgba(0,255,213,1)"
  ,boxWarningColor = "rgb(244,156,104)"
  ,boxDangerColor = "rgb(255,88,55)"
  
  ,tabBoxTabColor = "rgb(255,255,255)"
  ,tabBoxTabTextSize = 14
  ,tabBoxTabTextColor = "rgb(0,0,0)"
  ,tabBoxTabTextColorSelected = "rgb(0,0,0)"
  ,tabBoxBackColor = "rgb(255,255,255)"
  ,tabBoxHighlightColor = "rgba(44,222,235,1)"
  ,tabBoxBorderRadius = 5
  
  ### inputs
  ,buttonBackColor = "rgb(245,245,245)"
  ,buttonTextColor = "rgb(0,0,0)"
  ,buttonBorderColor = "rgb(200,200,200)"
  ,buttonBorderRadius = 5
  
  ,buttonBackColorHover = "rgb(235,235,235)"
  ,buttonTextColorHover = "rgb(100,100,100)"
  ,buttonBorderColorHover = "rgb(200,200,200)"
  
  ,textboxBackColor = "rgb(255,255,255)"
  ,textboxBorderColor = "rgb(200,200,200)"
  ,textboxBorderRadius = 5
  ,textboxBackColorSelect = "rgb(245,245,245)"
  ,textboxBorderColorSelect = "rgb(200,200,200)"
  
  ### tables
  ,tableBackColor = "rgb(255,255,255)"
  ,tableBorderColor = "rgb(240,240,240)"
  ,tableBorderTopSize = 1
  ,tableBorderRowSize = 1
  
)

body <- dashboardBody(
  customTheme,
  fluidRow(
    column(3,
           dropdownButton(
             inputId = "mydropdown",
             label = "Global Control",
             status = "primary",
             circle = FALSE,

             uiOutput("anno_dropbt"),             
             radioButtons("mode", "Display Mode",choiceNames = list("All","Target Genes For scPET format","Target One Gene For scPET Format"),choiceValues = list("All","TargetGenes","TargetOneGene")),

             uiOutput("gene"),
             #radioButtons("shareBarcode","Barcode Condition",
             #             choices = list("Don't Share Barcode","Share Barcode And Sort By PET","Share Barcode And Sort By ATAC","Share Barcode And Sort By RNA"),),
             uiOutput("shareBarcode"),
             textInput("maxCellNumber","Max Cell Number",value="800"),
             textInput("downsample","Downsampled Cell Number",value="-1"),
             textInput("petCount","PET Count Threshold For Loop Format",value="2"),
             textInput("hic_part","Depth Ratio For hic Format",value="0.25"),
             checkboxInput("adjustBinNumber","Recast Bin Number Based On Mean For Bedgraph Format",value = TRUE),
             checkboxInput("showLabels","Show Labels",value = TRUE),
             actionButton("refresh","",icon("refresh"))
           ),
          offset = 0.2,
    ),
    column(3,
           dropdownButton(
             inputId = "mydropdown2",
             label = "Setting",
             icon = icon("gear"),
             status = "primary",
             circle = FALSE,
             uiOutput("setting")
           ),
    ),
    column(3,
           tags$div(
             style = "
                  font-size: 24px;           /* 字体大小 */
                  font-weight: bold;         /* 加粗 */
                  padding: 5px;              /* 内边距 */
                ",
                         textOutput("range")  # 要添加样式的文本输出
            )
    ),
    column(3,
           downloadButton("downloadData", "Download"))
     ),

  br(),

  fluidRow(
    
    tags$style(HTML(".my-div3{
      #display: grid;
      # margin-left: 1.25%;
      # margin-right: 0.16%;
      margin-left: 0%;
      margin-right: 0.15%;
      height: auto;
    }")),
    
    column(3,
          #selectInput("genome",label=NULL,choices = genomes,selected = useGenome)),
          uiOutput("genome")),
    # column(3,
    #        uiOutput("chrom")),
    # column(3,
    #        uiOutput('x')),
    # div(class="my-div",               
    #     checkboxInput3(paste0("Plot_check","0"),label=NULL,height="5px"),
    #     plotOutput("chromband",height="120px"),
    # ),
    # 
    column(3,
           textInput('region',label=NULL,value="chr1:4600000-4900000")),
    column(3,
           div(
              actionButton("left","<"),
              actionButton("right",">"),
              actionButton("wholeGenome","<->"),
              actionButton("zoomOut","-"),
              actionButton("zoomIn","+"),
              actionButton("submit","submit!",style = "background-color: #2196F3; color: white; border: none; padding: 6px 12px; border-radius: 4px;"),
              )
    ),
    column(3,
           selectInput("local",NULL,choices = c("Local Control","Autoscale","Group Autoscale"),selected = "Local Control"),
    ),
          ),
  # div(
  #   class="my-div3",
  # 
  #   plotlyOutput("chromband",height="90px"),
  # ),
  # 
  
  fluidRow(

    tags$style(HTML(".my-div{
      display: grid;
      # margin-left: 0.5%;
      # margin-right: 0%;
      # grid-template-columns: 1.5% 97.5%;
      margin-left: 0.8%;
      margin-right: 0.4%;
      grid-template-columns: 1% 98%;
      height: auto;
    }")), 
    tags$head(tags$script('
                                var width = 0
                                $(document).on("shiny:connected", function(e) {
                                    width = $("#containerUI").width();
                                    Shiny.onInputChange("width", width);
                                });
                                $(window).resize(function(e) {
                                    width = $("#containerUI").width();
                                    Shiny.onInputChange("width", width);
                                });
                            ')),
    
    div(class="my-div",               
        checkboxInput3(paste0("Plot_check","0"),label=NULL,height="5px"),
        plotOutput("chromband",height="120px"),
    ),
    
    br(),
    uiOutput("containerUI"),
    sortable_js(css_id = "containerUI",
                options = sortable::sortable_options(
      onSort = sortable::sortable_js_capture_input("container_selected")
    )),

    
  )
)

ui <- dashboardPage(
  dashboardHeader(title = "scMV"
  ),
  dashboardSidebar(
    radioButtons("group", "Group Mode",choiceNames = list("Individual","Group_vertical","Group_horizontal"),choiceValues = list("Individual","Group_vertical","Group_horizontal")),
    sidebarMenuOutput("upData"),
    collapsed = TRUE
  ),
  body
)
server <- function(input, output,session) {
  
  ######################################
  ######        initialize          #####
  ######################################
  record <- read.table(paste0(tmppath,"/record/record.csv"))
  rownames(record) <- record$id
  useGenome <- record["genome",]$stat
  
  genomeRecord <- reactiveVal(useGenome) 
  gene <- reactiveVal("")
  mode <- reactiveVal("All")
  shareBarcode <- reactiveVal("Don't Share Barcode")
  downsample <- reactiveVal(-1)
  #sameCell <- reactiveVal(FALSE)
  #nclass <- reactiveVal(20)
  showLabels <- reactiveVal(TRUE)
  
  con <- dbConnect(MySQL(), username = 'root', password = 'kkltmax8', dbname = useGenome,host = 'localhost')
  tmpdata <- dbGetQuery(con,"SELECT * FROM header WHERE format = 'anno'")
  dir <- tmpdata$sample
  dbDisconnect(con)
  anno_dropbt <- reactiveVal(dir[1])
  
  petCount <- reactiveVal(2)
  adjustBinNumber <- reactiveVal(TRUE)
  maxCellNumber <- reactiveVal(800)
  trackInfo <- reactiveVal("")
  trackInfo_origin <- reactiveVal("")
  
  file_list <- reactiveVal(c(""))
  containerCount_react <- reactiveVal(0)
  hic_part <- reactiveVal(0.25)

  ######################################
  #####       Load Genome          #####
  ######################################
  
  chrom <- reactiveVal("chr1")
  chrom_start <- reactiveVal(0)
  chrom_end <- reactiveVal(248956422)
  # observeEvent(input$chrom,{
  #   chrom(input$chrom)
  #   con <- dbConnect(MySQL(), username = 'root', password = 'kkltmax8', dbname = useGenome,host = 'localhost')
  #   sample <- dbGetQuery(con,"SELECT sample FROM header WHERE format = 'size'")
  #   result <- dbGetQuery(con,paste0("SELECT * from `",sample,"`"))
  #   dbDisconnect(con)
  #   chr_size <- result
  #   rownames(chr_size) <- chr_size$chromosome
  #   chrom_end(chr_size[chrom(),]$length)
  # })
  
  observeEvent(input$genome,{
    
    if(input$genome!=genomeRecord()){
      dat <- data.frame(id=c("genome"),stat=c(input$genome))
      write.table(dat,file=paste0(tmppath,"/record/record.csv"))
      session$reload()
    }
    
  })
  ######################
  ### Menu button    ###
  ######################
  
  observeEvent(input$group,{
    
    if(input$group == "Individual"){
      output$upData <- renderMenu({
        con <- dbConnect(MySQL(), username = 'root', password = 'kkltmax8', dbname = useGenome,host = 'localhost')
        tmpdata <- dbGetQuery(con,"SELECT * FROM header WHERE format <> 'size' AND format <> 'chrband' AND format <> 'scRNAindex' AND format <> 'scATACindex' AND format <> 'scPETindex' AND format <> 'jointindex'")
        dir <- unique(tmpdata$project)
        dbDisconnect(con)

        list <- list()
        if(length(dir>=1)){
          for(i in 1:length(dir)){
            files <- tmpdata[tmpdata$project %in% dir[i],]$sample
            list[[i]] <- menuItem(dir[i],icon = icon("dice-d6"),startExpanded = TRUE,
                                  checkboxGroupInput(dir[i],label = NULL,choices = files))
          }
        }
        sidebarMenu(list)
      })      
    }else if(input$group == "Group_vertical"||input$group == "Group_horizontal"){
      output$upData <- renderMenu({
        con <- dbConnect(MySQL(), username = 'root', password = 'kkltmax8', dbname = useGenome,host = 'localhost')
        tmpdata <- dbGetQuery(con,"SELECT * FROM header WHERE format <> 'size' AND format <> 'chrband' AND format <> 'scRNAindex' AND format <> 'scATACindex' AND format <> 'scPETindex' AND format <> 'jointindex' AND `group` <> ''")
        dir <- unique(tmpdata$project)
        dbDisconnect(con)

        list <- list()
        if(length(dir>=1)){
          for(i in 1:length(dir)){
            groups <- unique(tmpdata[tmpdata$project %in% dir[i],]$group)
            files <- c()
            for(j in 1:length(groups)){
              file <- tmpdata[tmpdata$group %in% groups[j],]$sample[1]
              format <- tmpdata[tmpdata$group %in% groups[j],]$format[1]
              file_name <- gsub(paste0(".",format),"",file)
              files <- c(files,file_name)
            }
            
            list[[i]] <- menuItem(dir[i],icon = icon("dice-d6"),startExpanded = TRUE,
                                  checkboxGroupInput(dir[i],label = NULL,choices = files))
          }
        }
        sidebarMenu(list)
      })
    }
  })
  
  ######################
  ### Region button ####
  ######################
  
  obs_record = c(4600000,4900000)
  var <- reactiveVal(c(4600000,4900000))
  var2 <- reactiveVal(c(4600000,4900000))
  var3 <- reactiveVal(c(4600000,4900000))

  # observeEvent(input$region,{
  #   if(!is.null(input$region)){
  #     region = gsub("-|:"," ",input$region)
  #     regionarr = strsplit(trimws(region)," |\t")
  #     chrom = regionarr[[1]][1]
  #     start = as.numeric(regionarr[[1]][2])
  #     end = as.numeric(regionarr[[1]][3])
  #     if(is.na(start)) {start=chrom_start()}
  #     if(is.na(end)) {end=chrom_end()}
  #     obs_record = c(start,end)
  #     var2(obs_record)
  #     var3(obs_record)
  #     chrom(chrom)
  #   }
  # })



  
  # 
  # observeEvent(input$region, {
  #   if (!is.null(input$region)) {
  #     # 关键：在 reactive 上下文中先捕获 input$region 的当前值
  #     current_region <- input$region  # 保存到临时变量
  #     
  #     # 延迟执行（使用捕获的 current_region，而非直接访问 input$region）
  #     later::later(
  #       func = function() {
  #         # 使用临时变量 current_region，而非 input$region
  #         region <- gsub("-|:", " ", current_region)
  #         regionarr <- strsplit(trimws(region), " |\t")
  #         chrom_val <- regionarr[[1]][1]  # 避免与 reactive 变量 chrom 重名
  #         start <- as.numeric(regionarr[[1]][2])
  #         end <- as.numeric(regionarr[[1]][3])
  #         
  #         if (is.na(start)) { start <- chrom_start() }
  #         if (is.na(end)) { end <- chrom_end() }
  #         
  #         obs_record <- c(start, end)
  #         var2(obs_record)
  #         var3(obs_record)
  #         chrom(chrom_val)  # 假设 chrom 是 reactiveVal
  #       },
  #       delay = 1  # 延迟1秒
  #     )
  #   }
  # })
  # 
  # observeEvent(event_data("plotly_relayout"),{
  #   d <- event_data("plotly_relayout")
  #   if(!is.null(d$`xaxis.range[0]`)){
  #     a <- c(as.integer(d$`xaxis.range[0]`),as.integer(d$`xaxis.range[1]`))
  #     if(a[1]<chrom_start()) {a[1]=chrom_start()}
  #     if(a[2]>chrom_end()) {a[2]=chrom_end()}
  #     var(a)
  #     var3(obs_record)
  #   }else if(!is.null(d$`xaxis.range`)){
  #     a <- c(as.integer(d$xaxis.range[1]),as.integer(d$xaxis.range[2]))
  #     #print(paste0(a[1],"<",chrom_start()))
  #     if(a[1]<chrom_start()) {a[1]=chrom_start()}
  #     if(a[2]>chrom_end()) {a[2]=chrom_end()}
  #     var(a)
  #     var3(obs_record)
  #   }
  # 
  # })
  
  ###########################
  ###   refresh button    ###
  ###########################
  
  observeEvent(input$refresh,{
    adjustBinNumber(input$adjustBinNumber)
    petCount(input$petCount)
    maxCellNumber(as.integer(input$maxCellNumber))
    gene(input$gene)
    mode(input$mode)
    downsample(as.integer(input$downsample))
    shareBarcode(input$shareBarcode)
    anno_dropbt(input$anno_dropbt)
    showLabels(input$showLabels)
    hic_part(as.numeric(input$hic_part))
  })
  
  ##########################
  ###   submit button    ###
  ##########################
  observeEvent(input$left,{
    step = (var2()[2] - var2()[1]) / 3
    start_ = as.integer(var2()[1] - step)
    end_ = as.integer(var2()[2] - step)
    if(start_ < chrom_start() ) start_ = chrom_start()
    if(end_ > chrom_end()) end_ = chrom_end()
    var2(c(start_,end_))
    updateTextInput(
      session, 
      "region", 
      value = paste0(chrom(), ":", start_,"-",end_)
    )
  })
  observeEvent(input$right,{
    step = (var2()[2] - var2()[1]) / 3
    start_ = as.integer(var2()[1] + step)
    end_ = as.integer(var2()[2] + step)
    if(start_ < chrom_start() ) start_ = chrom_start()
    if(end_ > chrom_end()) end_ = chrom_end()
    var2(c(start_,end_))
    updateTextInput(
      session, 
      "region", 
      value = paste0(chrom(), ":", start_,"-",end_)
    )
  })
  
  observeEvent(input$zoomIn,{
    step = (var2()[2] - var2()[1]) / 4
    start_ = as.integer(var2()[1] + step)
    end_ = as.integer(var2()[2] - step)
    if(start_ < chrom_start() ) start_ = chrom_start()
    if(end_ > chrom_end()) end_ = chrom_end()
    var2(c(start_,end_))
    updateTextInput(
      session, 
      "region", 
      value = paste0(chrom(), ":", start_,"-",end_)
    )
  })
  observeEvent(input$zoomOut,{
    step = (var2()[2] - var2()[1]) / 2
    start_ = as.integer(var2()[1] - step)
    end_ = as.integer(var2()[2] + step)
    if(start_ < chrom_start() ) start_ = chrom_start()
    if(end_ > chrom_end()) end_ = chrom_end()
    var2(c(start_,end_))
    updateTextInput(
      session, 
      "region", 
      value = paste0(chrom(), ":", start_,"-",end_)
    )
  })
  observeEvent(input$wholeGenome,{
    start_ = chrom_start()
    end_ = chrom_end()
    var2(c(start_,end_))
    updateTextInput(
      session, 
      "region", 
      value = paste0(chrom(), ":", start_,"-",end_)
    )
  })
  
  observeEvent(input$submit,{
    if(!is.null(input$region)){
      region = gsub("-|:"," ",input$region)
      regionarr = strsplit(trimws(region)," |\t")
      chrom = regionarr[[1]][1]
      start = as.numeric(regionarr[[1]][2])
      end = as.numeric(regionarr[[1]][3])
      if(is.na(start)) {start=chrom_start()}
      if(is.na(end)) {end=chrom_end()}
      obs_record = c(start,end)
      var2(obs_record)
      var3(obs_record)
      chrom(chrom)
    }
    genomeRecord(input$genome)

    ############### LOAD AND PLOT
    containerCount <- 0

    con <- dbConnect(MySQL(), username = 'root', password = 'kkltmax8', dbname = useGenome,host = 'localhost')
    tmpdata <- dbGetQuery(con,"SELECT * FROM header WHERE format <> 'size' AND format <> 'chrband'")
    dir <- unique(tmpdata$project)
    
    file_tmp_new <- c()
    ## dir is projects 
    for(i in 1:length(dir)){
      # dir2 is a project 
      dir2 = dir[i]
      if(!is.null(input[[dir2]])){
        # file_tmp is ✔ files of a project
        file_tmp <- input[[dir2]]
        
        if(input$group == "Individual"){
          file_tmp_new <- c(file_tmp_new,file_tmp)
        }else if(input$group == "Group_vertical"||input$group == "Group_horizontal"){
          for(file_tmp_part in file_tmp){
            group_name <- strsplit(sub("\\.","\t",file_tmp_part),"\t")[[1]][2]
            file_tmp_part_dat <- dbGetQuery(con,paste0("SELECT * FROM header WHERE `group` = '",group_name,"'"))
            format_order <- c("loop","ATACbedgraph","RNAbedgraph","scPET","scATAC","scRNA")
            file_tmp_part_dat_mod <- file_tmp_part_dat[file_tmp_part_dat$format %in% format_order,]
            rownames(file_tmp_part_dat_mod) <- file_tmp_part_dat_mod$format
            file_tmp_part_dat_mod_sample <- file_tmp_part_dat_mod[format_order,]$sample
            #file_tmp_new <- c(file_tmp_new,file_tmp_part_dat_mod_sample)
            file_tmp_new <- c(file_tmp_new,file_tmp_part_dat_mod_sample[1:3],anno_dropbt(),file_tmp_part_dat_mod_sample[4:6])
            
          }
        }
      }
    }
    
    dbDisconnect(con)
    file_list(file_tmp_new)
    containerCount = length(file_tmp_new)
    containerCount_react(containerCount)      
  

    ### record Yaxis
    
    observe({
      trackInfo <- getTrackInfo(region=paste0(chrom(),":",var2()[1],"-",var2()[2]),db=useGenome,tbls=file_list(),adjustBinNumber=adjustBinNumber())

      if(input$group == "Group_horizontal"||input$group == "Group_vertical"){
        lapply(1:6, function(i){
          for(j in 1:(containerCount/6)){
            if(exists("max_vec")){
              max_vec <- c(max_vec,trackInfo[i+(j-1)*6,]$max_vec)
            }else{
              max_vec <- trackInfo[i+(j-1)*6,]$max_vec
            }
            if(exists("min_vec")){
              min_vec <- c(max_vec,trackInfo[i+(j-1)*6,]$min_vec)
            }else{
              min_vec <- trackInfo[i+(j-1)*6,]$min_vec
            }
          }
          max <- max(max_vec)
          min <- min(min_vec)
          for(j in 1:(containerCount/6)){
            trackInfo[i+(j-1)*6,]$max_vec <- max
            trackInfo[i+(j-1)*6,]$min_vec <- min
          }
          rm(max_vec)
          rm(min_vec)
        })
        shareBarcode("Share Barcode And Sort By PET")
      }  

      trackInfo <- unique(trackInfo)
      rownames(trackInfo) <- trackInfo$tbl_vec
      trackInfo(trackInfo)
      trackInfo_origin(trackInfo)
    })
    
    ### output plots
    
    lapply(1:containerCount,function(i){
      output[[paste0("Plot",i)]] <- renderPlot({
        run2(region=paste0(chrom(),":",var2()[1],"-",var2()[2]),db=useGenome,tbl=file_list()[i],petCount=petCount(),adjustBinNumber=adjustBinNumber(),maxCellNumber=maxCellNumber(),mode=mode(),downsample=downsample(),gene=gene(),shareBarcode=shareBarcode(),anno=anno_dropbt(),showLabels=showLabels(),trackInfo=trackInfo()[i,],hic_part=hic_part())
      },res=72*1.5)
    })
    print("+++++++++++++++++++++++")
    print(input$group)
    if(input$group == "Individual"||input$group == "Group_vertical"){
      output$containerUI <- renderUI({
        if(containerCount >= 1){
          containers <- lapply(1:containerCount,function(i){
            div(class="my-div",`data-rank-id` = file_tmp_new[i],
                checkboxInput2(paste0("Plot_check",i),label=NULL,height="5px"),
                plotOutput2(paste0("Plot",i),height = getFormatAndHeight(chrom(),var2()[1],var2()[2],db=useGenome,tbl=file_tmp_new[i],hic_height=as.numeric(input$width)*0.5*hic_part()))
            )
          })
          tagList(containers)
        }

      })
    }else{
      output$containerUI <- renderUI({
        
        style <- tags$style(HTML(paste0(".my-div2{
          display: grid;
          grid-template-columns: repeat(",containerCount/7,",",98/(containerCount/7),"%);
          height: auto;
        }")))
        containers <- lapply(1:7,function(i){
          list <- lapply(1:(containerCount/7), function(j){
            num = (j-1)*7 + i
            div(class="my-div",               
                checkboxInput2(paste0("Plot_check",num),label=NULL,height="5px"),
                plotOutput2(paste0("Plot",num),height = getFormatAndHeight(chrom(),var2()[1],var2()[2],db=useGenome,tbl=file_tmp_new[num]))
            )
          })
          div(class="my-div2",`data-rank-id` = file_tmp_new[i],list)
        })

        tagList(style,containers)

      })      
    }
  })  
  
  ##########################
  ###   yaxis button    ###
  ##########################
  
  observeEvent(input$mydropdown2_state,{
    
    checks <- c()
    if(containerCount_react()>=1 && !is.null(input[[paste0("Plot_check","1")]])){
      for(i in 1:containerCount_react()){
          checks <- c(checks,input[[paste0("Plot_check",i)]])
      }
    }
    if(sum(checks)==1){
      output$setting  <- renderUI({
        if(trackInfo()[checks,]$format_vec=="hic"){
          con <- dbConnect(MySQL(), username = 'root', password = 'kkltmax8', dbname = useGenome,host = 'localhost')
          tbl <- trackInfo()[checks,]$tbl_vec
          result <- dbGetQuery(con,paste0("SELECT * FROM `",tbl,"`"))
          file <- result$filePath
          norms <- readHicNormTypes(file)
          ress = readHicBpResolutions(file)
          dbDisconnect(con)
          
          tagList(
            textInput("Yaxis.Min","Yaxis.Min",value=trackInfo()[checks,]$min_vec),
            textInput("Yaxis.Max","Yaxis.Max",value=trackInfo()[checks,]$max_vec),
            textInput("Yaxis.Color1","Yaxis.Color1",value=trackInfo()[checks,]$color1_vec),
            textInput("Yaxis.Color2","Yaxis.Color2",value=trackInfo()[checks,]$color2_vec),
            selectInput("Yaxis.Norm","Yaxis.Norm",choices = norms,selected = trackInfo()[checks,]$norm_vec),
            selectInput("Yaxis.Res","Yaxis.Res",choices = ress,selected = trackInfo()[checks,]$res_vec),
            actionButton("setting_refresh","",icon("refresh"))          
          )          
        }else{
          tagList(
            textInput("Yaxis.Min","Yaxis.Min",value=trackInfo()[checks,]$min_vec),
            textInput("Yaxis.Max","Yaxis.Max",value=trackInfo()[checks,]$max_vec),
            textInput("Yaxis.Color1","Yaxis.Color1",value=trackInfo()[checks,]$color1_vec),
            textInput("Yaxis.Color2","Yaxis.Color2",value=trackInfo()[checks,]$color2_vec),
            actionButton("setting_refresh","",icon("refresh"))          
          )
        }

      })
    }else{
      output$setting  <- renderText({
        "Please select only one checkbox!"
      })
    }
  })
  
  observeEvent(input$setting_refresh,{
    
    checks <- c()
    if(containerCount_react()>=1 && !is.null(input[[paste0("Plot_check","1")]])){
      for(i in 1:containerCount_react()){
        checks <- c(checks,input[[paste0("Plot_check",i)]])
      }
    }
    dat <- trackInfo()
    dat[checks,]$min_vec <- as.numeric(input$Yaxis.Min)
    dat[checks,]$max_vec <- as.numeric(input$Yaxis.Max)
    dat[checks,]$color1_vec <- input$Yaxis.Color1
    dat[checks,]$color2_vec <- input$Yaxis.Color2
    
    if(is.null(input$Yaxis.Norm)){
      dat[checks,]$norm_vec <- NA
    }else{
      dat[checks,]$norm_vec <- input$Yaxis.Norm  
    }
    
    if(is.null(input$Yaxis.Res)){
      dat[checks,]$res_vec <- NA
    }else{
      dat[checks,]$res_vec <- as.numeric(input$Yaxis.Res)
    }
    
    trackInfo(dat)
    
  })
  
  ##########################
  ###   local button    ###
  ##########################
  
  observeEvent((input$local),{
    if(input$local=="Group Autoscale"){
      vec <- c()
      if(containerCount_react()>=1 && !is.null(input[[paste0("Plot_check","1")]])){
        for(i in 1:containerCount_react()){
          if(input[[paste0("Plot_check",i)]]&&input$local=="Group Autoscale") {
            vec[i] <- TRUE
            updateSelectInput(session, "local",
                              label = NULL,
                              choices = c("Local Control","Autoscale","Group Autoscale"),
                              selected = "Local Control"
            )
          }else{
            vec[i] <- FALSE
          }
        }
        
        dat <- trackInfo()

        if(length(unique(dat[vec,]$format_vec))!=1){
          shinyalert(text = "Please select same format!", size="xs",imageWidth = "10px",imageHeight = "10px")
        }else{
          max <- max(dat[vec,]$max_vec)
          min <- min(dat[vec,]$min_vec)
          dat[vec,]$max_vec <- max
          dat[vec,]$min_vec <- min
          trackInfo(dat)
        }
      }
    }else if(input$local=="Autoscale"){
      vec <- c()
      if(containerCount_react()>=1 && !is.null(input[[paste0("Plot_check","1")]])){
        for(i in 1:containerCount_react()){
          if(input[[paste0("Plot_check",i)]]&&input$local=="Autoscale") {
            vec[i] <- TRUE
            updateSelectInput(session, "local",
                              label = NULL,
                              choices = c("Local Control","Autoscale","Group Autoscale"),
                              selected = "Local Control"
            )
          }else{
            vec[i] <- FALSE
          }
        }
        
        dat <- trackInfo()
        dat_origin <- trackInfo_origin()
        dat[vec,]$max_vec <- dat_origin[vec,]$max_vec
        dat[vec,]$min_vec <- dat_origin[vec,]$min_vec
        trackInfo(dat)
      }
    }
  })
  
  
  ##########################
  ######  Basic UI    ######
  ##########################
  output$range <- renderText({
    as.numeric(var2()[2]) - as.numeric(var2()[1])
  })   

  
  output$genome <- renderUI({
    selectInput("genome",label=NULL,choices = genomes,selected = useGenome)
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste0("figure.pdf")
    },
    content = function(file) {
      
      if(input$group == "Individual"||input$group == "Group_vertical"){
        
        heights <- lapply(1:length(file_list()), function(i){
          height <- getFormatAndHeight(chrom(),var2()[1],var2()[2],db=useGenome,tbl=file_list()[i])
          as.integer(gsub("px","",height))
        })
        hts <- unlist(heights)
        sum <- sum(hts)
        hts <- hts / sum * 0.96
        hts_1 <- c(1-0.02)
        hts_2 <- hts_1[1] - c(hts[1])
        if(length(hts)>=2){
          for(i in 2:length(hts)){
            hts_1[i] <- hts_2[i-1]
            hts_2[i] <- hts_1[i] - hts[i]
          }
        }
        plotnew <- c(F,rep(T,length(hts)-1))
        
        pdf(file,height=sum/480*7,width=7)
        if(is.null(input$container_selected)){
          lapply(1:length(file_list()), function(i){
            run2(region=paste0(chrom(),":",var2()[1],"-",var2()[2]),db=useGenome,tbl=file_list()[i],petCount=petCount(),adjustBinNumber=adjustBinNumber(),maxCellNumber=maxCellNumber(),mode=mode(),downsample=downsample(),gene=gene(),shareBarcode=shareBarcode(),anno=anno_dropbt(),showLabels=showLabels(),trackInfo=trackInfo()[i,],parfig = c(0,1,hts_2[i],hts_1[i]),plotnew = plotnew[i])
          })        
        }else{
          lapply(1:length(file_list()), function(i){
            run2(region=paste0(chrom(),":",var2()[1],"-",var2()[2]),db=useGenome,tbl=input$container_selected[i],petCount=petCount(),adjustBinNumber=adjustBinNumber(),maxCellNumber=maxCellNumber(),mode=mode(),downsample=downsample(),gene=gene(),shareBarcode=shareBarcode(),anno=anno_dropbt(),showLabels=showLabels(),trackInfo=trackInfo()[input$container_selected[i],],parfig = c(0,1,hts_2[i],hts_1[i]),plotnew = plotnew[i])
          })         
        }        
      }else{

        heights <- lapply(1:length(file_list()), function(i){
          height <- getFormatAndHeight(chrom(),var2()[1],var2()[2],db=useGenome,tbl=file_list()[i])
          as.integer(gsub("px","",height))
        })
        hts <- unlist(heights)
        plotnew <- c(F,rep(T,length(hts)-1))
        
        hts <- hts[1:6]
        sum <- sum(hts)
        hts <- hts / sum * 0.96
        hts_1 <- c(1-0.02)
        hts_2 <- hts_1[1] - c(hts[1])
        if(length(hts)>=2){
          for(i in 2:length(hts)){
            hts_1[i] <- hts_2[i-1]
            hts_2[i] <- hts_1[i] - hts[i]
          }
        }
        
        file_group_num <- length(file_list())/6
        
        pdf(file,height=sum/480*7,width=7*file_group_num)
        
        x1 <- seq(0,1,(1/file_group_num))[1:(file_group_num)]
        x2 <- seq(0,1,(1/file_group_num))[2:(file_group_num+1)]

        if(is.null(input$container_selected)){
          lapply(1:length(file_list()), function(i){
            run2(region=paste0(chrom(),":",var2()[1],"-",var2()[2]),db=useGenome,tbl=file_list()[i],petCount=petCount(),adjustBinNumber=adjustBinNumber(),maxCellNumber=maxCellNumber(),mode=mode(),downsample=downsample(),gene=gene(),shareBarcode=shareBarcode(),anno=anno_dropbt(),showLabels=showLabels(),trackInfo=trackInfo()[i,],parfig = c(x1[(i-1)%/%6+1],x2[(i-1)%/%6+1],hts_2[(i-1)%%6+1],hts_1[(i-1)%%6+1]),plotnew = plotnew[i])
          })
        }else{

          
          tmp <- factor(file_list()[1:6],levels = input$container_selected[2:7])
          order_tmp <- order(tmp)
          
          tmp_list <- as.character(tmp[order(tmp)])

          if(file_group_num>=2){
            for(i in 2:file_group_num){
              tmp <- file_list()[((i-1)*6+1):((i-1)*6+6)]
              tmp_list <- c(tmp_list,tmp[order_tmp])
            }            
          }
          lapply(1:length(file_list()), function(i){
            run2(region=paste0(chrom(),":",var2()[1],"-",var2()[2]),db=useGenome,tbl=tmp_list[i],petCount=petCount(),adjustBinNumber=adjustBinNumber(),maxCellNumber=maxCellNumber(),mode=mode(),downsample=downsample(),gene=gene(),shareBarcode=shareBarcode(),anno=anno_dropbt(),showLabels=showLabels(),trackInfo=trackInfo()[tmp_list[i],],parfig = c(x1[(i-1)%/%6+1],x2[(i-1)%/%6+1],hts_2[(i-1)%%6+1],hts_1[(i-1)%%6+1]),plotnew = plotnew[i])
          })
        }
      }
      dev.off()
    }
  )
  
  output$shareBarcode <- renderUI({
    
    if(input$group == "Group_horizontal"||input$group == "Group_vertical"){
      radioButtons("shareBarcode","Barcode Condition",
                   choices = list("Share Barcode And Sort By PET","Share Barcode And Sort By ATAC","Share Barcode And Sort By RNA"),selected = "Share Barcode And Sort By PET") 
    }else{
      radioButtons("shareBarcode","Barcode Condition",
                   choices = list("Don't Share Barcode","Share Barcode And Sort By PET","Share Barcode And Sort By ATAC","Share Barcode And Sort By RNA"),selected = "Don't Share Barcode")
    }
  })

  
  output$anno_dropbt <- renderUI({
    con <- dbConnect(MySQL(), username = 'root', password = 'kkltmax8', dbname = useGenome,host = 'localhost')
    tmpdata <- dbGetQuery(con,"SELECT * FROM header WHERE format = 'anno'")
    dir <- tmpdata$sample
    dbDisconnect(con)

    anno_dropbt(dir[1])
    selectInput("anno_dropbt","Annotation",dir)
    
  })
  
  output$gene <- renderUI({
    transcriptinfo <- readAnno(chrom(),var2()[1],var2()[2],anno=input$anno_dropbt,db=useGenome)
    genenames <- transcriptinfo$names
    
    
    if(input$mode=="All" || input$mode=="TargetGenes"){
      NULL
    }else{
      selectInput("gene","* Gene (For target one gene mode)",choices=genenames)  
    }
    
  })
    
  # output$chrom <- renderUI({
  # 
  #   con <- dbConnect(MySQL(), username = 'root', password = 'kkltmax8', dbname = useGenome,host = 'localhost')
  #   sample <- dbGetQuery(con,"SELECT sample FROM header WHERE format = 'size'")
  #   result <- dbGetQuery(con,paste0("SELECT * from `",sample,"`"))
  #   dbDisconnect(con)
  #   selectInput("chrom",label=NULL,choices = result$chromosome,selected = chrom())
  # 
  # })
  
  # output$x <- renderUI({
  #   textInput("region",NULL,value=paste0(chrom(),":",var()[1],"-",var()[2]))  
  # })
  
  output$chromband <- renderPlot({
    con <- dbConnect(MySQL(), username = 'root', password = 'kkltmax8', dbname = useGenome,host = 'localhost')
    sample <- dbGetQuery(con,"SELECT sample FROM header WHERE format = 'chrband'")
    result <- dbGetQuery(con,paste0("SELECT * from `",sample,"`"))
    chr_band_tb <- result
    dbDisconnect(con)

    chr_band_tb <- chr_band_tb[chr_band_tb$chromosome %in% chrom(),]
    width <- (chr_band_tb$end - chr_band_tb$start)/2
    color <- chr_band_tb$color
    names(color) <- color
    min <- min(chr_band_tb$start)
    max <- max(chr_band_tb$end)
    
    chrom_start(min)
    chrom_end(max)
    
    varstart <- var2()[1]
    varend <- var2()[2]

    m=12
    theme = theme(title=element_text(size=12,hjust=0.2,lineheight=0.2),
                  axis.title.x=element_text(size=m,hjust=0.5),
                  axis.title.y=element_blank(),
                  axis.text.x=element_text(size=m),
                  axis.text.y=element_blank(),
                  axis.ticks.y = element_blank(),
                  axis.ticks.x = element_line(size=1),
                  #axis.line = element_line(colour = "black"),
                  panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.background = element_blank(),
                  panel.border = element_rect(color = "black", linewidth = 1, fill = NA),legend.position="none",
                  legend.text = element_blank(),
                  legend.title = element_blank(),
                  plot.title = element_blank())
 
    p <- ggplot(chr_band_tb, aes(x=start+(end-start)/2, y=1, group=color,fill=color)) +
      geom_bar(stat="identity", width=width) +
      scale_fill_manual(values = color) + theme + scale_x_continuous(limits=c(min,max),labels = comma,expand=c(0,0)) + xlab("")
    p1 <- ggplot(chr_band_tb, aes(x=start+(end-start)/2, y=1, group=color,fill=color)) +
      geom_bar(stat="identity", width=width) +
      scale_fill_manual(values = color) + theme + scale_x_continuous(limits=c(var2()[1],var2()[2]),labels = comma,expand=c(0,0)) + xlab("")
    # p_sty <- style(p, hoverinfo="none")
    # ggplotly(p_sty, dynamicTicks = TRUE) %>% 
    #   rangeslider(
    #     start=var2()[1],
    #     end=var2()[2],
    #     thickness=0.70,
    #     borderwidth=2) %>%
    #   layout(
    #   yaxis = list(
    #     title = NA,
    #     fixedrange = T,
    #     autorange = F,
    #     tickmode = 'auto',
    #     nticks = 10),
    #   xaxis = list(
    #     title = NA,
    #     ticktext = as.list(format(seq(var2()[1]+(var2()[2]-var2()[1])/8,var2()[2]-(var2()[2]-var2()[1])/8,(var2()[2]-var2()[1])/8),big.mark = ",")), 
    #     tickvals = as.list(seq(var2()[1]+(var2()[2]-var2()[1])/8,var2()[2]-(var2()[2]-var2()[1])/8,(var2()[2]-var2()[1])/8)),
    #     tickmode = "array"
    #   ),
    #   margin = list(
    #     l = 24,
    #     b = 40
    #   ))
    p / p1
      
  })
}
#shinyApp(ui, server)
shinyApp(ui, server,options = list(launch.browser = TRUE))
