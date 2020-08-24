#circleapp for shiny

library(shiny)
library(circlize)
library(ggplot2)
library(gridExtra)

celltype<-read.csv(file ="celltypes.csv", stringsAsFactors = FALSE)[,c(1,3,4)]
celltype$Cell.type[22]<-"Cdk1+"
celltype$Cell.type[14]<-"NK_Tcell_Bcell"
celltype$Cell.type[c(12,13,19)]<-"Extras"
map<-celltype[!duplicated(celltype[,2:3]),]
map<-map[order(map[,3]),]
map<-map[order(map$Cell.type),]
map<-map[which(!(map$Cell.type %in% "Extras")),]

celltypes<-as.list(map$Cell.type)


linkdata<-read.csv(file="linkData-All.csv")


timepoints<-list("LPS0hr","LPS1hr","LPS4hr","LPS16hr","LPS27hr","LPS36hr","LPS48hr")

#if(interactive()){
ui <- fluidPage( 
  titlePanel(h1("Ligand-Receptor Interactions",align = "center")),
  sidebarLayout(
  sidebarPanel(
    h4("Sepsis TimePoints",align="center"),
               selectInput("times",
                           label = "Timepoints",
                           choices = timepoints,
                           selected = timepoints[1]
                           ),
    selectInput("mainCell",
                        label = "Focused Cell Type",
                        choices =  append("All",celltypes),
                        selected = "All"
                        ),
               
    
    h4("Cell Types",align = "center"),
               p("Must select EXACTLY 8 cell types",align="center"),
            checkboxGroupInput("cells",
                        label = "Cell Types 1-8:",
                        choices = celltypes
                        ,selected = celltypes[c(13,15,12,8,6,5,3,2)]
                        ),
    #selectInput("cellPoint",label = "Cell Type",choices = append("Select One",celltypes),selected = "Select One"),
    #actionButton("generateButton","Generate Graph")
    submitButton("Generate Graph"),
    p(" ",align="center"),
    p(" ",align="center"),
    p("How to read graph:",align="center"),
    p(" ",align="center"),
     img(src='for_circle.jpg',align="center",height='200px',width='200px')
               ,width = 2),
  
  mainPanel(h2("Cell Interactions",align = "center"),
          
            textOutput("selected_var2"),
            fluidRow(
              style= "border: 15px dashed white",splitLayout(
                cellWidths = c("75%","25%"),plotOutput("circle",click = "plot_click",height = '60%',width = "60%"),
            fluidRow(style="overflow-y: scroll; height:600px; width:250px",tableOutput("table"))
              )
            )
            
            #verbatimTextOutput("info")
            #textOutput("info"),
            #tableOutput("pointsTable")
            
            )
  

  ),
  
  fluidRow(
    column(width = 1, align = "center", img(src="https://brand.iu.edu/images/resources-trident.jpg", height=50, width=50)),
    column(width = 11,
           p( "Jered Myslinski",
              br(),
              "Takashi Hato",
              a("thato@iu.edu", href="mailto:thato@iu.edu")
           )
    )
  )
  
  
)

server <- function(input, output, session) {
  
  outVar = reactive({df = input$cells; df})
#  observe({ updateSelectInput(session, "mainCell",choices = as.list(outVar),selected = "S1") })
  observe({ 
    dd = input$cells
    updateSelectInput(session, "mainCell",choices = append("All",as.list(dd)),selected = "All")
    updateSelectInput(session, "cellPoint",choices = append("Select One",as.list(dd)), selected = "Select One")
  })
  
    output$selected_var2<-renderText({
      input$cells
  
  
  
  })
  
    output$table<-renderTable({
      time<-which(timepoints %in% input$times)
      
      #groups8<-paste(input$cells,collapse = ",")
      groups8<-map$newnum[which(map$Cell.type %in% unlist(input$cells))]
      
      cellmapping<-data.frame(groups8)
      cellmapping$types<-map[which(map[,3] %in% groups8),2]
      interacts<-NULL
      for(q in 1:7){
        tmp<-linkdata[which(linkdata$timepoint==q),c(1:7)]
        tmp<-as.character(unique(tmp[which((tmp[,1] %in% c(groups8))&(tmp[,4] %in% c(groups8))),3]))
        interacts<-append(interacts,tmp)
      }
      interacts<-unique(interacts)
      interactionTable<-data.frame(interacts)
      interactionTable$number<-c(1:dim(interactionTable)[1])
      interactionTable<-interactionTable[,c(2,1)]
     interactionTable
    })
    
   output$circle <-renderPlot({
     #input$generateButton
     time<-which(timepoints %in% input$times)
     
     #groups8<-paste(input$cells,collapse = ",")
     groups8<-map$newnum[which(map$Cell.type %in% unlist(input$cells))]
    
    cellmapping<-data.frame(groups8)
    cellmapping$types<-map[which(map[,3] %in% groups8),2]
     
    
    interacts<-NULL
    for(q in 1:7){
      tmp<-linkdata[which(linkdata$timepoint==q),c(1:7)]
      tmp<-as.character(unique(tmp[which((tmp[,1] %in% c(groups8))&(tmp[,4] %in% c(groups8))),3]))
      interacts<-append(interacts,tmp)
    }
    interacts<-unique(interacts)
    interlens<-round(length(interacts)/10)
    
    
      linkData<-linkdata[linkdata$timepoint==time,c(1:8)]
      templink<-linkData[which((linkData[,1] %in% c(groups8))&(linkData[,4] %in% c(groups8))),]
     
      tempcircle<-templink
      colnames(tempcircle)[c(4:6)]<-colnames(tempcircle)[c(1:3)]
      tempcircle<-rbind(tempcircle[,c(1,3,7,8)],tempcircle[,c(4,6,7,8)])
      colnames(tempcircle)[2]<-"interactPositionA"
      tt<-tempcircle$groupA
      po<-as.character(tempcircle$interactPositionA)
      for(tr in 1:length(tt)){
        tt[tr]<-cellmapping[which(cellmapping[,1] %in% tt[tr]),2]
        po[tr]<-which(interacts %in% po[tr])
      }
      tempcircle$groupA<-tt
      tempcircle$interactPositionA<-as.numeric(po)
      #tempcircle<-tempcircle[which(tempcircle[,1] %in% 0:7),]
      #tempcircle$groupA<-paste0("Group ",tempcircle$groupA)
      tempcircle$strength<-as.numeric(tempcircle$strength)
      tempcircle<-data.frame(factors = tempcircle$groupA,
                             x = tempcircle$interactPositionA, y = tempcircle$strength, z=tempcircle$opacity)
      tempcircle<-tempcircle[order(tempcircle[,1]),]
      tem<-as.character(unique(tempcircle[,1]))
      for(u in 1:8){
        temm<-c(tem[u],1,-0.5,0.001)
        temm<-rbind(temm,c(tem[u],length(interacts),-0.5,0.001))
        temm<-as.data.frame(temm)
        colnames(temm)<-colnames(tempcircle)
        tempcircle<-rbind(tempcircle,temm)
      }
      tempcircle[,2]<-as.numeric(tempcircle[,2])
      tempcircle[,3]<-as.numeric(tempcircle[,3])
      tempcircle[,4]<-as.numeric(tempcircle[,4])
      
      grps<-groups8
      set.seed(999)
      n = dim(tempcircle)[1]
      df = data.frame(factors = tempcircle$factors,
                      x = rnorm(n), y = runif(n))
      df$x<-tempcircle$x
      df$y<-tempcircle$y
      
      
      if(input$mainCell=="All"){
      colors<-rainbow(8)
      circos.par("track.height" = 0.1)
      circos.initialize(factors = df$factors, x = df$x)
      circos.track(factors = df$factors, y = df$y,
                   panel.fun = function(x, y) {
                     circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + uy(7, "mm"), 
                                 CELL_META$sector.index)
                     circos.axis(labels.cex = 0.7,h = "bottom",direction = "inside",
                                 major.at = c(seq(1,(interlens*10)+1,10)),labels.facing = "clockwise")
                     circos.yaxis(side = "left",labels.cex = 0.01)
                     #circos.dendrogram(dend = df$x)
                     #get.cell.meta.data("sector.index")
                   })
      col = colors[c(1:8)]
      circos.trackPoints(df$factors, df$x, df$y, col = col, pch = 19, cex = 0.7)
      title(main=paste0(timepoints[time]))
      for(r in 1:length(templink$groupA)){
        
        circos.link(cellmapping[which(cellmapping[,1] %in% templink[r,1]),2],which(interacts %in% templink[r,3]), 
                    cellmapping[which(cellmapping[,1] %in% templink[r,4]),2],which(interacts %in% templink[r,3]),
                    col = alpha(colors[which(cellmapping$groups8 %in% templink[r,1])],
                                templink[r,8]), arr.lwd = 1)
      }
  
     
      }else{
        
        
        
          grps<-groups8
          v<- cellmapping[which(cellmapping$types %in% input$mainCell),1]
          set.seed(999)
          n = dim(tempcircle)[1]
          df = data.frame(factors = tempcircle$factors,
                          x = rnorm(n), y = runif(n))
          df$x<-tempcircle$x
          df$y<-tempcircle$y
          #df$factors<-tempcircle$factors
          
        
          colors<-rainbow(8)
          circos.par("track.height" = 0.1)
          circos.initialize(factors = df$factors, x = df$x)
          circos.track(factors = df$factors, y = df$y,
                       panel.fun = function(x, y) {
                         circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + uy(7, "mm"), 
                                     CELL_META$sector.index)
                         circos.axis(labels.cex = 0.7,h = "bottom",direction = "inside",
                                     major.at = c(seq(1,(interlens*10)+1,10)),labels.facing = "clockwise")
                         circos.yaxis(side = "left",labels.cex = 0.01)
                         #circos.dendrogram(dend = df$x)
                         #get.cell.meta.data("sector.index")
                       })
          col = colors[c(1:8)]
          circos.trackPoints(df$factors, df$x, df$y, col = col, pch = 19, cex = 0.7)#,track.index = interacts)
          title(main=timepoints[time])
          for(r in 1:length(templink$groupA)){
            if(templink[r,1]==v){
              circos.link(cellmapping[which(cellmapping[,1] %in% templink[r,1]),2],which(interacts %in% templink[r,3]), 
                          cellmapping[which(cellmapping[,1] %in% templink[r,4]),2],which(interacts %in% templink[r,3]),
                          col = alpha(colors[which(cellmapping$groups8 %in% templink[r,1])],
                                      templink[r,8]), arr.lwd = 1 )
            }else{
              circos.link(cellmapping[which(cellmapping[,1] %in% templink[r,1]),2],which(interacts %in% templink[r,3]), 
                          cellmapping[which(cellmapping[,1] %in% templink[r,4]),2],which(interacts %in% templink[r,3]),
                          #col = alpha(colors[which(levels(df$factors) %in% cellmapping[templink[r,1],2])],0.02)
                          col = alpha("grey",0.05)
                          , arr.lwd = 1 )
            }
            
          }
    
        
        
        
        }
      
      
      
      },height = 800, width = 800)
 testing<-FALSE
 if(testing==TRUE){
  output$info <-renderText({
    time<-which(timepoints %in% input$times)
    
    #groups8<-paste(input$cells,collapse = ",")
    groups8<-map$newnum[which(map$Cell.type %in% unlist(input$cells))]
    
    cellmapping<-data.frame(groups8)
    cellmapping$types<-map[which(map[,3] %in% groups8),2]
    
    
    interacts<-NULL
    for(q in 1:7){
      tmp<-linkdata[which(linkdata$timepoint==q),c(1:7)]
      tmp<-as.character(unique(tmp[which((tmp[,1] %in% c(groups8))&(tmp[,4] %in% c(groups8))),3]))
      interacts<-append(interacts,tmp)
    }
    interacts<-unique(interacts)
    interlens<-round(length(interacts)/10)
    
    
    linkData<-linkdata[linkdata$timepoint==time,c(1:8)]
    templink<-linkData[which((linkData[,1] %in% c(groups8))&(linkData[,4] %in% c(groups8))),]
    
    tempcircle<-templink
    colnames(tempcircle)[c(4:6)]<-colnames(tempcircle)[c(1:3)]
    tempcircle<-rbind(tempcircle[,c(1,3,7,8)],tempcircle[,c(4,6,7,8)])
    colnames(tempcircle)[2]<-"interactPositionA"
    tt<-tempcircle$groupA
    po<-as.character(tempcircle$interactPositionA)
    for(tr in 1:length(tt)){
      tt[tr]<-cellmapping[which(cellmapping[,1] %in% tt[tr]),2]
      po[tr]<-which(interacts %in% po[tr])
    }
    tempcircle$groupA<-tt
    tempcircle$interactPositionA<-as.numeric(po)
    #tempcircle<-tempcircle[which(tempcircle[,1] %in% 0:7),]
    #tempcircle$groupA<-paste0("Group ",tempcircle$groupA)
    tempcircle$strength<-as.numeric(tempcircle$strength)
    tempcircle<-data.frame(factors = tempcircle$groupA,
                           x = tempcircle$interactPositionA, y = tempcircle$strength, z=tempcircle$opacity)
    tempcircle<-tempcircle[order(tempcircle[,1]),]
    tem<-as.character(unique(tempcircle[,1]))
    for(u in 1:8){
      temm<-c(tem[u],1,-0.5,0.001)
      temm<-rbind(temm,c(tem[u],length(interacts),-0.5,0.001))
      temm<-as.data.frame(temm)
      colnames(temm)<-colnames(tempcircle)
      tempcircle<-rbind(tempcircle,temm)
    }
    tempcircle[,2]<-as.numeric(tempcircle[,2])
    tempcircle[,3]<-as.numeric(tempcircle[,3])
    tempcircle[,4]<-as.numeric(tempcircle[,4])
    
    grps<-groups8
    set.seed(999)
    n = dim(tempcircle)[1]
    df = data.frame(factors = tempcircle$factors,
                    x = rnorm(n), y = runif(n))
    df$x<-tempcircle$x
    df$y<-tempcircle$y
    #nearPoints(df[,2:3], input$plot_click,threshold = 10, maxpoints = 1, addDist = TRUE)
    
    print(input$plot_click)
  })
 
     
  output$pointsTable<-renderTable({ 
    if(input$pointCell=="Select One"){}else{
      time<-which(timepoints %in% input$times)
      groups8<-map$newnum[which(map$Cell.type %in% unlist(input$cells))]
      groups8
   #   print(time);print(input$cellPoint)
  #tab<-linkdata[which((linkdata$timepoint==time)&(linkdata$interactPositionA>=1)&(linkdata$interactPositionA<=171)&(linkdata$groupA %in% map$newnum[which(map$Cell.type %in% input$cellPoint)])&(linkdata$groupB %in% c(13,15,12,8,6,5,3,2))),c(1:4,7)]
    }
    
  })
}
  
  } 



#interactive
#}
# Create Shiny app ----
shinyApp(ui = ui, server = server)

