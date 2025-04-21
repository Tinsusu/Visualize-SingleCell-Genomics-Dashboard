
source('global.R')
ui <- dashboardPage(
  dashboardHeader(title = "scRNAseq Plot"),
  dashboardSidebar(
      sidebarMenu(id='tab',
                ## add useShinyjs, this one for anything in UI can be disable or enable dynamicly
                useShinyjs(),
                ## The first name is for human, the second name is the name for internal operation
                ## So the server is going to recognize home while human will see Kick thing off here!
                ## the third component is icon, here we use dna icon
                menuItem("Kick Things Off Here!", tabName = "home", icon = icon("dna")),
                menuItem("scRNAseq Investigator", tabName = "input", icon = icon("dot-circle")),
                ## This part will appear on the side bar when a certain condition is met.
                ## We only need it when we click scRNAseq icon >>> input.tab is kinda in-built input
                conditionalPanel(condition = "input.tab == 'input'",
                      div(
                        ### fileInput is also a buil-in function "file" is id for system but user goingto see Please upload
                          fileInput("file", "Please Upload File", multiple=TRUE, accept=c('.rds')),
                          actionButton("reset", "Reset", icon = icon("undo"), 
                                       style = "color: #fff; background-color: #C69C6D; border-radius: 50%; width: 87.25%; padding: 15px; font-size: 18px; text-align: center;"),
                          actionButton("run", "Run", icon = icon("play"), 
                                       style = "color: #fff; background-color: #77DD77; width: 87.25%; padding: 15px; font-size: 18px; text-align: center;")
                        )
            )
    )
  ), 
  
  ## large space that will eventually contain results for us
  ## We'r gonna set two menu items, which depend on what item we click on
  dashboardBody(
    tabItems(
      tabItem(tabName = "input", # tabItem refers to tab in sidebar (not main panel)
              #tabsetpanel for we can add some subtabs too.
              tabsetPanel(id = 'main_tabs',
                          tabPanel("Follow Me please",
                                   includeMarkdown("./markdown/instructions.md")
                          )
              )
      ),
      tabItem(tabName = "home", #when we select home from slide bar
              tags$h1(HTML("
              <div style='
                font-family: \"Comic Sans MS\", cursive, sans-serif;
                font-size: 36px;
                 text-align: center;
                background: linear-gradient(90deg, red, orange, yellow, green, blue, indigo, violet);
                -webkit-background-clip: text;
                -webkit-text-fill-color: transparent;
                font-weight: bold;
                padding: 20px;
                border-bottom: 5px solid black;
                display: inline-block;
              '>
                🦄 Welcome to the Single-cell RNAseq Seurat Analysis App 🌈
              </div>
            ")) ,
              
      )
    )
  )         
)

### input refers to anything that we interact with the app, click run, uploadfile ....
### output refers to anything display on our app
server <- function(input, output, session) {
     # we need to extend the size of large file because shiny has a limit size
     options(shiny.maxRequestSize = 300*1024^2) #shiny has a limit size, so do this to extend it.
     # By default run is clickable, we need to disable it until we upload a file.
     shinyjs::disable("run")  #disable run button diable(id of button)
     # we can use observe an event when a certain file uploadded, 
     ## then want to enable run button
     observe({
          if(is.null(input$file) != TRUE){
            shinyjs::enable("run")
          } else {
            shinyjs::disable("run")
          }
     })
      
     observeEvent(input$reset,{
         shinyjs::reset("file")
         ## when we remove the file run button should not be able to work
         shinyjs::diable("run")
     })
     
     observeEvent(input$run, {
       shinyjs::disable("run")
       
       # show loading widget using Rpack
       show_modal_spinner(text = "Now we are preparing for plots")
       # obj <- load_seurat_obj(input$file$datapath) # we need data path for that temp file
       ## so obj will either be seurat object or error vector so we can do simple check
       ## so if its a vector show popup
       
       obj <- load_seurat_obj(input$file$datapath)
       if (is.vector(obj)){
         showModal(modalDialog(
           title = "Error with file",
           HTML("<h5>There is an error with the your file.</h5><br>",
                paste(unlist(obj), collapse = "<br><br>")) ### error msg from load seurat function
         ))
         shinyjs::enable("run")
         
       } else {
         
         ######################################################################
         ##### create renderplot for UMAP ###########
         output$umap <- renderPlot({
           if (!is.null(input$metadata_col)) {
             create_metadata_UMAP(obj, input$metadata_col)
           }
         })
         
         
         ######################################################################
         ######## Create Render plot for gene #########
         output$featurePlot <- renderPlot({
           if (!is.null(input$gene)) {
             create_feature_plot(obj, input$gene)
           }
         })
         
         
         ######### Downloaad plot ###################
         output$downloadFeaturePlot <- downloadHandler(
           filename = function(){
             paste0(input$gene, '_feature_plot', '.png')
           },
           content = function(file){
             plot <- create_feature_plot(obj, input$gene)
             ggsave(filename=file, width = 10, height = 5, type = "cairo")
           }
         )
         output$download_umap <- downloadHandler(
           filename = function(){
             paste0(input$metadata_col, '_UMAP', '.png')
           },
           content = function(file){
             plot <- create_metadata_UMAP(obj, input$metadata_col)
             ggsave(filename=file, width = 10, height = 5, type = "cairo")
           }
         )
         
         output$download_tsne <- downloadHandler(
           filename = function(){
             paste0(input$metadata_col, '_TSNE', '.png')
           },
           content = function(file){
             plot <- create_metadata_TSNE(obj, input$metadata_col)
             ggsave(filename=file, width = 10, height = 5, type = "cairo")
           }
         )
         
         
         
         ######################################################################
         ##### create renderplot for TSNE ###########
         output$tsne <- renderPlot({
           if (!is.null(input$metadata_col)) {
             create_metadata_TSNE(obj, input$metadata_col)
           }
         })
         
         
         ####################################################################
         
         insertTab(
           inputId = "main_tabs",
           tabPanel(
             "UMAP",
             ####### this part will devide two part, chart and dropdown menu
             fluidRow(
               column(
                 width = 8,
                 plotOutput(outputId = 'umap'), ##buildin function plot umap
                 downloadButton("download_umap", "Download UMAP Plot") ## id = download_umap, name
               ),
               column(
                 width = 4,
                 selectizeInput("metadata_col", 
                                "Metadata Column", 
                                colnames(obj@meta.data) 
                 )
               )
             ),
             style = "height: 90%; width: 95%; padding-top: 5%;"
           ),
           select = TRUE
         )
         
         insertTab(
           inputId = "main_tabs",
           tabPanel(
             "Gene Expression",
             fluidRow(
               column(
                 width = 8,
                 plotOutput(outputId = 'featurePlot'),
                 downloadButton("downloadFeaturePlot", "Download Gene Plot")
               ),
               column(
                 width = 4,
                 selectizeInput("gene", 
                                "Genes", 
                                rownames(obj)
                 )
               )
             ),
             style = "height: 90%; width: 95%; padding-top: 5%;"
           )
         )
         ### Stop the spinner once the tabs created
         
#         insertTab(
#           inputId = "main_tabs",
#           tabPanel("t-sne")
#         )
         
         insertTab(
           inputId = "main_tabs",
           tabPanel(
             "TSNE",
             ####### this part will devide two part, chart and dropdown menu
             fluidRow(
               column(
                 width = 8,
                 plotOutput(outputId = 'tsne'), ##buildin function plot tsne
                 downloadButton("download_tsne", "Download TSNE Plot") ## id = download_tsne, name
               ),
               column(
                 width = 4,
                 selectizeInput("metadata_col", 
                                "Metadata Column", 
                                colnames(obj@meta.data) 
                 )
               )
             ),
             style = "height: 90%; width: 95%; padding-top: 5%;"
           ),
           select = TRUE
         )
         
         remove_modal_spinner() ### then eenable run button
         shinyjs::enable("run")
         
       }
     })
     
     # Clear all sidebar inputs when 'Reset' button is clicked
     observeEvent(input$reset, {
       shinyjs::reset("file")
       removeTab("main_tabs", "UMAP")
       removeTab("main_tabs", "Gene Expression Plot")
       shinyjs::disable("run")
     
     })
     
}


shinyApp(ui, server)



