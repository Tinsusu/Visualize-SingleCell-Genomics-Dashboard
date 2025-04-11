
source('global.R')
ui <- dashboardPage(
  dashboardHeader(title = "scRNAseq Plot"),
  dashboardSidebar(
      sidebarMenu(id='tab',
                useShinyjs(),
                menuItem("Kick Things Off Here!", tabName = "home", icon = icon("dna")),
                menuItem("scRNAseq Investigator", tabName = "input", icon = icon("dot-circle")),
                conditionalPanel(condition = "input.tab == 'input'",
                      div(
                          fileInput("file", "Please Upload File", multiple=TRUE, accept=c('.rds')),
                          actionButton("reset", "Reset", icon = icon("undo"), 
                                       style = "color: #fff; background-color: #C69C6D; border-radius: 50%; width: 87.25%; padding: 15px; font-size: 18px; text-align: center;"),
                          actionButton("run", "Run", icon = icon("play"), 
                                       style = "color: #fff; background-color: #77DD77; width: 87.25%; padding: 15px; font-size: 18px; text-align: center;")
                        )
            )
    )
  ), 
  
  dashboardBody(
    tabItems(
      tabItem(tabName = "input", # tabItem refers to tab in sidebar (not main panel)
              tabsetPanel(id = 'main_tabs',
                          tabPanel("Follow Me please",
                                   includeMarkdown("./markdown/instructions.md")
                          )
              )
      ),
      tabItem(tabName = "home",
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

server <- function(input, output, session) {
     options(shiny.maxRequestSize = 300*1024^2) #shiny has a limit size, so do this to extend it.
     
     shinyjs::disable("run") #disable run button 
     
     # we can use observe an event when a certain file uploadded, want to enable run button
     observe({
          if(is.null(input$file) != TRUE){
            shinyjs::enable("run")
          } else {
            shinyjs::disable("run")
          }
     })
     
     observeEvent(input$reset,{
         shinyjs::reset("file")
         shinyjs::diable("run")
     })
     
}

shinyApp(ui, server)



