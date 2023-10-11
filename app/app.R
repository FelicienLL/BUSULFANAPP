library(shiny)
library(rhandsontable)
library(mapbayr)
library(mrgsolve)
source("funs.R")
ui <- fluidPage(
  titlePanel("Model-informed precision dosing of Busulfan (Thai pediatric population)"),
  sidebarLayout(
    sidebarPanel(
      width = 4, 
      h4("Patient"),
      fluidRow(
        column(width = 6, numericInput("BW", label = "Body weight (kg)", value = 25)),
        column(width = 12, numericInput("AUCTARGET", label = "Cumulative target AUC (micromolar.h)", value = 16000))
      ), 
      h4("Dosing"),
      h6("Dose unit: mg."),
      fluidRow(
        rHandsontableOutput("dosingEdit")
      ), 
      h4("TDM Concentrations"),
      h6("Concentration unit: ng/mL"),
      fluidRow(
        rHandsontableOutput("concentrationEdit")
      ),
      h4("Estimate!"), 
      actionButton("go", "Go!")
    ), 
    mainPanel(
      width = 8,
      fluidRow(
        tabsetPanel(
          id = "tabset",
          tabPanel(
            title = "Data - table", 
            tableOutput("dataTable")
          ), 
          tabPanel(
            title = "Concentration vs time - plot", 
            plotOutput("plotest")
          )
        )
      ), 
      fluidRow(
        tabPanel(
          title = "Dose calculation - table", 
          tableOutput("repost")
        )
      )
    ) 
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  MWbusulfan <- 246.304
  mod <- mread("busulfan_shiny.cpp")
  
  dosingData <- reactive({
    if(is.null(input$dosingEdit)){
      DF <- data.frame(
        DAY = 0:4,
        Date = format(seq(from = Sys.Date(), by = 1, length.out = 5), "%d/%m/%Y"), 
        Dose = rep(0, 5), 
        Hour_start = rep(NA_character_, 5), 
        Hour_end = rep(NA_character_, 5)
      )
    } else {
      DF <- input$dosingEdit %>% 
        hot_to_r()
    }   
    DF
  })
  
  output$dosingEdit <- renderRHandsontable({
    rhandsontable(dosingData()) %>% 
      hot_col(col = "Date", type = "date", dateFormat = "DD/MM/YYYY")
  })
  
  concentrationData <- reactive({
    if(is.null(input$concentrationEdit)){
      DF <- data.frame(
        DAY = rep(c(1L,3L), each = 3), 
        Date = format(rep(c(Sys.Date() + 1, Sys.Date() + 3), each = 3), "%d/%m/%Y"), 
        Hour_sample = rep(NA_character_, 6),
        Concentration = rep(NA_real_, 6), 
        Exclude = rep(FALSE, 6)
      )
    } else {
      DF <- input$concentrationEdit %>% 
        hot_to_r()
    }   
    DF
  })
  
  output$concentrationEdit <- renderRHandsontable({
    rhandsontable(concentrationData()) %>% 
      hot_col(col = "Date", type = "date", dateFormat = "DD/MM/YYYY")
  })
  
  nmtranData <- reactive({
    make_data(ddata = dosingData(), cdata = concentrationData(), bw = input$BW)
  })
  
  re_est <- eventReactive(input$go, {
    mapbayest(mod, nmtranData())
  })
  
  re_post <- reactive({
    shiny::req(re_est())
    make_post(ddata = dosingData(), cdata = concentrationData(), est = re_est(), auc = input$AUCTARGET)
  })
  
  re_plot <- reactive({
    shiny::req(re_est())
    make_plot(est = re_est(), post = re_post())
  })
  
  output$dataTable <- renderTable({
    nmtranData() %>% 
      dplyr::mutate(.datehour = format(.datehour,'%d-%m-%Y %H:%M:%S'))
  })
  
  output$plotest <- renderPlot({
    re_plot()
  })
  
  output$repost <- renderTable({
    re_post() %>% 
      dplyr::select(-dplyr::starts_with("ETA"))
  })
  
  observeEvent(input$go, {
    updateTabsetPanel(inputId = "tabset", selected = "Concentration vs time - plot")
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
