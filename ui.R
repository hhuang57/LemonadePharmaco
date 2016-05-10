library(shiny)
library(plotly)

shinyUI(fluidPage(
  titlePanel("Data Visualization: Missed Dose Analysis"),
  
  wellPanel(
    selectInput("Dose", 
                "Select Dose data to highlight",
                choices = c("Dose 1", "Dose 2", "Dose 3", "Dose 4"),
                selected = "Dose 1"),

    br()
    
    
  ),
  fluidRow(
    column(8,
           plotlyOutput("plot",600,400),
           textOutput("text")
    )
  )
))