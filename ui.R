library(shiny)
library(plotly)

shinyUI(fluidPage(
  titlePanel("Data Visualization: Missed Dose Analysis"),
  fluidRow(
    
    column(12,
           fluidRow(
             column(6,
                    plotlyOutput("plot1",700,400),
                    textOutput("text1")),
             column(5,
                    plotlyOutput("plot2",600,400),
                    textOutput("text2")
             )
           )
    )),
  #slider widget to zoom
  fluidRow(
    column(10,
           sliderInput("slider", label = h3("adjust the weight scale to zoom"), min = 0, 
                       max = 1000, value = c(0, 400))
    )
  ),
  wellPanel(
    selectInput("Dose", 
                "Select Dose data to highlight",
                choices = c("Dose 1", "Dose 2", "Dose 3", "Dose 4"),
                selected = "Dose 1"),

    br(),
    selectInput("Comp","Choose Subpopulation of Interest:",
                choices = list("Male Age 12 1970","Male Age 12 2000","Female Age 12 1970",
                               "Female Age 12 2000","Male Age 20 1970", "Male Age 20 2000", 
                               "Female Age 12 1970","Female Age 20 1970","Male Age 50 1970",
                               "Male Age 50 2000","Female Age 50 1970","Female Age 50 2000"),selected="Male Age 12 1970") 
  )
    
  )



)
