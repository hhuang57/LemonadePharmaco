# server.R
# load the necessary libraries
library(R.matlab) 
library(gridExtra)
library(ggplot2)
library(reshape)
library(shiny)
library(plotly)
library(reshape2)
library(plyr)

Data <- readMat("data/diffdose.mat", header=T) # Loading in simulated data
df <- as.data.frame(Data) # convert to dataframe
df <- rename(df,c('out.1' = 'Dose.1', 'out.2' = 'Dose.2', 'out.3' = 'Dose.3', 'out.4' = 'Dose.4'))
S1 = subset(df,select=c(time.1,Dose.1)) # Different dataframes per Dose for simulated data
S1 <- rename(S1,c('time.1' = 'time')) # Rename to make column names the same for every Dose subset dataframe
S1 <- melt(S1,id.vars='time')
S2 = subset(df,select=c(time.2,Dose.2))
S2 <- rename(S2,c('time.2' = 'time'))
S2 <- melt(S2,id.vars='time')
S3 = subset(df,select=c(time.3,Dose.3))
S3 <- rename(S3,c('time.3' = 'time'))
S3 <- melt(S3,id.vars='time')
S4 = subset(df,select=c(time.4,Dose.4))
S4 <- rename(S4,c('time.4' = 'time'))
S4 <- melt(S4,id.vars='time')

# Import data.
Data_stats <- readMat("data/mu_sd.mat",header=T) 
df_stats <- as.data.frame(Data_stats) 
Data_subj <- readMat("data/subject.mat",header=T) 
df_subj <- as.data.frame(Data_subj)

#rename subpoulation
df_subj <-rename(df_subj,c('subject.1'='Male Age 12 1970','subject.7'='Male Age 12 2000','subject.4'='Female Age 12 1970',
                           'subject.10'='Female Age 12 2000','subject.2'='Male Age 20 1970','subject.8'='Male Age 20 2000',
                           'subject.5'='Female Age 20 1970','subject.11'='Female Age 20 2000','subject.3'='Male Age 50 1970',
                           'subject.9'='Male Age 50 2000','subject.6'='Female Age 50 1970','subject.12'='Female Age 50 2000'))

df_melt <- melt(df_subj,id.vars="subjectid") 


org_plot <- ggplot(NULL) + # Base plot: only experimental data, no simulated data
  geom_line(data=df,aes(x=time.1,y=Dose.1),color='red', size=0.5,alpha=0.3) +
  geom_line(data=df,aes(x=time.2,y=Dose.2),color='blue',size=0.5,alpha=0.3) +
  geom_line(data=df,aes(x=time.3,y=Dose.3),color='green',size=0.5,alpha=0.3) +
  geom_line(data=df,aes(x=time.4,y=Dose.4),color='yellow',size=0.5,alpha=0.3) +
  scale_x_continuous("Time (hrs)") +
  scale_y_continuous("TFV (nmol/L)") +
  theme(legend.text=element_text(size=10),
        legend.background=element_rect(fill='transparent')) + 
  scale_color_manual(name="",values=c('Dose 1' = "red", 'Dose 2' = "blue", 'Dose 3' = "green", 'Dose 4' = "yellow")) +
  ggtitle("Simulation on Different Doses")


shinyServer(function(input, output) {
  
  output$plot1 <- renderPlotly({ 
    
    #base plot for fig 1
    box <- ggplot(df_melt,aes(x=variable,y=value,fill='grey'),alpha=0.2) + 
      geom_boxplot() +
      scale_x_discrete(breaks=NULL) +
      
      #limit y range based on slider input
      scale_y_continuous("Weight (lbs)",limits=input$slider) +
      ggtitle("Weight Distribution") + 
      guides(fill=FALSE) 
    
    #highlight specified subpopulation
    addd <- df_melt[df_melt$variable %in% input$Comp,]
    boxx<-box
    boxx <- boxx +
      geom_boxplot(data = addd,aes(x=variable,y=value),fill='blue',alpha=0.6)
    p1 <- ggplotly(boxx)
    
  })

  output$plot2 <- renderPlotly({ # plot Figure 2: Optimized Patient Simulations
    dataInput_line <- switch(input$Dose, # Select simulation data subset by highlighted Dose
                             "Dose 1" = S1, "Dose 2" = S2, "Dose 3" = S3,
                             "Dose 4" = S4)

      gg <- org_plot
      gg <- gg +
        geom_line(data=dataInput_line,aes(x=time,y=value),color='red',alpha=0.8,size=0.5)

    
    p <- ggplotly(gg)
  })
  
  output$text1 <- renderText({ # Figure 2 caption
    paste("Dose 1=75 mg, Dose 2=150 mg, Dose 3=300 mg, Dose 4=600 mg. You are currently highlighting data in red for ", input$Dose)
  })
  
  output$text2 <- renderText({ # Figure 2 caption
    paste("Dose 1=75 mg, Dose 2=150 mg, Dose 3=300 mg, Dose 4=600 mg. You are currently highlighting data in red for ", input$Dose)
  })
}) 