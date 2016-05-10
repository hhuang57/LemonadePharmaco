library(R.matlab) 
library(gridExtra)
library(ggplot2)
library(reshape)

# Problem 1 Part b: Artesunate Concentrations Over Time
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

org_plot <- ggplot(NULL) + # Base plot: only experimental data, no simulated data
  geom_line(data=df,aes(x=time.1,y=Dose.1),color='red', size=0.5,alpha=0.5) +
  geom_line(data=df,aes(x=time.2,y=Dose.2),color='blue',size=0.5,alpha=0.5) +
  geom_line(data=df,aes(x=time.3,y=Dose.3),color='green',size=0.5,alpha=0.5) +
  geom_line(data=df,aes(x=time.4,y=Dose.4),color='yellow',size=0.5,alpha=0.5) +
  scale_x_continuous("Time (hrs)") +
  scale_y_continuous("TFV (nmol/L)") +
  #theme(legend.text=element_text(size=10),
  #      legend.background=element_rect(fill='transparent')) + 
  scale_color_manual(name="",values=c('Dose 1' = "red", 'Dose 2' = "blue", 'Dose 3' = "green", 'Dose 4' = "yellow")) +
  ggtitle("Simulation on Different Doses")

print(org_plot)