### UI for Shiny App: Circle Plots for Karin Meeker 
## Author: Diana A. Hobbs
# August 2022


#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinythemes)
library(chorddiag)

# Define UI for application that draws a histogram
shinyUI(fluidPage(theme=shinytheme("paper"),
  # Application title
  titlePanel("Correlation Matrices"),
  br(),
  br(),
  sidebarLayout(sidebarPanel(verticalLayout(radioButtons('select_group',"Group Category",inline = F,
                                          choices = c("All Samples","CDR=0/A-","CDR=0/A+","CDR>0/A+"),
                                          selected = 'All Samples'),
                                          sliderInput("slider","Threshold",min=.25, max=.75, value=.5, step=.25),
                                          tags$img(src='Legend.png', height="100px", width="100px", alt="Something went wrong."))),#, align = "right"
                mainPanel(chorddiagOutput("distPlot", height = 600)))
  
))

 