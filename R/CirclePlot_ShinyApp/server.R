### Data Processing & Server for Shiny App: Circle Plots for Karin Meeker 
## Author: Diana A. Hobbs
# August 2022

#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

#### Set Up Packages and Load Data ####
pacman::p_load(tidyverse,corrplot,circlize,devtools,chorddiag,htmlwidgets, igraph, tidygraph,shinythemes) 
#devtools::install_github("mattflor/chorddiag", build_vignettes = TRUE,force=TRUE)
library(chorddiag)
library(shiny)

# Load data
#setwd("/Users/hobbsd/Documents/GitHub/KarinMeeker/CirclePlots/")
data <- read.csv("./data_final.csv", header=T)

# Variables for correlation matrices 
vars <- data%>%select(group,22:48)%>%
  rename("Sum of Boxes"=sumbox,
         #MMSE,
         #PRS,
         "Cortical Signature" = cortSig,
         #Centiloid,
         #Tauopathy,
         "Tau/AB40" = CSF_Tau,
         "pTau/AB40" = CSF_pTau,
         "Plasma AB42/40" = PL_Ab42.AB40,
         "AB42/40" = CSF_Ab42.AB40,
         "NfL" = CSF_NfL,
         "Plasma NfL" = PL_NfL,
         "Ng" = CSF_Ng,
         "VILIP-1" = CSF_VILIP1,
         "SNAP-25" = CSF_SNAP25,
         "YKL-40" = CSF_YKL40,
         "sTREM2" = CSF_sTREM2,
         "pT111/T111" = T111,
         "pT153/T153" = T153,
         "pT175/T175" = T175,
         "pT181/T181" = T181,
         "pS199/S199" = S199,
         "pS202/S202" = S202,
         "pT205/T205" = T205,
         "pS208/S208" = S208,
         "pT217/T217" = T217,
         "pT231/T231" = T231)

#groupNames <- c("Sum of Boxes", "MMSE", "PRS", "Cortical Signature", "Centiloid", "Tauopathy",
#                "CSF Tau", "CSF pTau", "Plasma Ab42/40", "CSF Ab42/40", "CSF NfL", "Plasma NfL",
#                "CSF Ng", "CSF VILIP-1", "CSF SNAP-25", "CSF YKL-40", "CSF sTREM2", 
#                "T111", "T153", "T175", "T181", "S199", "S202", "T205", "S208", "T217", "T231")

### Data frames for matrices 
# All Samples  
all <- vars%>%select(-group)

# Group 1: Cognitively Normal, Amyloid Negative "CDR=0/A-"
group1 <- vars%>%filter(group==1)%>%select(-group)

# Group 2: Cognitively Normal, Amyloid Positive "CDR=0/A+" 
group2 <- vars%>%filter(group==2)%>%select(-group)

# Group 3: Cognitively Impaired, Amyloid Positive "CDR>0/A+"
group3 <- vars%>%filter(group==3)%>%select(-group)

### Correlation matrices 
# All Samples 
corr_matrix_all <- cor(na.omit(all), method="spearman", use = "everything")
#corr_matrix_all[upper.tri(corr_matrix_all, diag=TRUE)] <- 0 # Set upper triangle and diagonals to 0
corr_matrix_all<- abs(corr_matrix_all)
corr_matrix_all[corr_matrix_all < 0.5] <- 0 # filter out any correlation that is less than |.5|

# Group 1: Cognitively Normal, Amyloid Negative "CDR=0/A-"
corr_matrix_1 <- cor(na.omit(group1), method="spearman", use = "everything")
#corr_matrix_1[upper.tri(corr_matrix_1, diag=TRUE)] <- 0 # Set upper triangle and diagonals to 0
corr_matrix_1<- abs(corr_matrix_1)
corr_matrix_1[corr_matrix_1 < 0.5] <- 0 # filter out any correlation that is less than |.5|

# Group 2: Cognitively Normal, Amyloid Positive "CDR=0/A+" 
corr_matrix_2 <- cor(na.omit(group2), method="spearman", use = "everything")
c#orr_matrix_2[upper.tri(corr_matrix_2, diag=TRUE)] <- 0 # Set upper triangle and diagonals to 0
corr_matrix_2<- abs(corr_matrix_2)
corr_matrix_2[corr_matrix_2 < 0.5] <- 0 # filter out any correlation that is less than |.5|

# Group 3: Cognitively Impaired, Amyloid Positive "CDR>0/A+"
corr_matrix_3 <- cor(na.omit(group3), method="spearman", use = "everything")
#corr_matrix_3[upper.tri(corr_matrix_3, diag=TRUE)] <- 0 # Set upper triangle and diagonals to 0
corr_matrix_3<- abs(corr_matrix_3)
corr_matrix_3[corr_matrix_3 < 0.5] <- 0 # filter out any correlation that is less than |.5|


### Correlation matrices using a threshold of 0.25
# All Samples 
corr_matrix_all_25 <- cor(na.omit(all), method="spearman", use = "everything")
#corr_matrix_all_25[upper.tri(corr_matrix_all_25, diag=TRUE)] <- 0 # Set upper triangle and diagonals to 0
corr_matrix_all_25<- abs(corr_matrix_all_25)
corr_matrix_all_25[corr_matrix_all_25 < 0.25] <- 0 # filter out any correlation that is less than |.5|

# Group 1: Cognitively Normal, Amyloid Negative "CDR=0/A-"
corr_matrix_1_25 <- cor(na.omit(group1), method="spearman", use = "everything")
#corr_matrix_1_25[upper.tri(corr_matrix_1_25, diag=TRUE)] <- 0 # Set upper triangle and diagonals to 0
corr_matrix_1_25<- abs(corr_matrix_1_25)
corr_matrix_1_25[corr_matrix_1_25 < 0.25] <- 0 # filter out any correlation that is less than |.5|

# Group 2: Cognitively Normal, Amyloid Positive "CDR=0/A+" 
corr_matrix_2_25 <- cor(na.omit(group2), method="spearman", use = "everything")
#corr_matrix_2_25[upper.tri(corr_matrix_2_25, diag=TRUE)] <- 0 # Set upper triangle and diagonals to 0
corr_matrix_2_25<- abs(corr_matrix_2_25)
corr_matrix_2_25[corr_matrix_2_25 < 0.25] <- 0 # filter out any correlation that is less than |.5|

# Group 3: Cognitively Impaired, Amyloid Positive "CDR>0/A+"
corr_matrix_3_25 <- cor(na.omit(group3), method="spearman", use = "everything")
#corr_matrix_3_25[upper.tri(corr_matrix_3_25, diag=TRUE)] <- 0 # Set upper triangle and diagonals to 0
corr_matrix_3_25<- abs(corr_matrix_3_25)
corr_matrix_3_25[corr_matrix_3_25 < 0.25] <- 0 # filter out any correlation that is less than |.5|


### Correlation matrices using a threshold of 0.75
# All Samples 
corr_matrix_all_75 <- cor(na.omit(all), method="spearman", use = "everything")
#corr_matrix_all_75[upper.tri(corr_matrix_all_75, diag=TRUE)] <- 0 # Set upper triangle and diagonals to 0
corr_matrix_all_75<- abs(corr_matrix_all_75)
corr_matrix_all_75[corr_matrix_all_75 < 0.75] <- 0 # filter out any correlation that is less than |.5|

# Group 1: Cognitively Normal, Amyloid Negative "CDR=0/A-"
corr_matrix_1_75 <- cor(na.omit(group1), method="spearman", use = "everything")
#corr_matrix_1_75[upper.tri(corr_matrix_1_75, diag=TRUE)] <- 0 # Set upper triangle and diagonals to 0
corr_matrix_1_75<- abs(corr_matrix_1_75)
corr_matrix_1_75[corr_matrix_1_75 < 0.75] <- 0 # filter out any correlation that is less than |.5|

# Group 2: Cognitively Normal, Amyloid Positive "CDR=0/A+" 
corr_matrix_2_75 <- cor(na.omit(group2), method="spearman", use = "everything")
#corr_matrix_2_75[upper.tri(corr_matrix_2_75, diag=TRUE)] <- 0 # Set upper triangle and diagonals to 0
corr_matrix_2_75<- abs(corr_matrix_2_75)
corr_matrix_2_75[corr_matrix_2_75 < 0.75] <- 0 # filter out any correlation that is less than |.5|

# Group 3: Cognitively Impaired, Amyloid Positive "CDR>0/A+"
corr_matrix_3_75 <- cor(na.omit(group3), method="spearman", use = "everything")
#corr_matrix_3_75[upper.tri(corr_matrix_3_75, diag=TRUE)] <- 0 # Set upper triangle and diagonals to 0
corr_matrix_3_75<- abs(corr_matrix_3_75)
corr_matrix_3_75[corr_matrix_3_75 < 0.75] <- 0 # filter out any correlation that is less than |.5|

### Grouping by color
#FFFFFF = White 
#D0ECFB = Blue (Core AD)
#FFC6B2 = Red (Neuronal Dysfunction)
#F9EDAD = Yellow (Non-specifc Tau)
#7FE0DF = Teal (Neurodegeneration)

allColors <- c("#D0ECFB", # Sumbox
               "#FFFFFF","#FFFFFF", # MMSE, PRS
               "#7FE0DF", # cortSig
               "#D0ECFB","#D0ECFB","#D0ECFB","#D0ECFB","#D0ECFB","#D0ECFB", # Centiloid, Tauopathy, CSF Tau, CSF pTau, Plasma Ab42/40, CSF Ab42/40
               "#7FE0DF","#7FE0DF", # CSF NfL, Plasma NfL
               "#FFC6B2","#FFC6B2","#FFC6B2","#FFC6B2","#FFC6B2", # CSF Ng, CSF VILIP-1, CSF SNAP-25, CSF YLK-40, CSF sTREM2
               "#D0ECFB","#D0ECFB", # T111, T153
               "#FFFFFF", # T175
               "#D0ECFB", # T181
               "#F9EDAD","#F9EDAD","#F9EDAD", # S199, S202, T205
               "#D0ECFB","#D0ECFB","#D0ECFB") # S208, T217, T231

group1Colors <- c("#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF", # Sumbox, MMSE, PRS, cortSig
                  "#D0ECFB", # Centiloid
                  "#FFFFFF", # Tauopathy
                  "#D0ECFB","#D0ECFB","#D0ECFB","#D0ECFB", # CSF Tau, CSF pTau, Plasma Ab42/40, CSF Ab42/40
                  "#7FE0DF","#7FE0DF", # CSF NfL, Plasma NfL
                  "#FFC6B2","#FFC6B2","#FFC6B2","#FFC6B2","#FFC6B2", # CSF Ng, CSF VILIP-1, CSF SNAP-25, CSF YLK-40, CSF sTREM2
                  "#D0ECFB","#D0ECFB", # T111, T153
                  "#FFFFFF", # T175                
                  "#D0ECFB", # T181
                  "#F9EDAD","#F9EDAD","#F9EDAD", # S199, S202, T205
                  "#D0ECFB","#D0ECFB","#D0ECFB") # S208, T217, T231

group2Colors <- c("#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF", # Sumbox, MMSE, PRS, cortSig
                  "#D0ECFB","#D0ECFB","#D0ECFB","#D0ECFB","#D0ECFB","#D0ECFB", # Centiloid, Tauopathy, CSF Tau, CSF pTau, Plasma Ab42/40, CSF Ab42/40
                  "#7FE0DF","#7FE0DF", # CSF NfL, Plasma NfL
                  "#FFC6B2","#FFC6B2","#FFC6B2","#FFC6B2","#FFC6B2", # CSF Ng, CSF VILIP-1, CSF SNAP-25, CSF YLK-40, CSF sTREM2
                  "#D0ECFB","#D0ECFB", # T111, T153
                  "#FFFFFF", # T175                
                  "#D0ECFB", # T181
                  "#F9EDAD","#F9EDAD","#F9EDAD", # S199, S202, T205
                  "#D0ECFB","#D0ECFB","#D0ECFB") # S208, T217, T231

group3Colors <- c("#FFFFFF","#FFFFFF","#FFFFFF", # Sumbox, MMSE, PRS, 
                  "#7FE0DF", # cortSig
                  "#D0ECFB","#D0ECFB","#D0ECFB","#D0ECFB", # Centiloid, Tauopathy, CSF Tau, CSF pTau
                  "#7FE0DF", # Plasma Ab42/40
                  "#D0ECFB", #CSF Ab42/40
                  "#7FE0DF","#7FE0DF", # CSF NfL, Plasma NfL
                  "#FFC6B2","#FFC6B2","#FFC6B2","#FFC6B2", # CSF Ng, CSF VILIP-1, CSF SNAP-25, CSF YLK-40
                  "#FFFFFF", # CSF sTREM2
                  "#D0ECFB","#D0ECFB", # T111, T153
                  "#FFFFFF", # T175   
                  "#D0ECFB", # T181
                  "#FFFFFF","#FFFFFF", # S199, S202
                  "#D0ECFB","#D0ECFB","#D0ECFB","#D0ECFB") # 205, S208, T217, T231


shinyServer(function(input, output) {
  
  
  output$distPlot <- renderChorddiag({
    
    if(input$select_group =="All Samples" && input$slider == 0.5){
      chorddiag(corr_matrix_all,
                type="directional", width=850, height=850,
                groupnamePadding = 20,
                groupPadding = 3,
                groupColors = allColors,
                #groupNames = groupNames,
                groupThickness=0.05,
                groupnameFontsize = 14,
                showTicks = F,
                showTooltips = T,
                margin = 120,
                #chordedgeColor = "#808080",
                tooltipGroupConnector = "    &#x25B6;    ")
    }else if(
      input$select_group =="All Samples" && input$slider == 0.25){
      chorddiag(corr_matrix_all_25,
                type="directional", width=850, height=850,
                groupnamePadding = 20,
                groupPadding = 3,
                groupColors = allColors,
                #groupNames = groupNames,
                groupThickness=0.05,
                groupnameFontsize = 14,
                showTicks = F,
                showTooltips = T,
                margin = 120,
                #chordedgeColor = "#808080",
                tooltipGroupConnector = "    &#x25B6;    ")
    }else if(
      input$select_group =="All Samples" && input$slider == 0.75){
      chorddiag(corr_matrix_all_75,
                type="directional", width=850, height=850,
                groupnamePadding = 20,
                groupPadding = 3,
                groupColors = allColors,
                #groupNames = groupNames,
                groupThickness=0.05,
                groupnameFontsize = 14,
                showTicks = F,
                showTooltips = T,
                margin = 120,
                #chordedgeColor = "#808080",
                tooltipGroupConnector = "    &#x25B6;    ") 
     }else if(
      input$select_group =="CDR=0/A-" && input$slider == 0.5){
      chorddiag(corr_matrix_1,
                type="directional", width=850, height=850,
                groupnamePadding = 20,
                groupPadding = 3,
                groupColors = group1Colors,
                #groupNames = groupNames,
                groupThickness=0.05,
                groupnameFontsize = 14,
                showTicks = F,
                showTooltips = T,
                margin = 120,
                #chordedgeColor = "#808080",
                tooltipGroupConnector = "    &#x25B6;    ")
     }else if(
       input$select_group =="CDR=0/A-" && input$slider == 0.25){
       chorddiag(corr_matrix_1_25,
                 type="directional", width=850, height=850,
                 groupnamePadding = 20,
                 groupPadding = 3,
                 groupColors = group1Colors,
                 #groupNames = groupNames,
                 groupThickness=0.05,
                 groupnameFontsize = 14,
                 showTicks = F,
                 showTooltips = T,
                 margin = 120,
                 #chordedgeColor = "#808080",
                 tooltipGroupConnector = "    &#x25B6;    ")
     }else if(
       input$select_group =="CDR=0/A-" && input$slider == 0.75){
       chorddiag(corr_matrix_1_75,
                 type="directional", width=850, height=850,
                 groupnamePadding = 20,
                 groupPadding = 3,
                 groupColors = group1Colors,
                 #groupNames = groupNames,
                 groupThickness=0.05,
                 groupnameFontsize = 14,
                 showTicks = F,
                 showTooltips = T,
                 margin = 120,
                 #chordedgeColor = "#808080",
                 tooltipGroupConnector = "    &#x25B6;    ")
    } else if(
      input$select_group =="CDR=0/A+" && input$slider == 0.5){
      chorddiag(corr_matrix_2,
                type="directional", width=850, height=850,
                groupnamePadding = 20,
                groupPadding = 3,
                groupColors = group2Colors,
                #groupNames = groupNames,
                groupThickness=0.05,
                groupnameFontsize = 14,
                showTicks = F,
                showTooltips = T,
                margin = 120,
                #chordedgeColor = "#808080",
                tooltipGroupConnector = "    &#x25B6;    ")
    } else if(
      input$select_group =="CDR=0/A+" && input$slider == 0.25){
      chorddiag(corr_matrix_2_25,
                type="directional", width=850, height=850,
                groupnamePadding = 20,
                groupPadding = 3,
                groupColors = group2Colors,
                #groupNames = groupNames,
                groupThickness=0.05,
                groupnameFontsize = 14,
                showTicks = F,
                showTooltips = T,
                margin = 120,
                #chordedgeColor = "#808080",
                tooltipGroupConnector = "    &#x25B6;    ")
    } else if(
      input$select_group =="CDR=0/A+" && input$slider == 0.75){
      chorddiag(corr_matrix_2_75,
                type="directional", width=850, height=850,
                groupnamePadding = 20,
                groupPadding = 3,
                groupColors = group2Colors,
                #groupNames = groupNames,
                groupThickness=0.05,
                groupnameFontsize = 14,
                showTicks = F,
                showTooltips = T,
                margin = 120,
                #chordedgeColor = "#808080",
                tooltipGroupConnector = "    &#x25B6;    ")
    } else if(
      input$select_group =="CDR>0/A+" && input$slider == 0.25){
      chorddiag(corr_matrix_3_25,
                type="directional", width=850, height=850,
                groupnamePadding = 20,
                groupPadding = 3,
                groupColors = group3Colors,
                #groupNames = groupNames,
                groupThickness=0.05,
                groupnameFontsize = 14,
                showTicks = F,
                showTooltips = T,
                margin = 120,
                #chordedgeColor = "#808080",
                tooltipGroupConnector = "    &#x25B6;    ")
    } else if(
      input$select_group =="CDR>0/A+" && input$slider == 0.75){
      chorddiag(corr_matrix_3_75,
                type="directional", width=850, height=850,
                groupnamePadding = 20,
                groupPadding = 3,
                groupColors = group3Colors,
                #groupNames = groupNames,
                groupThickness=0.05,
                groupnameFontsize = 14,
                showTicks = F,
                showTooltips = T,
                margin = 120,
                #chordedgeColor = "#808080",
                tooltipGroupConnector = "    &#x25B6;    ")
      }else{
      chorddiag(corr_matrix_3,
                type="directional", width=850, height=850,
                groupnamePadding = 20,
                groupPadding = 3,
                groupColors = group3Colors,
                #groupNames = groupNames,
                groupThickness=0.05,
                groupnameFontsize = 14,
                showTicks = F,
                showTooltips = T,
                margin = 120,
                #chordedgeColor = "#808080",
                tooltipGroupConnector = "    &#x25B6;    ")
    }
    
  })
  
})

