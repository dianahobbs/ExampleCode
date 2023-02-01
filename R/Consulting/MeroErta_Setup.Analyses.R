#### Mero/Erta data project 
### Setting up the data and running analyses
## Author: Diana Hobbs
## Date: May 2022

#### Set Up Packages and Load Data ####
pacman::p_load(tidyverse,jtools,psych,cowplot,multcomp,gmodels,lubridate,gtsummary,ggplot2)

data <- read.csv("./OriginalData.csv")%>%
  rename(PrePost=X2020...Pre.intervention)%>%
  mutate(PrePost = str_replace_all(PrePost, c("0"="2022", "1"="2020")),
         Sex = str_replace_all(Sex, c("Male"="0", "Female"="1")),
         University = str_replace_all(University, c("0"="Other", "1"="University")))

Mero <- data%>%
  filter(Drug=="meropenem")
Erta <- data%>%
  filter(Drug=="ertapenem")
Erta_0 <- Erta%>%
  filter(OnceOrder==0)
Erta_1 <- Erta%>%
  filter(OnceOrder==1)

table(data$University)
table(Mero$University)
table(Erta$University)
table(Erta_0$University)
table(Erta_1$University)


table(data$University,data$PrePost)
table(Mero$University,Mero$PrePost)
table(Erta$University,Erta$PrePost)
table(Erta_0$University,Erta_0$PrePost)
table(Erta_1$University,Erta_1$PrePost)

CrossTable(data$University, data$PrePost)
chisq.test(data$University, data$PrePost)

CrossTable(Mero$University, Mero$PrePost)
chisq.test(Mero$University, Mero$PrePost)

CrossTable(Erta$University, Erta$PrePost)
chisq.test(Erta$University, Erta$PrePost)

CrossTable(Erta_0$University, Erta_0$PrePost)
chisq.test(Erta_0$University, Erta_0$PrePost)
CrossTable(Erta_1$University, Erta_1$PrePost)
chisq.test(Erta_1$University, Erta_1$PrePost)

table(data$Location,data$PrePost)
table(Mero$Location,Mero$PrePost)
table(Erta$Location,Erta$PrePost)
table(Erta_0$Location,Erta_0$PrePost)
table(Erta_1$Location,Erta_1$PrePost)

CrossTable(data$Location, data$PrePost)
chisq.test(data$Location, data$PrePost)

CrossTable(Mero$Location, Mero$PrePost)
chisq.test(Mero$Location, Mero$PrePost)

CrossTable(Erta$Location, Erta$PrePost)
chisq.test(Erta$Location, Erta$PrePost)

CrossTable(Erta_0$Location, Erta_0$PrePost)
chisq.test(Erta_0$Location, Erta_0$PrePost)
CrossTable(Erta_1$Location, Erta_1$PrePost)
chisq.test(Erta_1$Location, Erta_1$PrePost)

# Create Table
theme_gtsummary_compact()
#data$MMSE_TauCDR1 <- as.numeric(data$MMSE_TauCDR1)

All <- data%>%
  select(Age, Sex, Location, University, PrePost, Drug, LOS, OnceOrder)%>%
  tbl_summary(
    statistic=list(all_continuous()~ "{mean}({sd})",
                   all_categorical()~"{n}({p}%)"),
    digits = all_continuous()~2,
    label = c(Age ~ "Age, years", Sex ~ "Sex, M:F", Location ~ "Location", 
              University ~ "Intervention Hospital", PrePost ~ "Year", Drug = "Drug", 
              LOS = "LOS", OnceOrder = "OnceOrder"))%>%
  bold_labels()%>%
  italicize_levels()%>%
  as_gt()%>%
  gt::gtsave(filename = "./All.png")
All

Split_Location <- data%>%
  select(Age, Sex, Location, University, PrePost, Drug, LOS, OnceOrder)%>%
  tbl_summary(by=Location,
    statistic=list(all_continuous()~ "{mean}({sd})",
                   all_categorical()~"{n}({p}%)"),
    digits = all_continuous()~2,
    label = c(Age ~ "Age, years", Sex ~ "Sex, M:F", #Location ~ "Location", 
              University ~ "Intervention Hospital", PrePost ~ "Year", Drug = "Drug", 
              LOS = "LOS", OnceOrder = "OnceOrder"))%>%
  bold_labels()%>%
  italicize_levels()%>%
  as_gt()%>%
  gt::gtsave(filename = "./Split_Location.png")
Split_Location

Split_Location_Erta_1 <- Erta_1%>%
  select(Age, Sex, Location, University, PrePost, Drug, LOS, OnceOrder)%>%
  tbl_summary(by=Location,
              statistic=list(all_continuous()~ "{mean}({sd})",
                             all_categorical()~"{n}({p}%)"),
              digits = all_continuous()~2,
              label = c(Age ~ "Age, years", Sex ~ "Sex, M:F", #Location ~ "Location", 
                        University ~ "Intervention Hospital", PrePost ~ "Year", Drug = "Drug", 
                        LOS = "LOS", OnceOrder = "OnceOrder"))%>%
  bold_labels()%>%
  italicize_levels()%>%
  as_gt()%>%
  gt::gtsave(filename = "./Split_Location_Erta_1.png")
Split_Location_Erta_1


# Basic line plot with points
ggplot(data=data, aes(x=PrePost, y=LOS, group=University)) +
  geom_line(aes(linetype=University))
  #eom_point()
ggplot(data, aes(factor(University), LOS))+
geom_violin(aes(fill = factor(PrePost)))

ggplot(data, aes(factor(Location), LOS))+
  geom_violin(aes(fill = factor(PrePost)))

# Stacked barplot with multiple groups
ggplot(data=data, aes(x=University, y=LOS, fill=PrePost)) +
  geom_bar(stat="identity", position=position_dodge())
ggplot(data=data, aes(x=Location, y=LOS, fill=PrePost)) +
  geom_bar(stat="identity", position=position_dodge())
  
