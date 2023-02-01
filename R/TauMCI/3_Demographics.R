# Title: Demographics for ADRC Cohort for Tau/MCI Replication Study (Sylvia/Cherie)
# Changing the data to include parahipp in the metaROI
# Author: Diana Hobbs
# Date: May 2022

#Functions to set up environment
rm(list=ls()) #removes all variables from the current environment

# Load packages 
pacman::p_load(tidyverse,jtools,naniar,lubridate,janitor,lme4,sfsmisc,foreign,MASS,scales,ggseg,cowplot,gtsummary,ggplot2,RColorBrewer,webshot)
#webshot::install_phantomjs(force=T)

# Load data
data<-read.csv("./Data/N151_Wide_5-25-22.csv", header=TRUE)

###### Descriptive Statistics ######
theme_gtsummary_compact()
data$MMSE_TauCDR1 <- as.numeric(data$MMSE_TauCDR1)
#data$ETHNIC <- as.factor(data$ETHNIC)
#data$RACE <- as.factor(data$RACE)


demo <- data%>%
  select(3:14, -Tau_EC, -Tau_IT, -Tracer, -AmyPositivity, MMSE_TauCDR1, PACC_BaseCog1, ETHNIC, RACE)%>%
  tbl_summary(#type = list(MMSE_TauCDR1 ~ "continuous"),
              statistic=list(all_continuous()~ "{mean}({sd})",
                             all_categorical()~"{n}({p}%)"),
              digits = all_continuous()~2,
              label = c(AGE ~ "Age, years", SEX ~ "Sex,F:M", EDUC ~ "Education, years", 
                        ETHNIC ~ "Ethnicity", RACE ~ "Race", APOE4 = "APOE ɛ4 carriers, n", 
                        Centiloid = "Global Aβ SUVR/DVR", Tau_metaROI = "Temporal meta-ROI SUVR", 
                        PercentHIPP = "Hippocampal volume (% of TIV)", Cort_Thick = "Temporal cortical thickness",
                        MMSE_TauCDR1 = "MMSE (/30)", PACC_BaseCog1 = "PACC-5, baseline)"))%>%
  bold_labels()%>%
  italicize_levels()%>%
  as_gt()%>%             
  gt::gtsave(filename = "./5_25_22/AllDemographics_N151_6-14-22.png")
demo
mean(data$MMSE_TauCDR1)
sd(data$MMSE_TauCDR1)
table(data$ETHNIC)
mean(data$EDUC)
sd(data$EDUC)

split_demo <- data%>%
  select(3:14, -Tau_EC, -Tau_IT, -Tracer, -AmyPositivity,  MMSE_TauCDR1, PACC_BaseCog1, Biomarker_group_metaROI, RACE,ETHNIC)%>% #, 
  tbl_summary(by=Biomarker_group_metaROI,
              #type = list(MMSE_TauCDR1 ~ "continuous"),
              statistic=list(all_continuous()~ "{mean}({sd})",
                             all_categorical()~"{n}/{N}({p}%)"),
              digits = all_continuous()~2,
              label = c(AGE ~ "Age, years", SEX ~ "Sex,F:M", EDUC ~ "Education, years", APOE4 = "APOE ɛ4 carriers, n", 
                        ETHNIC ~ "Ethnicity", RACE ~ "Race", Centiloid = "Global Aβ SUVR/DVR", Tau_metaROI = "Temporal meta-ROI SUVR", 
                        PercentHIPP = "Hippocampal volume (% of TIV)", Cort_Thick = "Temporal cortical thickness",
                        MMSE_TauCDR1 = "MMSE (/30)", PACC_BaseCog1 = "PACC-5, baseline)"))%>% 
  bold_labels()%>%
  italicize_levels()%>%
  add_p()%>%
  as_gt()%>%             
  gt::gtsave(filename = "./5_25_22/AllDemographics_BiomarkerSplit_N151_6-14-22.png")
split_demo

mean(data$MMSE_TauCDR1[data$Biomarker_group_metaROI==1])
sd(data$MMSE_TauCDR1[data$Biomarker_group_metaROI==1])
table(data$ETHNIC[data$Biomarker_group_metaROI==1])

mean(data$MMSE_TauCDR1[data$Biomarker_group_metaROI==2])
sd(data$MMSE_TauCDR1[data$Biomarker_group_metaROI==2])
table(data$ETHNIC[data$Biomarker_group_metaROI==2])

mean(data$MMSE_TauCDR1[data$Biomarker_group_metaROI==3])
sd(data$MMSE_TauCDR1[data$Biomarker_group_metaROI==3])
table(data$ETHNIC[data$Biomarker_group_metaROI==3])

mean(data$MMSE_TauCDR1[data$Biomarker_group_metaROI==4])
sd(data$MMSE_TauCDR1[data$Biomarker_group_metaROI==4])
table(data$ETHNIC[data$Biomarker_group_metaROI==4])



NCvC_Biomarker_Groups <- data%>%
  select(Converter_after_tauPET, Biomarker_group_metaROI, Biomarker_group_EC, Biomarker_group_IT, Biomarker_group_ANY)%>%
  tbl_summary(by=Converter_after_tauPET,
              statistic=list(all_continuous()~ "{mean}({sd})",
                             all_categorical()~"{n}/{N}({p}%)"),
              digits = all_continuous()~2,
              label = c(Biomarker_group_metaROI="Biomarker group meta ROI", 
                        Biomarker_group_EC="Biomarker group Entorhinal",
                        Biomarker_group_IT="Biomarker group Inferior Temporal",
                        Biomarker_group_ANY="Biomarker group ANY"))%>%
  bold_labels()%>%
  italicize_levels()%>%
  add_p()%>%
  as_gt()%>%             
  gt::gtsave(filename = "./5_25_22/Nonconvert-v-Convert_Biomarker_N151_6-14-22.png")
NCvC_Biomarker_Groups

# converters and non-converters 1: A+T+, 2: A+T-, 3:A-T+, 4: A-T-
Biomarker_metaROI<-data.frame(tabyl(data, Converter_after_tauPET, Biomarker_group_metaROI))
colnames(Biomarker_metaROI)=c("Converter_after_tauPET", "A+T+", "A+T-", "A-T+", "A-T-")
rownames(Biomarker_metaROI)=c("Non-converter", "Converter")
write.csv(Biomarker_metaROI, "./5_25_22/Biomarker_metaROI_N151_6-14-22.csv")

## Reviewer revision
#data1<-read.csv("./Data/N151_Wide_5-25-22.csv", header=TRUE)
#Biomarker_metaROI1<-data.frame(tabyl(data1, Converter_after_tauPET, Biomarker_group_metaROI))
#colnames(Biomarker_metaROI1)=c("Converter_after_tauPET", "A+T+", "A+T-", "A-T+", "A-T-")
#rownames(Biomarker_metaROI1)=c("Non-converter", "Converter")
#write.csv(Biomarker_metaROI1, "./5_25_22/Biomarker_metaROI_N151_5-25-22.csv")

Biomarker_EC<-data.frame(tabyl(data, Converter_after_tauPET, Biomarker_group_EC))
colnames(Biomarker_EC)=c("Converter_after_tauPET", "A+T+", "A+T-", "A-T+", "A-T-")
rownames(Biomarker_EC)=c("Non-converter", "Converter")
write.csv(Biomarker_EC, "./5_25_22/Biomarker_EC_N151_6-14-22.csv")

Biomarker_IT<-data.frame(tabyl(data, Converter_after_tauPET, Biomarker_group_IT))
colnames(Biomarker_IT)=c("Converter_after_tauPET", "A+T+", "A+T-", "A-T+", "A-T-")
rownames(Biomarker_IT)=c("Non-converter", "Converter")
write.csv(Biomarker_IT, "./5_25_22/Biomarker_IT_N151_6-14-22.csv")

Biomarker_ANY<-data.frame(tabyl(data, Converter_after_tauPET, Biomarker_group_ANY))
colnames(Biomarker_ANY)=c("Converter_after_tauPET", "A+T+", "A+T-", "A-T+", "A-T-")
rownames(Biomarker_ANY)=c("Non-converter", "Converter")
write.csv(Biomarker_ANY, "./5_25_22/Biomarker_ANY_N151_6-14-22.csv")

Biomarker_NeuroPos_HippVol<-data.frame(tabyl(data, Converter_after_tauPET, Biomarker_group_metaROI_HippVol))
colnames(Biomarker_NeuroPos_HippVol)=c("Converter_after_tauPET", "A+T+N+", "A+T+N-", "A+T-N+", "A+T-N-",
                        "A-T+N+", "A-T+N-", "A-T-N+", "A-T-N-")
rownames(Biomarker_NeuroPos_HippVol)=c("Non-converter", "Converter")
write.csv(Biomarker_NeuroPos_HippVol, "./5_25_22/Biomarker_NeuroPos_HippVol_N151_6-14-22.csv")
# Reviewer revision - there are NO 5's 
#Biomarker_NeuroPos_HippVol1<-data.frame(tabyl(data1, Converter_after_tauPET, Biomarker_group_metaROI_HippVol))
#colnames(Biomarker_NeuroPos_HippVol1)=c("Converter_after_tauPET", "A+T+N+", "A+T+N-", "A+T-N+", "A+T-N-",
#                                        "A-T+N+", "A-T+N-", "A-T-N+", "A-T-N-") 
#rownames(Biomarker_NeuroPos_HippVol1)=c("Non-converter", "Converter")
#write.csv(Biomarker_NeuroPos_HippVol1, "./5_25_22/Biomarker_NeuroPos_HippVol_N151_5-25-22.csv")


Biomarker_NeuroPos_CortThick<-data.frame(tabyl(data, Converter_after_tauPET, Biomarker_group_metaROI_CortThick))
colnames(Biomarker_NeuroPos_CortThick)=c("Converter_after_tauPET", "A+T+N+", "A+T+N-", "A+T-N+", "A+T-N-",
                                       "A-T+N+", "A-T+N-", "A-T-N+", "A-T-N-")
rownames(Biomarker_NeuroPos_CortThick)=c("Non-converter", "Converter")
write.csv(Biomarker_NeuroPos_CortThick, "./5_25_22/Biomarker_NeuroPos_CortThick_N151_6-14-22.csv")
# Reviewer revision
#Biomarker_NeuroPos_CortThick1<-data.frame(tabyl(data1, Converter_after_tauPET, Biomarker_group_metaROI_CortThick))
#colnames(Biomarker_NeuroPos_CortThick1)=c("Converter_after_tauPET", "A+T+N+", "A+T+N-", "A+T-N+", "A+T-N-",
#                                         "A-T+N+", "A-T+N-", "A-T-N+", "A-T-N-")
#rownames(Biomarker_NeuroPos_CortThick1)=c("Non-converter", "Converter")
#write.csv(Biomarker_NeuroPos_CortThick1, "./5_25_22/Biomarker_NeuroPos_CortThick_N151_5-25-22.csv")

# Reviewer revision
Biomarker_Race<-data.frame(tabyl(data, RACE, Biomarker_group_metaROI))
colnames(Biomarker_Race)=c("Race", "A+T+", "A+T-","A-T+", "A-T-")
#rownames(Biomarker_NeuroPos_CortThick1)=c("Non-converter", "Converter")
write.csv(Biomarker_Race, "./5_25_22/Biomarker_Race_N151_6-14-22.csv")

#Ensure variables are in correct format
data$Cox_followuptime <- as.numeric(data$Cox_followuptime)
data$Biomarker_group_metaROI <- as.factor(data$Biomarker_group_metaROI)
data$Biomarker_group_EC <- as.factor(data$Biomarker_group_EC)
data$Biomarker_group_IT <- as.factor(data$Biomarker_group_IT)
data$Biomarker_group_ANY <- as.factor(data$Biomarker_group_ANY)
data$SEX <- as.factor(data$SEX)
data$APOE4 <-as.factor(data$APOE4)
data$Converter_after_tauPET <- as.factor(data$Converter_after_tauPET)

#a = significant difference between A+T+ and A+T- groups, group12
#b = significant difference between A+T+ and A-T-groups, group14
#c = significant difference between A+T- and A-T- groups, group24
#1: A+T+, 2: A+T-, EXCLUDE 3:A-T+, 4: A-T-
data_group12<-data[which(data$Biomarker_group_metaROI==c("1","2")),]
data_group14<-data[which(data$Biomarker_group_metaROI==c("1","4")),]
data_group24<-data[which(data$Biomarker_group_metaROI==c("2","4")),]

data_group12$MMSE_TauCDR1 <- as.numeric(data_group12$MMSE_TauCDR1)
data_group12$EDUC <- as.numeric(data_group12$EDUC)

data_group14$MMSE_TauCDR1 <- as.numeric(data_group14$MMSE_TauCDR1)
data_group14$EDUC <- as.numeric(data_group14$EDUC)

data_group24$MMSE_TauCDR1 <- as.numeric(data_group24$MMSE_TauCDR1)
data_group24$EDUC <- as.numeric(data_group24$EDUC)


demo_group12 <- data_group12%>%
  select(3:14, -Tau_EC, -Tau_IT, -Tracer, -AmyPositivity, MMSE_TauCDR1, PACC_BaseCog1, Biomarker_group_metaROI, RACE, ETHNIC)%>%
  tbl_summary(by=Biomarker_group_metaROI,
              #type = list(MMSE_TauCDR1 ~ "continuous", EDUC ~ "continuous"),
              statistic=list(all_continuous()~ "{mean}({sd})",
                             all_categorical()~"{n}/{N}({p}%)"),
              digits = all_continuous()~2,
              label = c(AGE ~ "Age, years", SEX ~ "Sex,F:M", EDUC ~ "Education, years", APOE4 = "APOE ɛ4 carriers, n", 
                        ETHNIC ~ "Ethnicity", RACE ~ "Race", Centiloid = "Global Aβ SUVR/DVR", Tau_metaROI = "Temporal meta-ROI SUVR", 
                        PercentHIPP = "Hippocampal volume (% of TIV)", Cort_Thick = "Temporal cortical thickness",
                        MMSE_TauCDR1 = "MMSE (/30)", PACC_BaseCog1 = "PACC-5, baseline)"))%>%
  bold_labels()%>%
  italicize_levels()%>%
  add_p()%>%
  as_gt()%>%             
  gt::gtsave(filename = "./5_25_22/demo_group_12_N151_6-14-22.png")
demo_group12
data_group1<-data[which(data$Biomarker_group_metaROI==c("1")),]
data_group2<-data[which(data$Biomarker_group_metaROI==c("2")),]
data_group4<-data[which(data$Biomarker_group_metaROI==c("4")),]

mean(data_group1$MMSE_TauCDR1)
sd(data_group1$MMSE_TauCDR1)
mean(data_group1$EDUC)
sd(data_group1$EDUC)

mean(data_group2$MMSE_TauCDR1)
sd(data_group2$MMSE_TauCDR1)
mean(data_group2$EDUC)
sd(data_group2$EDUC)

mean(data_group4$MMSE_TauCDR1)
sd(data_group4$MMSE_TauCDR1)
mean(data_group4$EDUC)
sd(data_group4$EDUC)


demo_group14 <- data_group14%>%
  select(3:14, -Tau_EC, -Tau_IT, -Tracer, -AmyPositivity, MMSE_TauCDR1, PACC_BaseCog1, Biomarker_group_metaROI, RACE, ETHNIC)%>%
  tbl_summary(by=Biomarker_group_metaROI,
              #type = list(MMSE_TauCDR1 ~ "continuous", EDUC ~ "continuous"),
              statistic=list(all_continuous()~ "{mean}({sd})",
                             all_categorical()~"{n}/{N}({p}%)"),
              digits = all_continuous()~2,
              label = c(AGE ~ "Age, years", SEX ~ "Sex,F:M", EDUC ~ "Education, years", APOE4 = "APOE ɛ4 carriers, n", 
                        ETHNIC ~ "Ethnicity", RACE ~ "Race", Centiloid = "Global Aβ SUVR/DVR", Tau_metaROI = "Temporal meta-ROI SUVR", 
                        PercentHIPP = "Hippocampal volume (% of TIV)", Cort_Thick = "Temporal cortical thickness",
                        MMSE_TauCDR1 = "MMSE (/30)", PACC_BaseCog1 = "PACC-5, baseline)"))%>%
  bold_labels()%>%
  italicize_levels()%>%
  add_p()%>%
  as_gt()%>%             
  gt::gtsave(filename = "./5_25_22/demo_group_14_N151_6-14-22.png")
demo_group14

demo_group24 <- data_group24%>%
  select(3:14, -Tau_EC, -Tau_IT, -Tracer, -AmyPositivity, MMSE_TauCDR1, PACC_BaseCog1, Biomarker_group_metaROI, RACE, ETHNIC)%>%
  tbl_summary(by=Biomarker_group_metaROI,
              #type = list(MMSE_TauCDR1 ~ "continuous", EDUC ~ "continuous"),
              statistic=list(all_continuous()~ "{mean}({sd})",
                             all_categorical()~"{n}/{N}({p}%)"),
              digits = all_continuous()~2,
              label = c(AGE ~ "Age, years", SEX ~ "Sex,F:M", EDUC ~ "Education, years", APOE4 = "APOE ɛ4 carriers, n", 
                        ETHNIC ~ "Ethnicity", RACE ~ "Race", Centiloid = "Global Aβ SUVR/DVR", Tau_metaROI = "Temporal meta-ROI SUVR", 
                        PercentHIPP = "Hippocampal volume (% of TIV)", Cort_Thick = "Temporal cortical thickness",
                        MMSE_TauCDR1 = "MMSE (/30)", PACC_BaseCog1 = "PACC-5, baseline)"))%>%
  bold_labels()%>%
  italicize_levels()%>%
  add_p()%>%
  as_gt()%>%             
  gt::gtsave(filename = "./5_25_22/demo_group_24_N151_6-14-22.png")
demo_group24



###### Load data and select columns of interest
data<-data%>%
  mutate(SuppTable1 = case_when(
    Converter_after_tauPET == 0 & Biomarker_group_metaROI == 1 ~ "1 CU",
    Converter_after_tauPET == 0 & Biomarker_group_metaROI == 2 ~ "2 CU",
    Converter_after_tauPET == 0 & Biomarker_group_metaROI == 3 ~ "3 CU",
    Converter_after_tauPET == 0 & Biomarker_group_metaROI == 4 ~ "4 CU",
    Converter_after_tauPET == 1 & Biomarker_group_metaROI == 1 ~ "1 MCI",
    Converter_after_tauPET == 1 & Biomarker_group_metaROI == 2 ~ "2 MCI",
    Converter_after_tauPET == 1 & Biomarker_group_metaROI == 3 ~ "3 MCI",
    Converter_after_tauPET == 1 & Biomarker_group_metaROI == 4 ~ "4 MCI"))

###### Descriptive Statistics ######
theme_gtsummary_compact()

SuppTable1_demo <- data%>%
  select(3:14, -Tau_EC, -Tau_IT, -Tracer, -AmyPositivity, MMSE_TauCDR1, PACC_BaseCog1, ETHNIC, RACE, SuppTable1)%>%
  tbl_summary(by=SuppTable1,
              #type = list(MMSE_TauCDR1 ~ "continuous"),
    statistic=list(all_continuous()~ "{mean}({sd})",
                   all_categorical()~"{n}({p}%)"),
    digits = all_continuous()~2,
    label = c(AGE ~ "Age, years", SEX ~ "Sex,F:M", EDUC ~ "Education, years", 
              ETHNIC ~ "Ethnicity", RACE ~ "Race", APOE4 = "APOE ɛ4 carriers, n", 
              Centiloid = "Global Aβ SUVR/DVR", Tau_metaROI = "Temporal meta-ROI SUVR", 
              PercentHIPP = "Hippocampal volume (% of TIV)", Cort_Thick = "Temporal cortical thickness",
              MMSE_TauCDR1 = "MMSE (/30)", PACC_BaseCog1 = "PACC-5, baseline)"))%>%
  bold_labels()%>%
  italicize_levels()%>%
  as_gt()%>%                     
  gt::gtsave(filename = "./5_25_22/SuppTable1Demographics_BiomarkerSplit_6-14-22.png")
SuppTable1_demo
CU1<-data[which(data$SuppTable1==c("1 CU")),]
CU2<-data[which(data$SuppTable1==c("2 CU")),]
CU3<-data[which(data$SuppTable1==c("3 CU")),]
CU4<-data[which(data$SuppTable1==c("4 CU")),]

MCI1<-data[which(data$SuppTable1==c("1 MCI")),]
MCI2<-data[which(data$SuppTable1==c("2 MCI")),]
MCI3<-data[which(data$SuppTable1==c("3 MCI")),]
MCI4<-data[which(data$SuppTable1==c("4 MCI")),]

mean(CU1$MMSE_TauCDR1)
sd(CU1$MMSE_TauCDR1)

mean(CU2$MMSE_TauCDR1)
sd(CU2$MMSE_TauCDR1)

mean(CU3$MMSE_TauCDR1)
sd(CU3$MMSE_TauCDR1)

mean(CU4$MMSE_TauCDR1)
sd(CU4$MMSE_TauCDR1)

mean(MCI1$MMSE_TauCDR1)
sd(MCI1$MMSE_TauCDR1)

mean(MCI2$MMSE_TauCDR1)
sd(MCI2$MMSE_TauCDR1)

mean(MCI3$MMSE_TauCDR1)
sd(MCI3$MMSE_TauCDR1)

mean(MCI4$MMSE_TauCDR1)
sd(MCI4$MMSE_TauCDR1)














