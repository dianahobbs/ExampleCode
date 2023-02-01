# Title: ADRC - Demographics Rubin Replication Study
# Author: Diana Hobbs
# Date: April 2022

#### Set Up ####
# Load packages 
pacman::p_load(tidyverse,jtools,naniar,lubridate,lme4,sfsmisc,foreign,MASS,scales,ggseg,cowplot,gtsummary,ggplot2,RColorBrewer,gt)

###### LOAD & SETUP DATA 
mri <- read.csv("./Data/Rabin_Replication_MRI_Normalized.csv")%>%
  mutate(CVD = str_replace_all(CVD_Colsp, c("<5%"="0", "5-10%"="1", "10-20%"="2",">20%"="3")))%>%
  select(AGE,SEX,BMI,BPSYS,SMOKE,EDUC,MMSE,APOE4,Centiloid,CVD,followup_time,caudate,caudal.middle.frontal)

tau <- read.csv("./Data/Rabin_Replication_TAU.csv")%>%
  mutate(CVD = str_replace_all(CVD_Colsp, c("<5%"="0", "5-10%"="1", "10-20%"="2",">20%"="3")))%>%
  select(AGE,SEX,BMI,BPSYS,SMOKE,EDUC,MMSE,APOE4,Centiloid,CVD,followup_time,caudate,caudal.middle.frontal)

PiB <- read.csv("./Data/Rabin_Replication_PIB.csv")%>%
  mutate(CVD = str_replace_all(CVD_Colsp, c("<5%"="0", "5-10%"="1", "10-20%"="2",">20%"="3")))%>%
  select(AGE,SEX,BMI,BPSYS,SMOKE,EDUC,MMSE,APOE4,Centiloid,CVD,followup_time,caudate,caudal.middle.frontal)

# Reassign variable types
mri$SEX <- as.factor(mri$SEX); tau$SEX <- as.factor(tau$SEX); PiB$SEX <- as.factor(PiB$SEX)
mri$APOE4 <- as.factor(mri$APOE4); tau$APOE4 <- as.factor(tau$APOE4); PiB$APOE4 <- as.factor(PiB$APOE4)
mri$CDR <- as.factor(mri$CDR); tau$CDR <- as.factor(tau$CDR); PiB$CDR <- as.factor(PiB$CDR)
mri$SMOKE <- as.factor(mri$SMOKE); tau$SMOKE <- as.factor(tau$SMOKE); PiB$SMOKE <- as.factor(PiB$SMOKE)

mri$CVD <- as.factor(mri$CVD)%>%
  recode("0"="<5%","1"="5-10%","2"="10-20%","3"=">20%")
mri$SMOKE <- as.factor(mri$SMOKE)%>%
  recode("0"="Non-Smoker", "1"="Smoker")
mri$APOE4 <- as.factor(mri$APOE4)%>%
  recode("0"="Non-Carrier", "1"="Carrier")

tau$CVD <- as.factor(tau$CVD)%>%
  recode("0"="<5%","1"="5-10%","2"="10-20%","3"=">20%")
tau$SMOKE <- as.factor(tau$SMOKE)%>%
  recode("0"="Non-Smoker", "1"="Smoker")
tau$APOE4 <- as.factor(tau$APOE4)%>%
  recode("0"="Non-Carrier", "1"="Carrier")

PiB$CVD <- as.factor(PiB$CVD)%>%
  recode("0"="<5%","1"="5-10%","2"="10-20%","3"=">20%")
PiB$SMOKE <- as.factor(PiB$SMOKE)%>%
  recode("0"="Non-Smoker", "1"="Smoker")
PiB$APOE4 <- as.factor(PiB$APOE4)%>%
  recode("0"="Non-Carrier", "1"="Carrier")

###### Demographics Table (Same in MR and Tau) ######
theme_gtsummary_compact()
demo <- mri%>%
  select(-c(followup_time,caudate,caudal.middle.frontal))%>%
  tbl_summary(type = list(MMSE ~ "continuous",
                          EDUC ~ "continuous"),
              statistic=list(all_continuous()~ "{mean}({sd})",
                             all_categorical()~"{n}({p}%)"),
              digits = all_continuous()~2,
              label = c(AGE ~ "Age (yrs)", SEX ~ "Sex", BMI = "Body Mass Index (kg/m^2)", BPSYS = "Systolic Blood Pressure", SMOKE = "Current Smoking Status", 
                        EDUC ~ "Education (yrs)", MMSE = "Mini-Mental State Examination",APOE4 = "APOE 𝜀4  Status", Centiloid = "Centiloid", 
                        CVD = "Cardiovascular Risk"))%>%
  bold_labels()%>%
  italicize_levels()
demo

demo_PiB <- PiB%>%
  select(-c(followup_time,caudate,caudal.middle.frontal))%>%
  tbl_summary(type = list(MMSE ~ "continuous",
                          EDUC ~ "continuous"),
              statistic=list(all_continuous()~ "{mean}({sd})",
                             all_categorical()~"{n}({p}%)"),
              digits = all_continuous()~2,
              label = c(AGE ~ "Age (yrs)", SEX ~ "Sex", BMI = "Body Mass Index (kg/m^2)", BPSYS = "Systolic Blood Pressure", SMOKE = "Current Smoking Status", 
                        EDUC ~ "Education (yrs)", MMSE = "Mini-Mental State Examination",APOE4 = "APOE 𝜀4  Status", Centiloid = "Centiloid", 
                        CVD = "Cardiovascular Risk"))%>%
  bold_labels()%>%
  italicize_levels()
demo_PiB

tbl_merge(tbls=list(demo, demo_PiB), tab_spanner=c("**MRI and Tau-PET**", "**PiB-PET**"))%>%
  as_gt() %>%
  gt::tab_header(title = "Table 1. Demographics")
tbl_merge

tiff('./Rabin_Replication_Output/Demographics.tiff', units="in", width=7, height=9, res=200,type="cairo")
dev.off()

# Followup time from baseline amyloid scan 
followMR <- mri%>%
  select(followup_time)%>%
  tbl_summary(statistic=list(all_continuous()~ "{mean}({sd})"),
              digits = all_continuous()~2,
              label = c(followup_time = "Years from Baseline beta-Amyloid PET"))%>%
  bold_labels()%>%
  italicize_levels()
followMR

followTau <- tau%>%
  select(followup_time)%>%
  tbl_summary(statistic=list(all_continuous()~ "{mean}({sd})"),
              digits = all_continuous()~2,
              label = c(followup_time = "Years from Baseline beta-Amyloid PET"))%>%
  bold_labels()%>%
  italicize_levels()
followTau

tbl_merge(tbls=list(followMR, followTau), tab_spanner=c("**MRI**", "**Tau-PET**"))%>%
  as_gt() %>%
  gt::tab_header(title = "Table XX. Follow-up Time Between Baseline beta-Amyloid-PET and MRI or Tau-PET")
tbl_merge
tiff('./Rabin_Replication_Output/Followup_time.tiff', units="in", width=7, height=9, res=200,type="cairo")
dev.off()

###### regions that survived robust threshold
CVD_Tau1 <- ggplot(data = tau, aes(x=CVD, y=caudal.middle.frontal))+
  geom_violin(trim = FALSE , alpha = 0.5, aes(fill = CVD), size = .5)+
  geom_boxplot(width = .10, size = 0.5, outlier.shape = NA, aes(fill = CVD), alpha = 0.5)+
  scale_fill_manual(values=c("#74BAE3","#982B87","#261175","#BC4E63"), guide="none")+
  labs(y = "Caudal Middle Frontal SUVR", x = "Cardiovascular Disease Risk", title = "Cardiovascular Disease Risk and Tau-PET Regional SUVR")+
  theme(plot.title = element_text(hjust = 0.5), panel.background = element_blank(),
        axis.text = element_text(size=14), axis.title = element_text(size = 14), title = element_text(size = 14))
tiff('./Rabin_Replication_Output/Tau/Cort_violin.tiff', units="in", width=7, height=9, res=200,type="cairo")
CVD_Tau1
dev.off()

CVD_Tau2 <- ggplot(data = tau, aes(x=CVD, y=caudate))+
  geom_violin(trim = FALSE , alpha = 0.5, aes(fill = CVD), size = .5)+
  geom_boxplot(width = .10, size = 0.5, outlier.shape = NA, aes(fill = CVD), alpha = 0.5)+
  scale_fill_manual(values=c("#74BAE3","#982B87","#261175","#BC4E63"), guide="none")+
  labs(y = "Caudate", x = "Cardiovascular Disease Risk")+
  theme_classic()
tiff('./Rabin_Replication_Output/Tau/Subcort_violin.tiff', units="in", width=7, height=9, res=200,type="cairo")
CVD_Tau2
dev.off()

