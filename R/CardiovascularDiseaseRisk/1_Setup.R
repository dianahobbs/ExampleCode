# Title: ADRC Cohort Organization for Framingham Replication Study (Rabin)
# Cleaning and merging files  
# Author: Diana Hobbs
# Date: April 2022

#Functions to set up environment
rm(list=ls()) #removes all variables from the current environment
options(scipen = 999) #forces R not to use scientific notation

###### Set Up ######
# Load packages
pacman::p_load(tidyverse,jtools,naniar,lubridate,janitor)

###### Load data and select columns of interest ###### 

# demographics !*!*!!*!*! First open in excel .xlsx & format BIRTH date to 3/14/2012 and THEN save to .csv !*!*!!*!*!
demo <- read.csv("./Data/demographics.csv")%>%
  replace_with_na(replace=list(BIRTH = "###############################################################################################################################################################################################################################################################"))%>%
  mutate(BIRTH = as.Date(BIRTH, "%m/%d/%Y"))%>%
  # recode 1:Male, 2:Female
  mutate(GENDER = str_replace_all(GENDER, c("1"="Male", "2"="Female")))%>%
  rename(SEX=GENDER)%>%
  dplyr::select(ID,BIRTH,SEX,EDUC)%>%
  na.omit()

# apoe 
gene <- read.csv("./Data/apoe.csv")%>%
  mutate(APOE4 = str_replace_all(apoe,c("22"="0", "33"="0", "23"="0", "32"="0",
                                        "24"="1", "42"="1", # at least one E4 allele coded as E4+**
                                        "34"="1", "43"="1", "44"="1")))%>%
  dplyr::select(ID,APOE4)%>%
  na.omit()

# merge demographics with apoe genetics
demo<-merge(demo,gene,by="ID")

# cdr
cdr <- read.csv("./Data/mod_b4_cdr.csv")%>%
  mutate(cdr = case_when(cdr == 0 ~ "Cognitively Unimpaired",
                         cdr>0 ~ "Cognitively Impaired"))%>%
  rename(CDR=cdr)%>%
  dplyr::select(ID,TIME,CDR,MMSE)%>%
  na.omit()

# physical eval
## ****first open in excel .xlsx and format TESTDATE to 3/14/2012 and THEN save to csv
vital <- read.csv("./Data/mod_b1_physical_eval.csv")%>%
  replace_with_na(replace=list(HEIGHT=c(9999,999,99.9,99.0,88.8,100,183),WEIGHT=c(999,888),BPSYS = c(999,888)))%>%
  mutate(TESTDATE = as.Date(TESTDATE, "%m/%d/%Y"))%>%
  group_by(ID)%>%
  fill(HEIGHT, .direction = "down")%>%
  fill(HEIGHT, .direction = "up")%>%
  ungroup(ID)%>%
  mutate(BMI = (WEIGHT/(HEIGHT^2)*703))%>%
  dplyr::select(ID,TESTDATE,TIME,BMI,BPSYS)%>%
  na.omit()

# health history 
smoke <- read.csv("./Data/mod_a5_health_history.csv")%>%
  mutate(SMOKE = case_when(TOBAC30 == 0 ~ 0,
                           TOBAC30 == 1 ~ 1))%>%
  dplyr::select(ID,TIME,SMOKE)%>%
  na.omit()

# merge vitals & cdr & smoke
data <- merge(cdr, vital, by=c("ID","TIME"), all=F)
data <- merge(data, smoke, by=c("ID","TIME"), all=F)
# merge demo & data
data <- merge(demo,data, by="ID")
rm(cdr,demo,gene,smoke,vital)

# inclusion criteria - CDR baseline MUST be 0 and cap age at 75  
data <- data%>%
  mutate(AGE = as.numeric(difftime(TESTDATE, BIRTH, unit="weeks"))/52.25)%>%
  #filter(TIME==1 & CDR=="Cognitively Unimpaired" & AGE<=75 & MMSE>=27)%>%
  dplyr::select(-BIRTH)

###### PET ######

# amyloid (PiB AND av45) !*!*!!*!*! First open in excel .xlsx & format PET date to 3/14/2012 and THEN save to .csv !*!*!!*!*!
# PiB & Centiloid Conversion
PiB <- read.csv("./Data/PiB.csv")%>%
  mutate(amy_Date = as.Date(PET_Date, "%m/%d/%Y"))%>%
  mutate(PiB_Cent = 45.0 * pib_fsuvr_rsf_tot_cortmean - 47.5)%>%
  dplyr::select(ID,amy_Date,PiB_Cent)%>%
  na.omit()

# av45 & Centiloid Conversion
av45 <- read.csv("./Data/av45.csv")%>%
  mutate(amy_Date = as.Date(PET_Date, "%m/%d/%Y"))%>%
  mutate(av45_Cent = 53.6 * av45_fsuvr_rsf_tot_cortmean - 43.2)%>%
  dplyr::select(ID,amy_Date,av45_Cent)%>%
  na.omit()

# Merge PiB and av45 
amyloid<-merge(PiB,av45,by=c("ID","amy_Date"),all=T)%>%
  mutate(PiB_Cent = replace_na(PiB_Cent,9999), # code NAs as 0
         av45_Cent = replace_na(av45_Cent,9999))%>% 
  mutate(Tracer = case_when(
    PiB_Cent != "9999" ~ "PiB", 
    PiB_Cent == "9999" ~ "av45",
    av45_Cent == "9999" ~ "PiB"))%>%
  mutate(Centiloid = case_when(
    PiB_Cent != "9999" ~ PiB_Cent,
    PiB_Cent == "9999" ~ av45_Cent,
    av45_Cent == "9999" ~ PiB_Cent))%>%
  mutate(AmyPositivity = case_when(
    Tracer=="PiB" & Centiloid >= 16.4 ~ "A+",
    Tracer=="PiB" & Centiloid < 16.4 ~ "A-",
    Tracer=="av45" & Centiloid >= 20.6 ~ "A+",
    Tracer=="av45" & Centiloid < 20.6 ~ "A-"))%>%
  dplyr::select(ID,amy_Date,Tracer,Centiloid,AmyPositivity)

data <- merge(data, amyloid, by="ID", keep.all=T)%>%
  mutate(demoamy_DIFF = as.numeric(difftime(TESTDATE, amy_Date, unit="weeks"))/52.25)%>%
  mutate(demoamy_DIFF = abs(demoamy_DIFF))%>%
  filter(between(demoamy_DIFF, 0, 1))%>%
  group_by(ID)%>%
  arrange(amy_Date)%>%
  mutate(TIME=case_when((row_number()==1) ~ "1"))%>%
  ungroup()%>%
  filter(TIME==1)%>%
  dplyr::select(-demoamy_DIFF)
rm(amyloid,av45,PiB)

# Create Framingham score 
write.csv(data,"./Data/Rabin_Replication_Demographics.csv",row.names = FALSE)
rm(data)
source("1.1_Rabin_Replication_FramChart.R")
data <- data%>%
  mutate(CVD_Colsp = case_when(CVD_Risk == ">=30%" ~ ">20%",
                               CVD_Risk == "20-30%" ~ ">20%",
                               CVD_Risk == "10-20%" ~ "10-20%",
                               CVD_Risk == "5-10%" ~ "5-10%",
                               CVD_Risk == "<5%" ~ "<5%"))

# Create tau and MR files 
source("1.2_Rabin_Replication_Imaging.R")

# Only keep IDs in tau and MR that are present in demographics file
ID <- as.data.frame(table(data$ID))
tau <- tau[(tau$ID %in% ID$Var1),]
mri <- mri[(mri$ID %in% ID$Var1),]
PiB <- PiB[(PiB$ID %in% ID$Var1),]

# Only keep the first tau and MR scan in these 
# first merge tau & MR by date within 1 yr 
mriDate <- mri%>%dplyr::select(ID,MR_Date)
tauDate <- tau%>%dplyr::select(ID,tau_Date)
imageDate <- merge(mriDate,tauDate,by="ID")%>%
  mutate(mritau_DIFF = as.numeric(difftime(MR_Date, tau_Date, unit="weeks"))/52.25)%>%
  mutate(mritau_DIFF = abs(mritau_DIFF))%>%
  filter(between(mritau_DIFF, 0, 1))%>%
  dplyr::select(-c(mritau_DIFF)) #, TIME

mri <- merge(mri,imageDate,by=c("ID","MR_Date"),keep.all=F)
  
tau <- merge(tau,imageDate,by=c("ID","tau_Date"),keep.all=F)
rm(ID,imageDate,mriDate,tauDate)

# merge demographics with MR and tau 
mri <- merge(data,mri,by="ID",keep.all=F)
tau <- merge(data,tau,by="ID",keep.all=F)
PiB <- merge(data,PiB,by="ID",keep.all=F)
PiB <- PiB%>%filter(Tracer != "av45")%>%
  mutate(followup_time = 0)%>%
  replace_with_na(replace = list(followup_time=c("0")))%>%
  group_by(ID)%>%
  arrange(PET_Date)%>%
  mutate(TIME=case_when((row_number()==1) ~ "1"))%>%
  ungroup()%>%
  filter(TIME==1 & CDR=="Cognitively Unimpaired" & AGE<=75 & MMSE>=27)

# followup_time = Time between baseline AB and tau/MR
tau <- tau%>%
  mutate(followup_time = as.numeric(abs(difftime(tau_Date,amy_Date,unit="weeks"))/52.25))%>%
  filter(between(followup_time, 0, 2))%>%
  group_by(ID)%>%
  arrange(tau_Date)%>%
  mutate(TIME=case_when((row_number()==1) ~ "1"))%>%
  ungroup()%>%
  filter(TIME==1 & CDR=="Cognitively Unimpaired" & AGE<=75 & MMSE>=27)
tauID<-as.data.frame(table(tau$ID))

mri <- mri%>%
  mutate(followup_time = as.numeric(abs(difftime(MR_Date,amy_Date,unit="weeks"))/52.25))%>%
  filter(between(followup_time, 0, 2))%>%
  group_by(ID)%>%
  arrange(MR_Date)%>%
  mutate(TIME=case_when((row_number()==1) ~ "1"))%>%
  ungroup()%>%
  filter(TIME==1 & CDR=="Cognitively Unimpaired" & AGE<=75 & MMSE>=27)%>%
  dplyr::select(-c(tau_Date,TIME))
# Keep only IDs that are also in tau file 
mri <- mri[(mri$ID %in% tauID$Var1),]
rm(data,tauID)

###### Exclusion Criteria: No history stroke, no active alcoholism, no major medical/psychiatric/neurological conditions ######

# diagnoses
dx <- read.csv("./Data/mod_b4_cdr.csv")%>%
  dplyr::select(ID,TIME,dx1,dx2)%>%
  replace_with_na(replace = list(dx1=c("Q","."),
                                 dx2=c("A","B",".","0")))%>%
  replace_na(replace = list(dx1="NONE",dx2="NONE"))%>%
  mutate(dx1 = str_replace_all(dx1, c("Vascular Demt, primary"="nonAD Dementia",
                                      "Vascular Demt, secondary"="nonAD Dementia",
                                      "Non AD dem, Other secondary"="nonAD Dementia",
                                      "Non AD dem, Other primary"="nonAD Dementia",
                                      "Incipient Non-AD dem"="nonAD Dementia",
                                      "DLBD, primary"="nonAD Dementia",
                                      "DLBD, secondary"="nonAD Dementia",
                                      "Frontotemporal demt. secn"="nonAD Dementia",
                                      "Frontotemporal demt. prim"="nonAD Dementia",
                                      "Dementia/PD, primary"="nonAD Dementia",
                                      "Dementia/PD, secondary"="nonAD Dementia",
                                      "Atypical AD dem/SAD dem"="AD Dementia",
                                      "AD dem/ProA Prior to AD dem"="AD Dementia",
                                      "AD dem/PCD prior to AD dem"="AD Dementia",
                                      "AD dem/PCD after AD dem"="AD Dementia",
                                      "AD dem/Other prior to AD dem"="AD Dementia",
                                      "AD dem/FLD prior to AD dem"="AD Dementia",
                                      "AD dem/FLD after AD dem"="AD Dementia",
                                      "AD dem w/ProAph w/subsequent demt"="AD Dementia",
                                      "AD dem w/ProAph w/demt at onset"="AD Dementia",
                                      "AD dem w/post cort dys/demt at onset"="AD Dementia",
                                      "AD dem w/post cort dys"="AD Dementia",
                                      "AD dem w/PDI precede AD dem contribut"="AD Dementia",
                                      "AD dem w/PDI after AD dem not contrib"="AD Dementia",
                                      "AD dem w/PDI after AD dem contribut"="AD Dementia",
                                      "AD dem w/oth unusual features/demt on"="AD Dementia",
                                      "AD dem w/oth unusual features"="AD Dementia",
                                      "AD dem w/oth unusual feat/subs demt"="AD Dementia",
                                      "AD dem w/oth \\(\\list B\\)\\ not contrib"="AD Dementia",
                                      "AD dem w/oth \\(\\list B\\)\\ contribut"="AD Dementia",
                                      "AD dem w/Hallucinosis"="AD Dementia",
                                      "AD dem w/dispropor language dys"="AD Dementia",
                                      "AD dem w/DIP after AD dem not contrib"="AD Dementia",
                                      "AD dem w/depresss, not contribut"="AD Dementia",
                                      "AD dem w/depresss, contribut"="AD Dementia",
                                      "AD dem w/CVD not contrib"="AD Dementia",
                                      "AD dem w/CVD contribut"="AD Dementia",
                                      "AD dem w/abnorm behav syndrome"="AD Dementia",
                                      "AD dem w/abnorm behav syn/demt at ons"="AD Dementia",
                                      "AD dem visuospatial, with"="AD Dementia",
                                      "AD dem visuospatial, prior"="AD Dementia",
                                      "AD dem visuospatial, after"="AD Dementia",
                                      "AD dem Language dysf with"="AD Dementia",
                                      "AD dem Language dysf prior"="AD Dementia",
                                      "AD dem Language dysf after"="AD Dementia",
                                      "AD dem distrubed social, with"="AD Dementia",
                                      "AD dem distrubed social, prior"="AD Dementia",
                                      "AD dem distrubed social, after"="AD Dementia",
                                      "AD dem cannot be primary"="AD Dementia",
                                      "AD dem w/Frontal lobe/demt at onset"="AD Dementia",
                                      "AD dem w/DIP after AD dem contribut"="AD Dementia",
                                      "uncertain, possible NON AD dem"="Uncertain Dementia",
                                      "Unc: ques. Impairment"="Uncertain Impairment",
                                      "Unc: impair resolved"="Uncertain Impairment",
                                      "Unc: impair reversible"="Uncertain Impairment",
                                      "0.5 in memory only"="Uncertain Impairment",
                                      "Incipient demt PTP"="nonAD Dementia",
                                      "CDR 0.5 PTP"="nonAD Dementia",
                                      "Frontotemporal Dementia"="nonAD Dementia",
                                      "Parkinson Dementia"="nonAD Dementia",
                                      "Vascular Demt"="nonAD Dementia",
                                      "PTP"="nonAD Dementia",
                                      "Diffuse Lewy Body Disease"="nonAD Dementia",
                                      "uncertain impairment"="Uncertain Dementia",
                                      "uncertain dementia"="Uncertain Dementia",
                                      "Cognitively normal"="Cognitively Normal",
                                      "non AD dementia"="nonAD Dementia")))%>%
  mutate(dx2 = str_replace_all(dx2, c("Active Affective disorder"="Exclude",
                                      "Active Alcoholism"="Exclude",
                                      "Active B-12 Deficiency"="Exclude",
                                      "Active Bereavement"="Exclude",
                                      "Active CVD"="Exclude",
                                      "Active Global Cerebral Hypop"="Exclude",
                                      "Active Hydrocephalus"="Exclude",
                                      "Active Hypothyroidism"="Exclude",
                                      "Active Major Head Trauma"="Exclude",
                                      "Active Med-induced Cog Dys"="Exclude",
                                      "Active Mood disorder"="Exclude",
                                      "Active Other neurol/med diagnoses"="Exclude",
                                      "Active PD drug-induced"="Exclude",
                                      "Active PD idiopathic"="Exclude",
                                      "Active Seizure Disorder"="Exclude",
                                      "AD dem w/CVD not contrib"="Exclude",
                                      "AD dem w/depresss, not contribut"="Exclude",
                                      "affective disorder"="Exclude",
                                      "cerebrovascular dis"="Exclude",
                                      "Frontotemporal demt. prim"="Exclude",
                                      "Hypothyroidism"="Exclude",
                                      "med induced cog. dys"="Exclude",
                                      "other neurol/med diagnoses"="Exclude",
                                      "Park Idiopathic"="Exclude",
                                      "ProAph w/o dement"="Exclude",
                                      "Remote Affective disorder"="Remote",
                                      "Remote Alcoholism"="Remote",
                                      "Remote Amnestic syndrome"="Remote",
                                      "Remote B-12 Deficiency"="Remote",
                                      "Remote Bereavement"="Remote",
                                      "Remote CVD"="Remote",
                                      "Remote Global Cerebral Hypop"="Remote",
                                      "Remote Hydrocephalus"="Remote",
                                      "Remote Hypothyroidism"="Remote",
                                      "Remote Major Head Trauma"="Remote",
                                      "Remote Med-induced Cog Dys"="Remote",
                                      "Remote Mood disorder"="Remote",
                                      "Remote Other neurol/med diagnoses"="Remote",
                                      "Remote PD drug-induced"="Remote",
                                      "Remote PD idiopathic"="Remote",
                                      "Remote Seizure Disorder"="Remote",
                                      "Remote Exclude"="Remote")))%>%
  filter(dx1 != "nonAD Dementia", dx1 != "Uncertain Impairment", dx1 != "Uncertain Dementia", dx1 != "AD Dementia")%>%
  filter(dx2 != "Exclude")

stroke <- read.csv("./Data/mod_d1_diagnoses.csv")%>%
  dplyr::select(ID,TIME,STROKE)%>%
  replace_na(replace = list(STROKE="0"))%>%
  filter(STROKE != "1")

Include <- merge(dx,stroke,by=c("ID", "TIME"))%>%
  filter(TIME==1)
IncludeIDs <- as.data.frame(table(Include$ID))

# Keep only IDs that are also in inclusion file 
tau <- tau[(tau$ID %in% IncludeIDs$Var1),]
mri <- mri[(mri$ID %in% IncludeIDs$Var1),]
PiB <- PiB[(PiB$ID %in% IncludeIDs$Var1),]
rm(dx,Include,IncludeIDs,stroke)

write.csv(mri,"./Data/Rabin_Replication_MRI.csv",row.names = FALSE)
write.csv(tau,"./Data/Rabin_Replication_TAU.csv",row.names = FALSE)
write.csv(PiB,"./Data/Rabin_Replication_PIB.csv",row.names = FALSE)

# normalize the MR data
source("1.3_Rabin_Replication_MR-Normalize.R")









