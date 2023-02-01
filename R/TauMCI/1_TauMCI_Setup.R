# Title: ADRC Cohort Organization for Tau/MCI Replication Study (Sylvia/Cherie)
# Testing new way of redefining groups of neurodegeneration and cognitive decline  
# Author: Diana Hobbs
# Date: May 13, 2022

###### Set Up ######
# Load packages
pacman::p_load(tidyverse,jtools,naniar,lubridate,janitor,lme4,lmerTest)

###### Load data and select columns of interest ###### 
df <- read.csv("./Data/N151_Wide_5-25-22.csv")
# demographics !*!*!!*!*! First open in excel .xlsx & format BIRTH date to 3/14/2012 and THEN save to .csv !*!*!!*!*!
demo <- read.csv("./Data/demographics.csv")%>%
  replace_with_na(replace=list(BIRTH = "###############################################################################################################################################################################################################################################################"))%>%
  mutate(BIRTH = as.Date(BIRTH, "%m/%d/%Y"))%>%
  # recode 1:Male, 2:Female
  mutate(SEX = str_replace_all(GENDER, c("1"="Male", "2"="Female")))%>%
  mutate(RACE = str_replace_all(race,c("AIAN"="Asian","ASIAN"="Asian")))%>% #"ASIAN"="non-White","Black"="non-White","NHPI"="non-White"
  dplyr::select(ID,BIRTH,SEX,RACE,ETHNIC,EDUC)

# apoe 
gene <- read.csv("./Data/apoe.csv")%>%
  mutate(APOE4 = str_replace_all(apoe,c("22"="0", "33"="0", "23"="0", "32"="0",
                                        "24"="1", "42"="1", # double check that we want 2/4's to be coded as APOE4 carrier**
                                        "34"="1", "43"="1", "44"="1")))%>%
  dplyr::select(ID,APOE4)

# merge demographics with apoe genetics
demo<-merge(demo,gene,by="ID")

# tau !*!*!!*!*! First open in excel .xlsx & format PET date to 3/14/2012 and THEN save to .csv !*!*!!*!*!
# Create metaROI and VisualRead - do *NOT* use PVC 
# metaROI = mean(bilateral entorhinal, amygdala, fusiform, inferior temporal, middle temporal, and parahippocampal)
# VisualRead = mean(frontal, parietal, and occipital regions as defined in freesurfer https://surfer.nmr.mgh.harvard.edu/fswiki/CorticalParcellation)
tau<-read.csv("./Data/tau.csv")%>%
  mutate(Tau_Date = as.Date(PET_Date, "%m/%d/%Y"))%>%
  mutate(Tau_metaROI=(av1451_fsuvr_tot_ctx_entorhinal+av1451_fsuvr_tot_amygdala+av1451_fsuvr_tot_ctx_fusiform+av1451_fsuvr_tot_ctx_infertmp+av1451_fsuvr_tot_ctx_midtmp+av1451_fsuvr_tot_ctx_parahpcmpl)/6)%>%
  mutate(VisualRead = ((av1451_fsuvr_tot_ctx_superfrn + av1451_fsuvr_tot_ctx_rosmidfrn + av1451_fsuvr_tot_ctx_caudmidfrn + av1451_fsuvr_tot_ctx_parsopclrs + 
                          av1451_fsuvr_tot_ctx_parsorbls + av1451_fsuvr_tot_ctx_parstrngls + av1451_fsuvr_tot_ctx_latorbfrn + av1451_fsuvr_tot_ctx_medorbfrn + 
                          av1451_fsuvr_tot_ctx_precntrl + av1451_fsuvr_tot_ctx_paracntrl + av1451_fsuvr_tot_ctx_frnpole + av1451_fsuvr_tot_ctx_superprtl +
                          av1451_fsuvr_tot_ctx_inferprtl + av1451_fsuvr_tot_ctx_supramrgnl + av1451_fsuvr_tot_ctx_postcntrl + av1451_fsuvr_tot_ctx_precuneus +
                          av1451_fsuvr_tot_ctx_latocc + av1451_fsuvr_tot_ctx_lingual + av1451_fsuvr_tot_ctx_cuneus + av1451_fsuvr_tot_ctx_periclcrn +
                          av1451_fsuvr_tot_ctx_rosantcng + av1451_fsuvr_tot_ctx_caudantcng + av1451_fsuvr_tot_ctx_postcng + av1451_fsuvr_tot_ctx_isthmuscng)/24))%>%
  rename(Tau_EC=av1451_fsuvr_tot_ctx_entorhinal, Tau_IT=av1451_fsuvr_tot_ctx_infertmp)%>%
  dplyr::select(ID,Tau_Date,Tau_metaROI,Tau_EC,Tau_IT,VisualRead)%>%
  group_by(ID)%>% # Grouping by ID & arranging by tau date to find the first tau visit for each participant
  arrange(Tau_Date)%>%
  mutate(VISIT = case_when(
    (row_number()==1) ~ "Time1",
    (row_number()==2) ~ "Time2",
    (row_number()==3) ~ "Time3"))%>%
  ungroup()%>%
  filter(VISIT=="Time1")%>% # Select only the first tau visit
  na.omit() 

# Merge tau with demographics and select only those >60 yrs 
tau<-merge(demo,tau,by="ID",all=F)#%>%
  #mutate(AGE = as.numeric(difftime(Tau_Date, BIRTH, unit="weeks"))/52.25)%>%
  #filter(AGE >= 60)
# write.csv(tau, "./Data/tau_for_visualread.csv", row.names = FALSE)

###### MRI ###### !*!*!!*!*! First open in excel .xlsx & format MR date to 3/14/2012 and THEN save to .csv !*!*!!*!*!
# create cortical thickness measure as mean(bilateral entorhinal, fusiform, inferior temporal, middle temporal)
mri <- read.csv("./Data/mri_t3.csv", header=TRUE)%>%
  mutate(MR_Date = as.Date(MR_Date, "%m/%d/%Y"))%>%
  mutate(HIPP=MR_LV_HIPPOCAMPUS+MR_RV_HIPPOCAMPUS)%>%
  mutate(Cort_Thick = ((((MR_LT_ENTORHINAL + MR_RT_ENTORHINAL) +
                           (MR_LT_FUSIFORM + MR_RT_FUSIFORM) +
                           (MR_LT_INFRTMP + MR_RT_INFRTMP) +
                           (MR_LT_MIDTMP + MR_RT_MIDTMP)) / 8)))%>%
  rename(ICV=MR_TOTV_INTRACRANIAL)%>%
  mutate(PercentHIPP = ((HIPP * 100)/ICV))%>% # Hipp volume = sum as percentage of ICV
  dplyr::select(ID,MR_Date,Scanner,PercentHIPP,Cort_Thick)

# Merge tau/amyloid with MR
TauMR<-merge(tau,mri,by="ID")%>%
  mutate(AGE = as.numeric(difftime(MR_Date, BIRTH, unit="weeks"))/52.25)%>%
  filter(AGE >= 60)
TauMR<-TauMR%>%
  mutate(TauMR_DIFF = as.numeric(difftime(MR_Date, Tau_Date, unit="weeks"))/52.25)%>%
  mutate(TauMR_DIFF = abs(TauMR_DIFF))%>%
  filter(between(TauMR_DIFF, 0, 1.5))%>% # scans must be within 1.5 yrs of one another
  group_by(ID)%>% # Grouping by ID & arranging by tau date to find the first tau/amyloid/MR visit for each participant
  arrange(Tau_Date)%>%
  filter(row_number()==1)%>% # Select only the first tau visit
  ungroup()%>%
  dplyr::select(ID,AGE,SEX,RACE,ETHNIC,EDUC,APOE4,Tau_metaROI,Tau_EC,Tau_IT,PercentHIPP,Cort_Thick,Tau_Date,MR_Date,Scanner,TauMR_DIFF,VisualRead)
unique<-data.frame(table(TauMR$ID)) # all unique ID values, n = 440


###### Clinical and Cognitive ###### 

# CDR !*!*!!*!*! First open in excel .xlsx & format test date to 3/14/2012 and THEN save to .csv !*!*!!*!*!
CDR <- read.csv("./Data/mod_b4_cdr.csv")%>%
  mutate(testdate = as.Date(testdate, "%m/%d/%Y"))%>%
  rename(CDR_Date=testdate,CDR=cdr)%>%
  dplyr::select(ID,CDR,CDR_Date,MMSE,sumbox)

# Cognitive composite - PACC - !*!*!!*!*! First open in excel .xlsx & format test date to 3/14/2012 and THEN save to .csv !*!*!!*!*!
cog<-read.csv("./Data/psychometrics.csv")%>%
  mutate(psy_date = as.Date(psy_date, "%m/%d/%Y"))%>%
  dplyr::select(ID,psy_date,srtfree,lettnum,ANIMALS,tmb)%>%
  na.omit()

# Merge CDR with Cog 
CDRCog<-merge(CDR,cog,by="ID")%>%
  mutate(CDRCog_DIFF = as.numeric(difftime(psy_date, CDR_Date, unit="weeks"))/52.25)%>%
  mutate(CDRCog_DIFF = abs(CDRCog_DIFF))%>%
  filter(between(CDRCog_DIFF, 0, 1.5))%>% # testing dates must be within 1.5 yrs of one another
  group_by(ID)%>% # Grouping by ID & keeping only distinct CDR visits for each participant
  distinct(CDR_Date, .keep_all = T)%>%
  ungroup()

###### Revision Analyses##### 
# Merge CDRCog with MRI to do revision analyses
CDRCog_MR<-merge(CDRCog,mri,by="ID")%>%
  mutate(CDRCogMR_DIFF = as.numeric(difftime(MR_Date, CDR_Date, unit="weeks"))/52.25)%>%
  mutate(CDRCogMR_DIFF = abs(CDRCogMR_DIFF))%>%
  filter(between(CDRCogMR_DIFF,0,1.5))%>%
  group_by(ID)%>% 
  distinct(CDR_Date, .keep_all = T)%>%
  arrange(CDR_Date)%>%
  mutate(Visit= case_when(
    (row_number()==1) ~ "Cog1",
    (row_number()==2) ~ "Cog2",
    (row_number()==3) ~ "Cog3",
    (row_number()==4) ~ "Cog4",
    (row_number()==5) ~ "Cog5",
    (row_number()==6) ~ "Cog6",
    (row_number()==7) ~ "Cog7",
    (row_number()==8) ~ "Cog8",
    (row_number()==9) ~ "Cog9",
    (row_number()==10) ~ "Cog10",
    (row_number()==11) ~ "Cog11",
    (row_number()==12) ~ "Cog12",
    (row_number()==13) ~ "Cog13"))%>%
  ungroup()%>%
  dplyr::select(ID,Visit,psy_date,CDR,srtfree,lettnum,ANIMALS,tmb,PercentHIPP,Cort_Thick)
CDRCog_MR<-CDRCog_MR[!(CDRCog_MR$Visit=="Cog1" & CDRCog_MR$CDR>0),] # only keep people with baseline CDR 0

CDRCog_MR_Wide<-CDRCog_MR%>%
  pivot_wider(id_cols = 1,
              names_from = Visit,
              names_prefix = "",
              values_from = 3:10)%>%
  drop_na(CDR_Cog1)

data1<-CDRCog_MR_Wide # Determine who is a converter vs non-converter (Go back to the variables we used with MCI conversion setup)
data1$CDR_Cog1[is.na(data1$CDR_Cog1)] <- "NA"
data1$CDR_Cog2[is.na(data1$CDR_Cog2)] <- "NA"
data1$CDR_Cog3[is.na(data1$CDR_Cog3)] <- "NA"
data1$CDR_Cog4[is.na(data1$CDR_Cog4)] <- "NA"
data1$CDR_Cog5[is.na(data1$CDR_Cog5)] <- "NA"
data1$CDR_Cog6[is.na(data1$CDR_Cog6)] <- "NA"
data1$CDR_Cog7[is.na(data1$CDR_Cog7)] <- "NA"
data1$CDR_Cog8[is.na(data1$CDR_Cog8)] <- "NA"
data1$CDR_Cog9[is.na(data1$CDR_Cog9)] <- "NA"
data1$CDR_Cog10[is.na(data1$CDR_Cog10)] <- "NA"
data1$CDR_Cog11[is.na(data1$CDR_Cog11)] <- "NA"
data1$CDR_Cog12[is.na(data1$CDR_Cog12)] <- "NA"
data1$CDR_Cog13[is.na(data1$CDR_Cog13)] <- "NA"
# Converter_after_tauPET variable: 1: converter to MCI, 0: didn't convert
data1<-data1%>%
  mutate(Converter = case_when(
        CDR_Cog1==0 & CDR_Cog2==0 & CDR_Cog3==0 & CDR_Cog4==0 & CDR_Cog5==0 & CDR_Cog6==0 & CDR_Cog7==0 & CDR_Cog8==0 & CDR_Cog9==0 & CDR_Cog10==0 & CDR_Cog11==0 & CDR_Cog12==0 & CDR_Cog13==0~ "0", # if CDR at all 13 time-points is 0, then non-converter 
        CDR_Cog1==0 & CDR_Cog2==0 & CDR_Cog3==0 & CDR_Cog4==0 & CDR_Cog5==0 & CDR_Cog6==0 & CDR_Cog7==0 & CDR_Cog8==0 & CDR_Cog9==0 & CDR_Cog10==0 & CDR_Cog11==0 & CDR_Cog12==0 & CDR_Cog13=="NA"~ "0", # if CDR at all 13 time-points is 0, then non-converter 
        CDR_Cog1==0 & CDR_Cog2==0 & CDR_Cog3==0 & CDR_Cog4==0 & CDR_Cog5==0 & CDR_Cog6==0 & CDR_Cog7==0 & CDR_Cog8==0 & CDR_Cog9==0 & CDR_Cog10==0 & CDR_Cog11==0 & CDR_Cog12=="NA" & CDR_Cog13=="NA"~ "0", # if CDR at all 13 time-points is 0, then non-converter 
        CDR_Cog1==0 & CDR_Cog2==0 & CDR_Cog3==0 & CDR_Cog4==0 & CDR_Cog5==0 & CDR_Cog6==0 & CDR_Cog7==0 & CDR_Cog8==0 & CDR_Cog9==0 & CDR_Cog10==0 & CDR_Cog11=="NA" & CDR_Cog12=="NA" & CDR_Cog13=="NA"~ "0", # if CDR at all 13 time-points is 0, then non-converter 
        CDR_Cog1==0 & CDR_Cog2==0 & CDR_Cog3==0 & CDR_Cog4==0 & CDR_Cog5==0 & CDR_Cog6==0 & CDR_Cog7==0 & CDR_Cog8==0 & CDR_Cog9==0 & CDR_Cog10=="NA" & CDR_Cog11=="NA" & CDR_Cog12=="NA" & CDR_Cog13=="NA"~ "0", # if CDR at all 13 time-points is 0, then non-converter 
        CDR_Cog1==0 & CDR_Cog2==0 & CDR_Cog3==0 & CDR_Cog4==0 & CDR_Cog5==0 & CDR_Cog6==0 & CDR_Cog7==0 & CDR_Cog8==0 & CDR_Cog9=="NA" & CDR_Cog10=="NA" & CDR_Cog11=="NA" & CDR_Cog12=="NA" & CDR_Cog13=="NA"~ "0", # if CDR at all 13 time-points is 0, then non-converter 
        CDR_Cog1==0 & CDR_Cog2==0 & CDR_Cog3==0 & CDR_Cog4==0 & CDR_Cog5==0 & CDR_Cog6==0 & CDR_Cog7==0 & CDR_Cog8=="NA" & CDR_Cog9=="NA" & CDR_Cog10=="NA" & CDR_Cog11=="NA" & CDR_Cog12=="NA" & CDR_Cog13=="NA"~ "0", # if CDR at all 13 time-points is 0, then non-converter 
        CDR_Cog1==0 & CDR_Cog2==0 & CDR_Cog3==0 & CDR_Cog4==0 & CDR_Cog5==0 & CDR_Cog6=="NA" & CDR_Cog7=="NA" & CDR_Cog8=="NA" & CDR_Cog9=="NA" & CDR_Cog10=="NA" & CDR_Cog11=="NA" & CDR_Cog12=="NA" & CDR_Cog13=="NA"~ "0", # "" etc ""
        CDR_Cog1==0 & CDR_Cog2==0 & CDR_Cog3==0 & CDR_Cog4==0 & CDR_Cog5=="NA" & CDR_Cog6=="NA" & CDR_Cog7=="NA" & CDR_Cog8=="NA" & CDR_Cog9=="NA" & CDR_Cog10=="NA" & CDR_Cog11=="NA" & CDR_Cog12=="NA" & CDR_Cog13=="NA"~ "0",
        CDR_Cog1==0 & CDR_Cog2==0 & CDR_Cog3==0 & CDR_Cog4=="NA" & CDR_Cog5=="NA" & CDR_Cog6=="NA" & CDR_Cog7=="NA" & CDR_Cog8=="NA" & CDR_Cog9=="NA" & CDR_Cog10=="NA" & CDR_Cog11=="NA" & CDR_Cog12=="NA" & CDR_Cog13=="NA"~ "0",
        CDR_Cog1==0 & CDR_Cog2==0 & CDR_Cog3=="NA" & CDR_Cog4=="NA" & CDR_Cog5=="NA" & CDR_Cog6=="NA" & CDR_Cog7=="NA" & CDR_Cog8=="NA" & CDR_Cog9=="NA" & CDR_Cog10=="NA" & CDR_Cog11=="NA" & CDR_Cog12=="NA" & CDR_Cog13=="NA"~ "0",
        CDR_Cog1==0 & CDR_Cog2=="NA" & CDR_Cog3=="NA" & CDR_Cog4=="NA" & CDR_Cog5=="NA" & CDR_Cog6=="NA" & CDR_Cog7=="NA" & CDR_Cog8=="NA" & CDR_Cog9=="NA" & CDR_Cog10=="NA" & CDR_Cog11=="NA" & CDR_Cog12=="NA" & CDR_Cog13=="NA"~ "0",
        CDR_Cog1>=0.5 | CDR_Cog2>=0.5 | CDR_Cog3>=0.5 | CDR_Cog4>=0.5 | CDR_Cog5>=0.5 | CDR_Cog6>=0.5 | CDR_Cog7>=0.5 | CDR_Cog8>=0.5| CDR_Cog9>=0.5| CDR_Cog10>=0.5| CDR_Cog11>=0.5| CDR_Cog12>=0.5| CDR_Cog13>=0.5 ~ "1")) # if CDR >= 0.5 at ANY time-point, then converter

data1<-data1[!(data1$Converter==1),] # only keep people who are non converters
RevisionIDs_Keep<-data.frame(table(data1$ID))

CDRCog_MR_Revision<-CDRCog_MR[(CDRCog_MR$ID %in% RevisionIDs_Keep$Var1), ] # only keep long version data from people who are CDR0 at baseline and do NOT convert
CDRCog_MR_Revision_Wide<-data1 # rename this as wide data 



###### Revision Analyses ##### 
# For revision analyses with CDRCog_MR - only keep those IDs that are NOT in n=251 subset 
CDRCog_MR_Revision_Keep<-CDRCog_MR_Revision[!(CDRCog_MR_Revision$ID %in% df$ID), ]
CDRCog_MR_Revision_Keep_Wide<-CDRCog_MR_Revision_Wide[!(CDRCog_MR_Revision_Wide$ID %in% df$ID), ] #n=575

mean(CDRCog_MR_Revision_Keep_Wide$AGE)
sd(CDRCog_MR_Revision_Keep_Wide$AGE)

## Demographics - Age, Sex, Educ
CDRCog_MR_Revision_Keep_Wide <- merge(demo, CDRCog_MR_Revision_Keep_Wide, by="ID",all=F)%>%
  mutate(AGE = as.numeric(difftime(psy_date_Cog1, BIRTH, unit="weeks"))/52.25)%>%
  filter(AGE >= 60) # n=438 ** 

mean(CDRCog_MR_Revision_Keep_Wide$AGE)
sd(CDRCog_MR_Revision_Keep_Wide$AGE)

table(CDRCog_MR_Revision_Keep_Wide$SEX)

mean(CDRCog_MR_Revision_Keep_Wide$EDUC)
sd(CDRCog_MR_Revision_Keep_Wide$EDUC)

# Keep only 1 instance of CDR to get N+ values - 20% below 
CDRCog_MR_Revision_Keep<-CDRCog_MR_Revision_Keep[(CDRCog_MR_Revision_Keep$ID %in% CDRCog_MR_Revision_Keep_Wide$ID), ] #only keep those with IDs in keep_wide
CDRCog_MR_Neuropos<-CDRCog_MR_Revision_Keep%>% # Use these to create the neurodegenerative positivity scores!!!!!!!! 
  filter(Visit=="Cog1")

# Calculate neurodegration positivity - <20th percentile of that in group of individuals NOT in this paper that are cognitively normal
# Hippocampal volume 
HippoVol20percentile<-quantile(CDRCog_MR_Neuropos$PercentHIPP, probs=seq(0,1,0.2), na.rm=TRUE)
HippoVol20percentile # USE 0.459198 as the neurodegenerative cutoff for volume
# Cortical thickness 
HippoThick20percentile<-quantile(CDRCog_MR_Neuropos$Cort_Thick, probs=seq(0,1,0.2), na.rm=TRUE)
HippoThick20percentile # USE 2.741050  as the neurodegenerative cutoff for cortical thickness

# Create baseline PACC score to get m/sd to use in subsequent PACC analyses
msrtfree_Cog <- mean(CDRCog_MR_Revision_Keep_Wide$srtfree_Cog1, na.rm = T); sdsrtfree_Cog <- sd(CDRCog_MR_Revision_Keep_Wide$srtfree_Cog1, na.rm = T);CDRCog_MR_Revision_Keep_Wide$zsrtfree <- (CDRCog_MR_Revision_Keep_Wide$srtfree_Cog1 - msrtfree_Cog)/sdsrtfree_Cog
mlettnum_Cog <- mean(CDRCog_MR_Revision_Keep_Wide$lettnum_Cog1, na.rm = T); sdlettnum_Cog <- sd(CDRCog_MR_Revision_Keep_Wide$lettnum_Cog1, na.rm = T);CDRCog_MR_Revision_Keep_Wide$zlettnum <- (CDRCog_MR_Revision_Keep_Wide$lettnum_Cog1 - mlettnum_Cog)/sdlettnum_Cog
mANIMALS_Cog <- mean(CDRCog_MR_Revision_Keep_Wide$ANIMALS_Cog1, na.rm = T); sdANIMALS_Cog <- sd(CDRCog_MR_Revision_Keep_Wide$ANIMALS_Cog1, na.rm = T);CDRCog_MR_Revision_Keep_Wide$zANIMALS <- (CDRCog_MR_Revision_Keep_Wide$ANIMALS_Cog1 - mANIMALS_Cog)/sdANIMALS_Cog
mtmb_Cog <- mean(CDRCog_MR_Revision_Keep_Wide$tmb_Cog1, na.rm = T); sdtmb_Cog <- sd(CDRCog_MR_Revision_Keep_Wide$tmb_Cog1, na.rm = T);CDRCog_MR_Revision_Keep_Wide$ztmb <- ((CDRCog_MR_Revision_Keep_Wide$tmb_Cog1 - mtmb_Cog)/sdtmb_Cog)*(-1)
CDRCog_MR_Revision_Keep_Wide <- CDRCog_MR_Revision_Keep_Wide %>%
  mutate(PACC_Cog1 = (zsrtfree + zlettnum + zANIMALS + ztmb)/4)

# Create PACC score for subsequent time points using the m/sd from the first baseline visit 
# Time 2 
CDRCog_MR_Revision_Keep_Wide$zsrtfree <- (CDRCog_MR_Revision_Keep_Wide$srtfree_Cog2 - msrtfree_Cog)/sdsrtfree_Cog
CDRCog_MR_Revision_Keep_Wide$zlettnum <- (CDRCog_MR_Revision_Keep_Wide$lettnum_Cog2 - mlettnum_Cog)/sdlettnum_Cog
CDRCog_MR_Revision_Keep_Wide$zANIMALS <- (CDRCog_MR_Revision_Keep_Wide$ANIMALS_Cog2 - mANIMALS_Cog)/sdANIMALS_Cog
CDRCog_MR_Revision_Keep_Wide$ztmb <- ((CDRCog_MR_Revision_Keep_Wide$tmb_Cog2 - mtmb_Cog)/sdtmb_Cog)*(-1)
CDRCog_MR_Revision_Keep_Wide <- CDRCog_MR_Revision_Keep_Wide %>%
  mutate(PACC_Cog2 = (zsrtfree + zlettnum + zANIMALS + ztmb)/4)

# Time 3 
CDRCog_MR_Revision_Keep_Wide$zsrtfree <- (CDRCog_MR_Revision_Keep_Wide$srtfree_Cog3 - msrtfree_Cog)/sdsrtfree_Cog
CDRCog_MR_Revision_Keep_Wide$zlettnum <- (CDRCog_MR_Revision_Keep_Wide$lettnum_Cog3 - mlettnum_Cog)/sdlettnum_Cog
CDRCog_MR_Revision_Keep_Wide$zANIMALS <- (CDRCog_MR_Revision_Keep_Wide$ANIMALS_Cog3 - mANIMALS_Cog)/sdANIMALS_Cog
CDRCog_MR_Revision_Keep_Wide$ztmb <- ((CDRCog_MR_Revision_Keep_Wide$tmb_Cog3 - mtmb_Cog)/sdtmb_Cog)*(-1)
CDRCog_MR_Revision_Keep_Wide <- CDRCog_MR_Revision_Keep_Wide %>%
  mutate(PACC_Cog3 = (zsrtfree + zlettnum + zANIMALS + ztmb)/4)

# Time 4 
CDRCog_MR_Revision_Keep_Wide$zsrtfree <- (CDRCog_MR_Revision_Keep_Wide$srtfree_Cog4 - msrtfree_Cog)/sdsrtfree_Cog
CDRCog_MR_Revision_Keep_Wide$zlettnum <- (CDRCog_MR_Revision_Keep_Wide$lettnum_Cog4 - mlettnum_Cog)/sdlettnum_Cog
CDRCog_MR_Revision_Keep_Wide$zANIMALS <- (CDRCog_MR_Revision_Keep_Wide$ANIMALS_Cog4 - mANIMALS_Cog)/sdANIMALS_Cog
CDRCog_MR_Revision_Keep_Wide$ztmb <- ((CDRCog_MR_Revision_Keep_Wide$tmb_Cog4 - mtmb_Cog)/sdtmb_Cog)*(-1)
CDRCog_MR_Revision_Keep_Wide <- CDRCog_MR_Revision_Keep_Wide %>%
  mutate(PACC_Cog4 = (zsrtfree + zlettnum + zANIMALS + ztmb)/4)

# Time 5 
CDRCog_MR_Revision_Keep_Wide$zsrtfree <- (CDRCog_MR_Revision_Keep_Wide$srtfree_Cog5 - msrtfree_Cog)/sdsrtfree_Cog
CDRCog_MR_Revision_Keep_Wide$zlettnum <- (CDRCog_MR_Revision_Keep_Wide$lettnum_Cog5 - mlettnum_Cog)/sdlettnum_Cog
CDRCog_MR_Revision_Keep_Wide$zANIMALS <- (CDRCog_MR_Revision_Keep_Wide$ANIMALS_Cog5 - mANIMALS_Cog)/sdANIMALS_Cog
CDRCog_MR_Revision_Keep_Wide$ztmb <- ((CDRCog_MR_Revision_Keep_Wide$tmb_Cog5 - mtmb_Cog)/sdtmb_Cog)*(-1)
CDRCog_MR_Revision_Keep_Wide <- CDRCog_MR_Revision_Keep_Wide %>%
  mutate(PACC_Cog5 = (zsrtfree + zlettnum + zANIMALS + ztmb)/4)

# Time 6
CDRCog_MR_Revision_Keep_Wide$zsrtfree <- (CDRCog_MR_Revision_Keep_Wide$srtfree_Cog6 - msrtfree_Cog)/sdsrtfree_Cog
CDRCog_MR_Revision_Keep_Wide$zlettnum <- (CDRCog_MR_Revision_Keep_Wide$lettnum_Cog6 - mlettnum_Cog)/sdlettnum_Cog
CDRCog_MR_Revision_Keep_Wide$zANIMALS <- (CDRCog_MR_Revision_Keep_Wide$ANIMALS_Cog6 - mANIMALS_Cog)/sdANIMALS_Cog
CDRCog_MR_Revision_Keep_Wide$ztmb <- ((CDRCog_MR_Revision_Keep_Wide$tmb_Cog6 - mtmb_Cog)/sdtmb_Cog)*(-1)
CDRCog_MR_Revision_Keep_Wide <- CDRCog_MR_Revision_Keep_Wide %>%
  mutate(PACC_Cog6 = (zsrtfree + zlettnum + zANIMALS + ztmb)/4)

# Time 7
CDRCog_MR_Revision_Keep_Wide$zsrtfree <- (CDRCog_MR_Revision_Keep_Wide$srtfree_Cog7 - msrtfree_Cog)/sdsrtfree_Cog
CDRCog_MR_Revision_Keep_Wide$zlettnum <- (CDRCog_MR_Revision_Keep_Wide$lettnum_Cog7 - mlettnum_Cog)/sdlettnum_Cog
CDRCog_MR_Revision_Keep_Wide$zANIMALS <- (CDRCog_MR_Revision_Keep_Wide$ANIMALS_Cog7 - mANIMALS_Cog)/sdANIMALS_Cog
CDRCog_MR_Revision_Keep_Wide$ztmb <- ((CDRCog_MR_Revision_Keep_Wide$tmb_Cog7 - mtmb_Cog)/sdtmb_Cog)*(-1)
CDRCog_MR_Revision_Keep_Wide <- CDRCog_MR_Revision_Keep_Wide %>%
  mutate(PACC_Cog7 = (zsrtfree + zlettnum + zANIMALS + ztmb)/4)

# Time 8
CDRCog_MR_Revision_Keep_Wide$zsrtfree <- (CDRCog_MR_Revision_Keep_Wide$srtfree_Cog8 - msrtfree_Cog)/sdsrtfree_Cog
CDRCog_MR_Revision_Keep_Wide$zlettnum <- (CDRCog_MR_Revision_Keep_Wide$lettnum_Cog8 - mlettnum_Cog)/sdlettnum_Cog
CDRCog_MR_Revision_Keep_Wide$zANIMALS <- (CDRCog_MR_Revision_Keep_Wide$ANIMALS_Cog8 - mANIMALS_Cog)/sdANIMALS_Cog
CDRCog_MR_Revision_Keep_Wide$ztmb <- ((CDRCog_MR_Revision_Keep_Wide$tmb_Cog8 - mtmb_Cog)/sdtmb_Cog)*(-1)
CDRCog_MR_Revision_Keep_Wide <- CDRCog_MR_Revision_Keep_Wide %>%
  mutate(PACC_Cog8 = (zsrtfree + zlettnum + zANIMALS + ztmb)/4)

# Time 9
CDRCog_MR_Revision_Keep_Wide$zsrtfree <- (CDRCog_MR_Revision_Keep_Wide$srtfree_Cog9 - msrtfree_Cog)/sdsrtfree_Cog
CDRCog_MR_Revision_Keep_Wide$zlettnum <- (CDRCog_MR_Revision_Keep_Wide$lettnum_Cog9 - mlettnum_Cog)/sdlettnum_Cog
CDRCog_MR_Revision_Keep_Wide$zANIMALS <- (CDRCog_MR_Revision_Keep_Wide$ANIMALS_Cog9 - mANIMALS_Cog)/sdANIMALS_Cog
CDRCog_MR_Revision_Keep_Wide$ztmb <- ((CDRCog_MR_Revision_Keep_Wide$tmb_Cog9 - mtmb_Cog)/sdtmb_Cog)*(-1)
CDRCog_MR_Revision_Keep_Wide <- CDRCog_MR_Revision_Keep_Wide %>%
  mutate(PACC_Cog9 = (zsrtfree + zlettnum + zANIMALS + ztmb)/4)

# Time 10
CDRCog_MR_Revision_Keep_Wide$zsrtfree <- (CDRCog_MR_Revision_Keep_Wide$srtfree_Cog10 - msrtfree_Cog)/sdsrtfree_Cog
CDRCog_MR_Revision_Keep_Wide$zlettnum <- (CDRCog_MR_Revision_Keep_Wide$lettnum_Cog10 - mlettnum_Cog)/sdlettnum_Cog
CDRCog_MR_Revision_Keep_Wide$zANIMALS <- (CDRCog_MR_Revision_Keep_Wide$ANIMALS_Cog10 - mANIMALS_Cog)/sdANIMALS_Cog
CDRCog_MR_Revision_Keep_Wide$ztmb <- ((CDRCog_MR_Revision_Keep_Wide$tmb_Cog10 - mtmb_Cog)/sdtmb_Cog)*(-1)
CDRCog_MR_Revision_Keep_Wide <- CDRCog_MR_Revision_Keep_Wide %>%
  mutate(PACC_Cog10 = (zsrtfree + zlettnum + zANIMALS + ztmb)/4)

# Time 11
CDRCog_MR_Revision_Keep_Wide$zsrtfree <- (CDRCog_MR_Revision_Keep_Wide$srtfree_Cog11 - msrtfree_Cog)/sdsrtfree_Cog
CDRCog_MR_Revision_Keep_Wide$zlettnum <- (CDRCog_MR_Revision_Keep_Wide$lettnum_Cog11 - mlettnum_Cog)/sdlettnum_Cog
CDRCog_MR_Revision_Keep_Wide$zANIMALS <- (CDRCog_MR_Revision_Keep_Wide$ANIMALS_Cog11 - mANIMALS_Cog)/sdANIMALS_Cog
CDRCog_MR_Revision_Keep_Wide$ztmb <- ((CDRCog_MR_Revision_Keep_Wide$tmb_Cog11 - mtmb_Cog)/sdtmb_Cog)*(-1)
CDRCog_MR_Revision_Keep_Wide <- CDRCog_MR_Revision_Keep_Wide %>%
  mutate(PACC_Cog11 = (zsrtfree + zlettnum + zANIMALS + ztmb)/4)

# Time 12
CDRCog_MR_Revision_Keep_Wide$zsrtfree <- (CDRCog_MR_Revision_Keep_Wide$srtfree_Cog12 - msrtfree_Cog)/sdsrtfree_Cog
CDRCog_MR_Revision_Keep_Wide$zlettnum <- (CDRCog_MR_Revision_Keep_Wide$lettnum_Cog12 - mlettnum_Cog)/sdlettnum_Cog
CDRCog_MR_Revision_Keep_Wide$zANIMALS <- (CDRCog_MR_Revision_Keep_Wide$ANIMALS_Cog12 - mANIMALS_Cog)/sdANIMALS_Cog
CDRCog_MR_Revision_Keep_Wide$ztmb <- ((CDRCog_MR_Revision_Keep_Wide$tmb_Cog12 - mtmb_Cog)/sdtmb_Cog)*(-1)
CDRCog_MR_Revision_Keep_Wide <- CDRCog_MR_Revision_Keep_Wide %>%
  mutate(PACC_Cog12 = (zsrtfree + zlettnum + zANIMALS + ztmb)/4)

# Time 13
CDRCog_MR_Revision_Keep_Wide$zsrtfree <- (CDRCog_MR_Revision_Keep_Wide$srtfree_Cog13 - msrtfree_Cog)/sdsrtfree_Cog
CDRCog_MR_Revision_Keep_Wide$zlettnum <- (CDRCog_MR_Revision_Keep_Wide$lettnum_Cog13 - mlettnum_Cog)/sdlettnum_Cog
CDRCog_MR_Revision_Keep_Wide$zANIMALS <- (CDRCog_MR_Revision_Keep_Wide$ANIMALS_Cog13 - mANIMALS_Cog)/sdANIMALS_Cog
CDRCog_MR_Revision_Keep_Wide$ztmb <- ((CDRCog_MR_Revision_Keep_Wide$tmb_Cog13 - mtmb_Cog)/sdtmb_Cog)*(-1)
CDRCog_MR_Revision_Keep_Wide <- CDRCog_MR_Revision_Keep_Wide %>%
  mutate(PACC_Cog13 = (zsrtfree + zlettnum + zANIMALS + ztmb)/4)

# Time lag from baseline cognitive visit to each subsequent cognitive visit date 
CDRCog_MR_Revision_Keep_Wide<-CDRCog_MR_Revision_Keep_Wide%>%
  mutate(YrBl_Cog1 = as.numeric(difftime(psy_date_Cog1, psy_date_Cog1, unit="weeks"))/52.25,
         YrBl_Cog2 = as.numeric(difftime(psy_date_Cog2, psy_date_Cog1, unit="weeks"))/52.25, 
         YrBl_Cog3 = as.numeric(difftime(psy_date_Cog3, psy_date_Cog1, unit="weeks"))/52.25,
         YrBl_Cog4 = as.numeric(difftime(psy_date_Cog4, psy_date_Cog1, unit="weeks"))/52.25,
         YrBl_Cog5 = as.numeric(difftime(psy_date_Cog5, psy_date_Cog1, unit="weeks"))/52.25,
         YrBl_Cog6 = as.numeric(difftime(psy_date_Cog6, psy_date_Cog1, unit="weeks"))/52.25,
         YrBl_Cog7 = as.numeric(difftime(psy_date_Cog7, psy_date_Cog1, unit="weeks"))/52.25,
         YrBl_Cog8 = as.numeric(difftime(psy_date_Cog8, psy_date_Cog1, unit="weeks"))/52.25,
         YrBl_Cog9 = as.numeric(difftime(psy_date_Cog9, psy_date_Cog1, unit="weeks"))/52.25,
         YrBl_Cog10 = as.numeric(difftime(psy_date_Cog10, psy_date_Cog1, unit="weeks"))/52.25,
         YrBl_Cog11 = as.numeric(difftime(psy_date_Cog11, psy_date_Cog1, unit="weeks"))/52.25,
         YrBl_Cog12 = as.numeric(difftime(psy_date_Cog12, psy_date_Cog1, unit="weeks"))/52.25,
         YrBl_Cog13 = as.numeric(difftime(psy_date_Cog13, psy_date_Cog1, unit="weeks"))/52.25)

# PACC Revision
PACCRevision<-CDRCog_MR_Revision_Keep_Wide%>%
  dplyr::select(ID,118:130)%>%
  pivot_longer(cols=2:14,
               names_to="Visit",
               names_prefix = "Visit",
               values_to = "PACC")%>%
  na.omit()
PACCRevision$Visit<-gsub("[a-zA-Z ]", "", PACCRevision$Visit)
# YrBl 
YrBlRevision<-CDRCog_MR_Revision_Keep_Wide%>%
  dplyr::select(ID,131:143)%>%
  pivot_longer(cols=2:14,
               names_to="Visit",
               names_prefix = "Visit",
               values_to = "YrBl")%>%
  na.omit()
YrBlRevision$Visit<-gsub("[a-zA-Z ]", "", YrBlRevision$Visit)
CDRCog_MR_Revision_Keep_Long<-merge(PACCRevision,YrBlRevision,by=c("ID","Visit"), all=T)
write.csv(CDRCog_MR_Revision_Keep_Long, "./Data/N151_Revision_Long_5-13-22.csv", row.names=FALSE)

CDRCog_MR_Revision_Keep_Long <- lmer(PACC ~  (1 + YrBl|ID), data=CDRCog_MR_Revision_Keep_Long, contrasts=c("contr.sum", "contr.poly"), REML=FALSE)
summary(CDRCog_MR_Revision_Keep_Long)
write.csv(coef(CDRCog_MR_Revision_Keep_Long)$ID, './Tables/long2_coefficients_PACC_YrBl_5-13-22.csv', row.names = TRUE)

CogNormalCoef<-read.csv("./Tables/long2_coefficients_PACC_YrBl_5-13-22.csv")%>%
  rename(ID=X)
mCogNormalCoef<-mean(CogNormalCoef$YrBl)
mCogNormalCoef
sdCogNormalCoef<-sd(CogNormalCoef$YrBl)
sdCogNormalCoef

####### Cognitive decliners with new dataset that is not included in this group 
mCogNormalCoef
sdCogNormalCoef
msd<-mCogNormalCoef + sdCogNormalCoef
msd

CogNormalCoef_Below<-mCogNormalCoef - sdCogNormalCoef
CogNormalCoef_Below # -0.0406735
CogNormalCoef_Above<-mCogNormalCoef + sdCogNormalCoef
CogNormalCoef_Above

