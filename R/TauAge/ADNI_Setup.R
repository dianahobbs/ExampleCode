# Title: ADNI Cohort Organization for Brian's Tau Grant 
# Author: Diana Hobbs
# Date: November 2022

#### Set Up ####
# Load packages  
pacman::p_load(tidyverse,jtools,naniar,lubridate,zoo)

###### Load data and select columns of interest
# ***first open in excel and format dates as 3/14/2012 and THEN save to csv

# Demographics 
demo <- read.csv("./data_files/PTDEMOG.csv")%>% # n = 4636
  replace_with_na(replace=list(PTGENDER=c(-4), PTEDUCAT=c(-1,-4), PTRACCAT=c(-4)))%>%
  mutate(PTDOBMM = sprintf("%02.f", PTDOBMM))%>%
  mutate(BIRTH = as.yearmon(paste(PTDOBYY, PTDOBMM), "%Y %m"))%>%
  #PTGENDER: 1=Male; 2=Female
  mutate(PTGENDER = str_replace_all(PTGENDER, c("1"="Male", "2"="Female")))%>%
  #PTRACCAT: 1=American Indian or Alaskan Native; 2=Asian; 3=Native Hawaiian or Other Pacific Islander;
  #          4=Black or African American; 5=White; 6=More than one race; 7=Unknown
  mutate(PTRACCAT = str_replace_all(PTRACCAT,c("1"="American Indian or Alaskan Native","2"="Asian","3"="Native Hawaiian or Other Pacific Islander",
                                               "4"="Black or African American","5"="White","6"="More than one race","7"="Unknown")))%>%
  rename(sex=PTGENDER, educ=PTEDUCAT, race=PTRACCAT)%>%
  dplyr::select(RID,BIRTH,sex,educ,race)%>%
  na.omit()
demo$race[is.na(demo$race)] <- "Unknown"
demo <- demo%>%mutate(race.binary = case_when(race == "White" ~ 0,
                                              race != "White" ~ 1),
                      Female = case_when(sex == "Male" ~ 0,
                                         sex == "Female" ~ 1))
demo <- unique(demo) # n = 4046

### Imaging data 
# Tau
tau <- read.csv("./data_files/UCBERKELEYAV1451_04_26_22.csv")%>% # n = 1450
  mutate(tau.date = as.Date(EXAMDATE, "%m/%d/%Y"))%>%
  rowwise()%>%
  mutate("bankssts"=mean(c(CTX_LH_BANKSSTS_SUVR,CTX_RH_BANKSSTS_SUVR)),
         "caudal anterior cingulate"=mean(c(CTX_LH_CAUDALANTERIORCINGULATE_SUVR,CTX_RH_CAUDALANTERIORCINGULATE_SUVR)),
         "caudal middle frontal"=mean(c(CTX_LH_CAUDALMIDDLEFRONTAL_SUVR,CTX_RH_CAUDALMIDDLEFRONTAL_SUVR)),
         "cuneus"=mean(c(CTX_LH_CUNEUS_SUVR,CTX_RH_CUNEUS_SUVR)),
         "entorhinal"=mean(c(CTX_LH_ENTORHINAL_SUVR,CTX_RH_ENTORHINAL_SUVR)),
         "frontal pole"=mean(c(CTX_LH_FRONTALPOLE_SUVR,CTX_RH_FRONTALPOLE_SUVR)),
         "fusiform"=mean(c(CTX_LH_FUSIFORM_SUVR,CTX_RH_FUSIFORM_SUVR)),
         "inferior parietal"=mean(c(CTX_LH_INFERIORPARIETAL_SUVR,CTX_RH_INFERIORPARIETAL_SUVR)),
         "inferior temporal"=mean(c(CTX_LH_INFERIORTEMPORAL_SUVR,CTX_RH_INFERIORTEMPORAL_SUVR)),
         "insula"=mean(c(CTX_LH_INSULA_SUVR,CTX_RH_INSULA_SUVR)),
         "isthmus cingulate"=mean(c(CTX_LH_ISTHMUSCINGULATE_SUVR,CTX_RH_ISTHMUSCINGULATE_SUVR)),
         "lateral occipital"=mean(c(CTX_LH_LATERALOCCIPITAL_SUVR,CTX_RH_LATERALOCCIPITAL_SUVR)),
         "lateral orbitofrontal"=mean(c(CTX_LH_LATERALORBITOFRONTAL_SUVR,CTX_RH_LATERALORBITOFRONTAL_SUVR)),
         "lingual"=mean(c(CTX_LH_LINGUAL_SUVR,CTX_RH_LINGUAL_SUVR)),
         "medial orbitofrontal"=mean(c(CTX_LH_MEDIALORBITOFRONTAL_SUVR,CTX_RH_MEDIALORBITOFRONTAL_SUVR)),
         "middle temporal"=mean(c(CTX_LH_MIDDLETEMPORAL_SUVR,CTX_RH_MIDDLETEMPORAL_SUVR)),
         "paracentral"=mean(c(CTX_LH_PARACENTRAL_SUVR,CTX_RH_PARACENTRAL_SUVR)),
         "parahippocampal"=mean(c(CTX_LH_PARAHIPPOCAMPAL_SUVR,CTX_RH_PARAHIPPOCAMPAL_SUVR)),
         "pars opercularis"=mean(c(CTX_LH_PARSOPERCULARIS_SUVR,CTX_RH_PARSOPERCULARIS_SUVR)),
         "pars orbitalis"=mean(c(CTX_LH_PARSORBITALIS_SUVR,CTX_RH_PARSORBITALIS_SUVR)),
         "pars triangularis"=mean(c(CTX_LH_PARSTRIANGULARIS_SUVR,CTX_RH_PARSTRIANGULARIS_SUVR)),
         "pericalcarine"=mean(c(CTX_LH_PERICALCARINE_SUVR,CTX_RH_PERICALCARINE_SUVR)),
         "postcentral"=mean(c(CTX_LH_POSTCENTRAL_SUVR,CTX_RH_POSTCENTRAL_SUVR)),
         "posterior cingulate"=mean(c(CTX_LH_POSTERIORCINGULATE_SUVR,CTX_RH_POSTERIORCINGULATE_SUVR)),
         "precentral"=mean(c(CTX_LH_PRECENTRAL_SUVR,CTX_RH_PRECENTRAL_SUVR)),
         "precuneus"=mean(c(CTX_LH_PRECUNEUS_SUVR,CTX_RH_PRECUNEUS_SUVR)),
         "rostral anterior cingulate"=mean(c(CTX_LH_ROSTRALANTERIORCINGULATE_SUVR,CTX_RH_ROSTRALANTERIORCINGULATE_SUVR)),
         "rostral middle frontal"=mean(c(CTX_LH_ROSTRALMIDDLEFRONTAL_SUVR,CTX_RH_ROSTRALMIDDLEFRONTAL_SUVR)),
         "superior frontal"=mean(c(CTX_LH_SUPERIORFRONTAL_SUVR,CTX_RH_SUPERIORFRONTAL_SUVR)),
         "superior parietal"=mean(c(CTX_LH_SUPERIORPARIETAL_SUVR,CTX_RH_SUPERIORPARIETAL_SUVR)),
         "superior temporal"=mean(c(CTX_LH_SUPERIORTEMPORAL_SUVR,CTX_RH_SUPERIORTEMPORAL_SUVR)),
         "supramarginal"=mean(c(CTX_LH_SUPRAMARGINAL_SUVR,CTX_RH_SUPRAMARGINAL_SUVR)),
         "temporal pole"=mean(c(CTX_LH_TEMPORALPOLE_SUVR,CTX_RH_TEMPORALPOLE_SUVR)),
         "transverse temporal"=mean(c(CTX_LH_TRANSVERSETEMPORAL_SUVR,CTX_RH_TRANSVERSETEMPORAL_SUVR)),
         "amygdala"=mean(c(LEFT_AMYGDALA_SUVR,RIGHT_AMYGDALA_SUVR)),
         "caudate"=mean(c(LEFT_CAUDATE_SUVR,RIGHT_CAUDATE_SUVR)),
         "hippocampus"=mean(c(LEFT_HIPPOCAMPUS_SUVR,RIGHT_HIPPOCAMPUS_SUVR)),
         "pallidum"=mean(c(LEFT_PALLIDUM_SUVR,RIGHT_PALLIDUM_SUVR)),
         "putamen"=mean(c(LEFT_PUTAMEN_SUVR,RIGHT_PUTAMEN_SUVR)),
         "thalamus proper"=mean(c(LEFT_THALAMUS_PROPER_SUVR,RIGHT_THALAMUS_PROPER_SUVR)),
         "ventral DC"=mean(c(LEFT_VENTRALDC_SUVR,RIGHT_VENTRALDC_SUVR)))%>%
  # Meta-temporal ROI: L/R Amygdala, L/R Entorhinal, L/R Fusiform, L/R Inferior Temporal, L/R Middle Temporal
  dplyr::select(RID,VISCODE2,tau.date,"bankssts","caudal anterior cingulate","caudal middle frontal","cuneus","entorhinal",
                "frontal pole","fusiform","inferior parietal","inferior temporal","insula","isthmus cingulate","lateral occipital",
                "lateral orbitofrontal","lingual","medial orbitofrontal","middle temporal","paracentral","parahippocampal",
                "pars opercularis","pars orbitalis","pars triangularis","pericalcarine","postcentral","posterior cingulate",
                "precentral","precuneus","rostral anterior cingulate","rostral middle frontal","superior frontal",
                "superior parietal","superior temporal","supramarginal","temporal pole","transverse temporal","amygdala",
                "caudate","hippocampus","pallidum","putamen","thalamus proper","ventral DC")

           
# Amyloid 
av45 <- read.csv("./data_files/UCBERKELEYAV45_04_26_22.csv")%>% # n = 3086
  mutate(av45.date = as.Date(EXAMDATE, "%m/%d/%Y"),
         SUMMARYSUVR_WHOLECEREBNORM_1.11CUTOFF = case_when(
           SUMMARYSUVR_WHOLECEREBNORM_1.11CUTOFF == 0 ~ "A-",
           SUMMARYSUVR_WHOLECEREBNORM_1.11CUTOFF == 1 ~ "A+"),
         SUMMARYSUVR_COMPOSITE_REFNORM_0.78CUTOFF = case_when(
           SUMMARYSUVR_COMPOSITE_REFNORM_0.78CUTOFF == 0 ~ "A-",
           SUMMARYSUVR_COMPOSITE_REFNORM_0.78CUTOFF == 1 ~ "A+"))%>%
  # Cross-sectional: SUMMARYSUVR_WHOLECEREBNORM: normalized whole cerebellum with 1.11 threshold
  # Longitudinal: SUMMARYSUVR_COMPOSITE_REFNORM: already FS normalized with 0.78 threshold
  rename(av45_cross.sectional = SUMMARYSUVR_WHOLECEREBNORM,
         av45_APos_cross = SUMMARYSUVR_WHOLECEREBNORM_1.11CUTOFF,
         av45_longitudinal = SUMMARYSUVR_COMPOSITE_REFNORM,
         av45_APos_long = SUMMARYSUVR_COMPOSITE_REFNORM_0.78CUTOFF)%>%
  dplyr::select(RID,av45.date,av45_cross.sectional,av45_APos_cross,av45_longitudinal,av45_APos_long) # TP, VISCODE2

# CDR - USERDATE has more data AND the dates are within a year of EXAMDATE !!*** USE USERDATE !!!!*****
cdr <- read.csv("./data_files/CDR.csv")%>% # n = 12731
  mutate(cdr.date = as.Date(USERDATE, "%m/%d/%Y"),
         CDGLOBAL = case_when(CDGLOBAL == 0 ~ 0,
                              CDGLOBAL > 0 ~ 1))%>%
  replace_with_na(replace=list(CDGLOBAL=c(-1)))%>%
  rename(cdr=CDGLOBAL)%>%
  dplyr::select(RID,cdr.date,cdr)%>%
  na.omit()

### Cognitive tests
# Animal: fluency animals (total correct) "CATANIMSC" from Neurobat.csv (REMOVE -1 FOR NO DATA)
# TrailsB: time to complete "TRABSCOR" from Neurobat.csv
# WAIS: digit symbol substitute "WAISR_Score from Item.csv
# Memu (????) 
cog <- read.csv("./NEUROBAT.csv")%>% # n = 10949
  replace_with_na(replace=list(CATANIMSC=-1, TRABSCOR=-1))%>%
  mutate(psy.date = as.Date(USERDATE, "%m/%d/%Y"),
         tmb = case_when(TRABSCOR > 180 ~ 180,
                         TRABSCOR <= 180 ~ TRABSCOR))%>%
  rename(ANIMALS=CATANIMSC)%>%
  dplyr::select(RID, psy.date,ANIMALS,tmb)%>%
  na.omit()

# mem <- read.csv("./UWNPSYCHSUM_12_13_21.csv")%>% # n = 11457
#   mutate(USERDATE = as.Date(USERDATE, "%m/%d/%Y"))%>%
#   select(RID, USERDATE, ADNI_MEM)
# 
# wais <- read.csv("./ITEM.csv")%>% # n = 1395
#   select(RID, WAISR_ExamDate, WAISR_Score)%>%
#   mutate(WAISR_Score = str_replace_all(WAISR_Score, c("NULL" = "")),
#          WAISR_Score = as.numeric(WAISR_Score),
#          WAISR_ExamDate = as.Date(WAISR_ExamDate, "%m/%d/%Y"))%>%
#   na.omit()
# 
# merger <- merge(animtrail, mem, by = c("RID","USERDATE"), all=TRUE) # n = 12494
# psych <- merge(merger, wais, by = "RID", all=TRUE) # n = 18271
# psych <- psych%>% # n = 7152
#   mutate(datediff = as.numeric(difftime(USERDATE, WAISR_ExamDate, unit="weeks"))/52.25)%>%
#   filter(between(datediff,-3,3))%>%
#   select(-datediff, -WAISR_ExamDate)

CDR_cog <- merge(cdr,cog,by="RID")%>% # n = 2396
  mutate(cdr.cog.date.diff = abs(as.numeric(difftime(cdr.date, psy.date, unit="weeks"))/52.25))%>%
  group_by(RID)%>%
  arrange(cdr.cog.date.diff)%>%
  filter((row_number()==1))%>%
  dplyr::select(-cdr.cog.date.diff, -psy.date)%>%
  na.omit()

##### Merge files together 
demo_tau <- merge(demo,tau,by.y="RID")%>% # n = 1450
  mutate(age = as.numeric(difftime(tau.date, BIRTH, unit="weeks"))/52.25)%>%
  dplyr::select(RID,age,sex,Female,educ,race,race.binary,VISCODE2,tau.date,
                "bankssts","caudal anterior cingulate","caudal middle frontal","cuneus","entorhinal",
                "frontal pole","fusiform","inferior parietal","inferior temporal","insula","isthmus cingulate","lateral occipital",
                "lateral orbitofrontal","lingual","medial orbitofrontal","middle temporal","paracentral","parahippocampal",
                "pars opercularis","pars orbitalis","pars triangularis","pericalcarine","postcentral","posterior cingulate",
                "precentral","precuneus","rostral anterior cingulate","rostral middle frontal","superior frontal",
                "superior parietal","superior temporal","supramarginal","temporal pole","transverse temporal","amygdala",
                "caudate","hippocampus","pallidum","putamen","thalamus proper","ventral DC")

av45_tau <- merge(demo_tau, av45, by="RID")%>% # n = 946
  group_by(RID)%>%
  distinct(tau.date, .keep_all = T)%>%
  ungroup()%>%
  mutate(av45.tau.date.diff = abs(as.numeric(difftime(av45.date, tau.date, unit="weeks"))/52.25)) 

cdr_tau <- merge(demo_tau, CDR_cog, by="RID")%>% # n = 1432
  group_by(RID)%>%
  distinct(tau.date, .keep_all = T)%>%
  ungroup()%>%
  mutate(cdr.tau.date.diff = abs(as.numeric(difftime(cdr.date, tau.date, unit="weeks"))/52.25))
cdr.sub <- cdr_tau%>%dplyr::select(RID,VISCODE2,cdr.date,cdr,cdr.tau.date.diff,ANIMALS,tmb)

data <- merge(av45_tau,cdr.sub,by=c("RID","VISCODE2"))%>% # n = 933
  mutate(av45.cdr.date.diff = abs(as.numeric(difftime(av45.date, cdr.date, unit="weeks"))/52.25))%>%
  mutate(group.cross.sectional = case_when(cdr == 0 & av45_APos_cross == "A-" ~ "Control",
                                           cdr == 0 & av45_APos_cross == "A+" ~ "Preclinical",
                                           cdr == 1 & av45_APos_cross == "A-" ~ "Cognitively Impaired, A-",
                                           cdr == 1 & av45_APos_cross == "A+" ~ "Symptomatic"),
         group.longitudinal = case_when(cdr == 0 & av45_APos_long == "A-" ~ "Control",
                                        cdr == 0 & av45_APos_long == "A+" ~ "Preclinical",
                                        cdr == 1 & av45_APos_long == "A-" ~ "Cognitively Impaired, A-",
                                        cdr == 1 & av45_APos_long == "A+" ~ "Symptomatic"))

filter_av45.tau <- data%>%filter(between(av45.tau.date.diff,-3,3)) # n = 610
#filter_cdr.tau <- data%>%filter(between(cdr.tau.date.diff,-3,3)) # n = 536
#filter_av45.cdr <- data%>%filter(between(av45.cdr.date.diff,-3,3)) # n = 613

# keep only baseline 
bl_data <- data%>% # n = 539
  group_by(RID)%>%
  arrange(tau.date)%>%
  filter(row_number()==1)%>%
  ungroup()

bl_filter_av45.tau <- filter_av45.tau%>% # n = 407
  group_by(RID)%>%
  arrange(tau.date)%>%
  filter(row_number()==1)%>%
  ungroup()
# bl_filter_cdr.tau <- filter_cdr.tau%>% # n = 358
#   group_by(RID)%>%
#   arrange(tau.date)%>%
#   filter(row_number()==1)%>%
#   ungroup()
# bl_filter_av45.cdr <- filter_av45.cdr%>% # n = 365
#   group_by(RID)%>%
#   arrange(tau.date)%>%
#   filter(row_number()==1)%>%
#   ungroup()

filter_all <-filter_av45.tau%>%filter(between(av45.cdr.date.diff,-3,3))
filter_all <-filter_all%>%filter(between(cdr.tau.date.diff,-3,3))
bl_filter_all <- filter_all%>% # n = 264
  group_by(RID)%>%
  arrange(tau.date)%>%
  filter(row_number()==1)%>%
  ungroup()

########### NEXT STEPS 
# write.csv(data,"./all_data_n946.csv",row.names = FALSE)
# write.csv(bl_data,"./baseline_data_n546.csv",row.names = FALSE)
# write.csv(filter_av45.tau,"./filter.av45.tau_n623.csv",row.names = FALSE)
# write.csv(bl_filter_av45.tau,"./baseline.filter.av45.tau_n414.csv",row.names = FALSE)
write.csv(bl_filter_av45.tau,"./neuropsych_baseline.filter.av45.tau_n407.csv",row.names = FALSE)
# write.csv(filter_cdr.tau,"./filter.cdr.tau_n536.csv",row.names = FALSE)
# write.csv(bl_filter_cdr.tau,"./baseline.filter.cdr.tau_n358.csv",row.names = FALSE)
# write.csv(filter_av45.cdr,"./filter.av45.cdr_n613.csv",row.names = FALSE)
# write.csv(bl_filter_av45.cdr,"./baseline.filter.av45.cdr_n365.csv",row.names = FALSE)
# write.csv(filter_all,"./filter_all_n399.csv",row.names = FALSE)
# write.csv(bl_filter_all,"./baseline.filter_all_n271.csv",row.names = FALSE)
