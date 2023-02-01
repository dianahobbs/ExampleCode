# Title: ADRC Cohort Organization for Brian's Tau Grant 
# Author: Diana Hobbs
# Date: November 2022

#### Set Up ####
# Load packages  
pacman::p_load(tidyverse,jtools,naniar,lubridate,zoo)

###### Load data and select columns of interest
# ***first open in excel and format dates as 3/14/2012 and THEN save to csv

# Demographics 
demo <- read.csv("./datafiles/mod_demographics.csv")%>% # n = 3632
  replace_with_na(replace=list(BIRTH = "###############################################################################################################################################################################################################################################################"))%>%
  mutate(birth = as.Date(BIRTH, "%m/%d/%Y"),
         race = str_replace_all(race,c("AIAN"="Asian","ASIAN"="Asian")),
         sex = str_replace_all(sex, c("M"="Male", "F"="Female")))%>%
  rename(educ=EDUC)%>%
  dplyr::select(ID,birth, sex, race)
demo <- demo%>%mutate(race.binary = case_when(race == "White" ~ 0,
                                              race != "White" ~ 1),
                      Female = case_when(sex == "Male" ~ 0,
                                         sex == "Female" ~ 1))
demo <- unique(demo) # n = 3632

### Imaging data 
# Tau
tau<-read.csv("./datafiles/tau.csv")%>% # n = 611
  mutate(tau.date = as.Date(PET_Date, "%m/%d/%Y"))%>%
  #mutate(Tau_metaROI=(av1451_fsuvr_tot_ctx_entorhinal+av1451_fsuvr_tot_amygdala+av1451_fsuvr_tot_ctx_fusiform+av1451_fsuvr_tot_ctx_infertmp+av1451_fsuvr_tot_ctx_midtmp+av1451_fsuvr_tot_ctx_parahpcmpl)/6)%>%
  rename("bankssts"=av1451_fsuvr_rsf_tot_ctx_sstsban,
         "caudal anterior cingulate"=av1451_fsuvr_rsf_tot_ctx_caudant,
         "caudal middle frontal"=av1451_fsuvr_rsf_tot_ctx_caudmid,
         "corpus callosum"=av1451_fsuvr_rsf_tot_ctx_crpclm,
         "cuneus"=av1451_fsuvr_rsf_tot_ctx_cuneus,
         "entorhinal"=av1451_fsuvr_rsf_tot_ctx_entorhi,
         "frontal pole"=av1451_fsuvr_rsf_tot_ctx_frnpole,
         "fusiform"=av1451_fsuvr_rsf_tot_ctx_fusifor,
         "inferior parietal"=av1451_fsuvr_rsf_tot_ctx_inferpr,
         "inferior temporal"=av1451_fsuvr_rsf_tot_ctx_infertm,
         "insula"=av1451_fsuvr_rsf_tot_ctx_insula,
         "isthmus cingulate"=av1451_fsuvr_rsf_tot_ctx_isthmus,
         "lateral occipital"=av1451_fsuvr_rsf_tot_ctx_latocc,
         "lateral orbitofrontal"=av1451_fsuvr_rsf_tot_ctx_latorbf,
         "lingual"=av1451_fsuvr_rsf_tot_ctx_lingual,
         "medial orbitofrontal"=av1451_fsuvr_rsf_tot_ctx_medorbf,
         "middle temporal"=av1451_fsuvr_rsf_tot_ctx_midtmp,
         "paracentral"=av1451_fsuvr_rsf_tot_ctx_paracnt,
         "parahippocampal"=av1451_fsuvr_rsf_tot_ctx_parahpc,
         "pars opercularis"=av1451_fsuvr_rsf_tot_ctx_parsopc,
         "pars orbitalis"=av1451_fsuvr_rsf_tot_ctx_parsorb,
         "pars triangularis"=av1451_fsuvr_rsf_tot_ctx_parstrn,
         "pericalcarine"=av1451_fsuvr_rsf_tot_ctx_periclc,
         "postcentral"=av1451_fsuvr_rsf_tot_ctx_postcnt,
         "posterior cingulate"=av1451_fsuvr_rsf_tot_ctx_postcng,
         "precentral"=av1451_fsuvr_rsf_tot_ctx_precntr,
         "precuneus"=av1451_fsuvr_rsf_tot_ctx_precune,
         "rostral anterior cingulate"=av1451_fsuvr_rsf_tot_ctx_rosantc,
         "rostral middle frontal"=av1451_fsuvr_rsf_tot_ctx_rosmidf,
         "superior frontal"=av1451_fsuvr_rsf_tot_ctx_superfr,
         "superior parietal"=av1451_fsuvr_rsf_tot_ctx_superpr,
         "superior temporal"=av1451_fsuvr_rsf_tot_ctx_supertm,
         "supramarginal"=av1451_fsuvr_rsf_tot_ctx_supramr,
         "temporal pole"=av1451_fsuvr_rsf_tot_ctx_tmppole,
         "transverse temporal"=av1451_fsuvr_rsf_tot_ctx_transtm,
         "amygdala"=av1451_fsuvr_rsf_tot_amygdala,
         "caudate"=av1451_fsuvr_rsf_tot_caud,
         "hippocampus" =av1451_fsuvr_rsf_tot_hippocampus,
         "pallidum"=av1451_fsuvr_rsf_tot_pallidum,
         "putamen"=av1451_fsuvr_rsf_tot_putamen,
         "thalamus proper"=av1451_fsuvr_rsf_tot_thalamus_pr,
         "ventral DC"=av1451_fsuvr_rsf_tot_ventraldc)%>%
  dplyr::select(ID, tau.date, "bankssts","caudal anterior cingulate","caudal middle frontal","cuneus","entorhinal",
         "frontal pole","fusiform","inferior parietal","inferior temporal","insula","isthmus cingulate","lateral occipital",
         "lateral orbitofrontal","lingual","medial orbitofrontal","middle temporal","paracentral","parahippocampal",
         "pars opercularis","pars orbitalis","pars triangularis","pericalcarine","postcentral","posterior cingulate",
         "precentral","precuneus","rostral anterior cingulate","rostral middle frontal","superior frontal",
         "superior parietal","superior temporal","supramarginal","temporal pole","transverse temporal","amygdala",
         "caudate","hippocampus","pallidum","putamen","thalamus proper","ventral DC")%>%
    na.omit()

# Amyloid 
# PiB & Centiloid Conversion, eq: 111.8*PIB_3060_SUVR_RSF - 119.3
pib <- read.csv("./datafiles/pib.csv")%>% # n = 1235
  mutate(pib.date = as.Date(PET_Date, "%m/%d/%Y"))%>%
  mutate(pib_cent=111.8*pib_fsuvr_tot_cortmean-119.3)%>%
  dplyr::select(ID,Tracer,pib.date,pib_cent)%>%
  na.omit()

# av45 & Centiloid Conversion, eq: 163.6*PIB_3060_SUVR_RSF - 181.0
av45 <- read.csv("./datafiles/av45.csv")%>% # n = 739
  mutate(av45.date = as.Date(PET_Date, "%m/%d/%Y"))%>%
  mutate(av45_cent = 163.6*av45_fsuvr_tot_cortmean - 181.0)%>%
  dplyr::select(ID,Tracer,av45.date,av45_cent)%>%
  na.omit()

# Merge Pib and av45 
amyloid<-merge(pib,av45,by=c("ID", "Tracer"), all=T)%>% # n = 1974
  mutate(pib_cent = replace_na(pib_cent,0), # code NAs as 0
         av45_cent = replace_na(av45_cent,0))%>% 
  mutate(Tracer = case_when(
    pib_cent != "0" ~ "pib", 
    pib_cent == "0" ~ "av45",
    av45_cent == "0" ~ "pib"),
    amy.date = case_when(
      Tracer == "pib" ~ pib.date,
      Tracer == "av45" ~ av45.date))%>%
  mutate(centiloid = case_when(
    pib_cent != "0" ~ pib_cent,
    pib_cent == "0" ~ av45_cent,
    av45_cent == "0" ~ pib_cent))%>%
  mutate(amy_pos = case_when(
    Tracer=="pib" & centiloid >= 27.1 ~ "A+",
    Tracer=="pib" & centiloid < 27.1 ~ "A-",
    Tracer=="av45" & centiloid >= 21.9 ~ "A+",
    Tracer=="av45" & centiloid < 21.9 ~ "A-"))%>%
  dplyr::select(ID,amy.date,Tracer,centiloid,amy_pos)

### CDR 
cdr <- read.csv("./datafiles/mod_b4_cdr.csv")%>% # n = 18455
  mutate(cdr.date = as.Date(TESTDATE, "%m/%d/%Y"),
         cdr = case_when(cdr == 0 ~ 0,
                         cdr > 0 ~ 1))%>%
  dplyr::select(ID,cdr.date,cdr)%>%
  na.omit()

### Cognitive tests
# Knight ADRC: ANIMALS, TrailsB (tmb), Digit Symbol (digsym), srtfree 
# Removing digsym and srtfree
cog <- read.csv("./datafiles/mod_psychometrics.csv")%>% # n = 16709
  mutate(psy.date = as.Date(psy_date, "%m/%d/%Y"))%>%
  dplyr::select(ID,psy.date,ANIMALS,tmb)

CDR_cog <- merge(cdr,cog,by="ID")%>% # n = 3625
  mutate(cdr.cog.date.diff = abs(as.numeric(difftime(cdr.date, psy.date, unit="weeks"))/52.25))%>%
  group_by(ID)%>%
  arrange(cdr.cog.date.diff)%>%
  filter((row_number()==1))%>%
  dplyr::select(-cdr.cog.date.diff, -psy.date)
CDR_cog <- CDR_cog[!rowSums(is.na(CDR_cog)) >= 2, ] # n = 3011
# CDR_cog <- CDR_cog%>%mutate(wais = NA, memu = NA)%>% # n = 2709
#                      group_by(ID)%>%
#                      arrange(cdr.cog.date.diff)%>%
#                      filter((row_number()==1))%>%
#                      dplyr::select(-cdr.cog.date.diff, -psy.date)
         
##### Merge files 
# demo & tau
demo_tau <- merge(demo,tau,by.y="ID")%>% # n = 611
  mutate(age = as.numeric(difftime(tau.date, birth, unit="weeks"))/52.25)%>%
  dplyr::select(ID, age, sex, Female, race, race.binary, tau.date, 
                "bankssts","caudal anterior cingulate","caudal middle frontal",
                "cuneus","entorhinal","frontal pole","fusiform","inferior parietal","inferior temporal",
                "insula","isthmus cingulate","lateral occipital","lateral orbitofrontal","lingual",
                "medial orbitofrontal","middle temporal","paracentral","parahippocampal","pars opercularis",
                "pars orbitalis","pars triangularis","pericalcarine","postcentral","posterior cingulate",
                "precentral","precuneus","rostral anterior cingulate","rostral middle frontal","superior frontal",
                "superior parietal","superior temporal","supramarginal","temporal pole","transverse temporal",
                "amygdala","caudate","hippocampus","pallidum","putamen","thalamus proper","ventral DC")

# demo_tau and amy
tau_amy <- merge(demo_tau, amyloid, by="ID")%>% # n = 580
  group_by(ID)%>%
  distinct(tau.date, .keep_all = T)%>%
  ungroup()%>%
  mutate(amy.tau.date.diff = abs(as.numeric(difftime(amy.date, tau.date, unit="weeks"))/52.25))


# demo_tau_amy and cdr
cdr_tau <- merge(tau_amy, CDR_cog, by="ID")%>% # n = 578 
  group_by(ID)%>%
  distinct(tau.date, .keep_all = T)%>%
  ungroup()%>%
  mutate(cdr.tau.date.diff = abs(as.numeric(difftime(cdr.date, tau.date, unit="weeks"))/52.25))

##### %%%%%%%%%%%%%%%%%%%%%% Grabbing Brian's Data %%%%%%%%%%%%%%%%%%%%%% #####
brian <- read.csv("./datafiles/Brian_Orig_Data.csv")%>% # n = 144
  mutate(sex = str_replace_all(sex, c("female"="Female", "male"="Male", "FeMale"="Female")),
         race.binary = case_when(Race == "Caucasian" ~ 0,
                                 Race != "Caucasian" ~ 1),
         amy_pos = case_when(Tracer=="PIB" & centiloids >= 27.1 ~ "A+",
                             Tracer=="PIB" & centiloids < 27.1 ~ "A-",
                             Tracer=="AV45" & centiloids >= 21.9 ~ "A+",
                             Tracer=="AV45" & centiloids < 21.9 ~ "A-"),
         Female = case_when(sex == "Male" ~ 0, sex == "Female" ~ 1),
         cdr.date = NA, 
         cdr = 0,
         ANIMALS = NA, 
         tmb = NA, 
         digsym = NA, 
         srtfree = NA, 
         wais = NA, 
         memu = NA,
         cdr.tau.date.diff = NA,
         amy.date = as.Date(AB.PET.Date, "%m/%d/%Y"),
         tau.date = as.Date(tau_date, "%m/%d/%Y"),
         amy.tau.date.diff = abs(as.numeric(difftime(amy.date, tau.date, unit="weeks"))/52.25))%>%
  rename(race=Race, ID=subject,
         "bankssts"=SSTSBANK,
         "caudal anterior cingulate"=CAUDANTCNG,
         "caudal middle frontal"=CAUDMIDFRN,
         #"corpus callosum"=av1451_fsuvr_rsf_tot_ctx_crpclm,
         "cuneus"=CUNEUS,
         "entorhinal"=ENTORHINAL,
         "frontal pole"=FRNPOLE,
         "fusiform"=FUSIFORM,
         "inferior parietal"=INFERPRTL,
         "inferior temporal"=INFERTMP,
         "insula"=INSULA,
         "isthmus cingulate"=ISTHMUSCNG,
         "lateral occipital"=LATOCC,
         "lateral orbitofrontal"=LATORBFRN,
         "lingual"=LINGUAL,
         "medial orbitofrontal"=MEDORBFRN,
         "middle temporal"=MIDTMP,
         "paracentral"=PARACNTRL,
         "parahippocampal"=PARAHPCMPL,
         "pars opercularis"=PARSOPCLRS,
         "pars orbitalis"=PARSORBLS,
         "pars triangularis"=PARSTRNGLS,
         "pericalcarine"=PERICLCRN,
         "postcentral"=POSTCNTRL,
         "posterior cingulate"=POSTCNG,
         "precentral"=PRECNTRL,
         "precuneus"=PRECUNEUS,
         "rostral anterior cingulate"=ROSANTCNG,
         "rostral middle frontal"=ROSMIDFRN,
         "superior frontal"=SUPERFRN,
         "superior parietal"=SUPERPRTL,
         "superior temporal"=SUPERTMP,
         "supramarginal"=SUPRAMRGNL,
         "temporal pole"=TMPPOLE,
         "transverse temporal"=TRANSTMP,
         "amygdala"=AMYGDALA,
         "caudate"=CAUD,
         "hippocampus" =HIPPOCAMPUS,
         "pallidum"=PALLIDUM,
         "putamen"=PUTAMEN,
         "thalamus proper"=THALAMUS_PRPR,
         "ventral DC"=VENTRALDC,
         centiloid=centiloids)%>%
  dplyr::select(ID, cohort, age, sex, Female, race, race.binary,tau.date,
                "bankssts","caudal anterior cingulate","caudal middle frontal","cuneus","entorhinal",
                "frontal pole","fusiform","inferior parietal","inferior temporal","insula","isthmus cingulate","lateral occipital",
                "lateral orbitofrontal","lingual","medial orbitofrontal","middle temporal","paracentral","parahippocampal",
                "pars opercularis","pars orbitalis","pars triangularis","pericalcarine","postcentral","posterior cingulate",
                "precentral","precuneus","rostral anterior cingulate","rostral middle frontal","superior frontal",
                "superior parietal","superior temporal","supramarginal","temporal pole","transverse temporal","amygdala",
                "caudate","hippocampus","pallidum","putamen","thalamus proper","ventral DC",
                amy.date,Tracer,centiloid,amy_pos,amy.tau.date.diff, cdr.date, cdr, ANIMALS, tmb, digsym, srtfree, wais, memu, cdr.tau.date.diff)
brian_ADRC <- brian%>%filter(cohort == "LOAD") # n = 108
brian_DIAN <- brian%>%filter(cohort != "LOAD")%>% # n = 36
  dplyr::select(-cohort, -ANIMALS, -tmb, -wais, -memu)

# Cognitive tests
# DIAN OBS/TU: ANIMALS, TrailsB (TRAILB???), WAIS, MEMUNITS 
DIANID <- read.csv("./datafiles/df15_dianid_site_newid15.csv")
cog1 <- read.csv("./datafiles/psychometric_DCF.csv")%>% # n = 1338
  mutate(psy.date = as.Date(visit_date, "%m/%d/%Y"))%>%
  rename(tmb = TRAILB)%>%
  dplyr::select(newid15, visit, psy.date, ANIMALS, tmb) 
# cog1 <- read.csv("./datafiles/psychometric_DCF.csv")%>% # n = 1338
#   mutate(psy.date = as.Date(visit_date, "%m/%d/%Y"))%>%
#   rename(tmb = TRAILB, wais = WAIS, memu = MEMUNITS)%>%
#   dplyr::select(newid15, visit, psy.date, ANIMALS, tmb, wais, memu) 
dian_cog <- merge(DIANID, cog1, by="newid15")%>% # n = 1338
  rename(ID = DIANID)%>%
  dplyr::select(-newid15, -site, -visit)
dian_tu <- read.csv("./datafiles/DIAN_TU_FromShaney.csv")%>% # n = 27
  rename(tmb=TRAILB)%>%
  mutate(psy.date=NA)%>%
  dplyr::select(ID,psy.date,ANIMALS,tmb)
# dian_tu <- read.csv("./datafiles/DIAN_TU_FromShaney.csv")%>%
#   rename(tmb=TRAILB, wais=WAIS, memu=MEMUNITS)%>%
#   mutate(psy.date=NA)%>%
#   dplyr::select(ID,psy.date,ANIMALS,tmb,wais,memu)
dian1 <- rbind(dian_cog,dian_tu) # n = 1365

dian <- merge(brian_DIAN, dian1, by="ID")%>% # n = 35
  mutate(cdr.tau.date.diff = abs(as.numeric(difftime(psy.date, tau.date, unit="weeks")/52.25)))%>%
  group_by(ID)%>%
  distinct(tau.date, .keep_all = T)%>%
  ungroup()%>%
  dplyr::select(1:55, ANIMALS, tmb, cdr.tau.date.diff)
  # dplyr::select(1:55, ANIMALS, tmb, digsym, srtfree, wais, memu, cdr.tau.date.diff)

brianIDswithinFull <- as.data.frame(brian_ADRC$ID %in% cdr_tau$ID) # all of brian's ADRC IDs are in full ADRC, so all we need to merge are the DIAN Ids 
data <- rbind(cdr_tau, dian)%>% # n = 613 
  mutate(amy.cdr.date.diff = as.numeric(difftime(amy.date, cdr.date, unit="weeks"))/52.25)%>%
  mutate(group.cross.sectional = case_when(cdr == 0 & amy_pos == "A-" ~ "Control",
                                           cdr == 0 & amy_pos == "A+" ~ "Preclinical",
                                           cdr == 1 & amy_pos == "A-" ~ "Cognitively Impaired, A-",
                                           cdr == 1 & amy_pos == "A+" ~ "Symptomatic"))

filter_amy.tau <- data%>%filter(between(amy.tau.date.diff,-3,3)) # n = 444
filter_cdr.tau <- data%>%filter(between(cdr.tau.date.diff,-3,3)) # n = 363 (EXCLUDES ANYONE WITHOUT CDR INFO, IE DIAN!!~!!)
filter_amy.cdr <- data%>%filter(between(amy.cdr.date.diff,-3,3)) # n = 347 (EXCLUDES ANYONE WITHOUT CDR INFO, IE DIAN!!~!!)

# keep only baseline 
bl_data <- data%>% # n = 527
  group_by(ID)%>%
  arrange(tau.date)%>%
  filter(row_number()==1)%>%
  ungroup()

bl_filter_amy.tau <- filter_amy.tau%>% # n = 419
  group_by(ID)%>%
  arrange(tau.date)%>%
  filter(row_number()==1)%>%
  ungroup()
bl_filter_cdr.tau <- filter_cdr.tau%>% # n = 341 (EXCLUDES ANYONE WITHOUT CDR INFO!!~!!)
  group_by(ID)%>%
  arrange(tau.date)%>%
  filter(row_number()==1)%>%
  ungroup()
bl_filter_amy.cdr <- filter_amy.cdr%>% # n = 309 (EXCLUDES ANYONE WITHOUT CDR INFO!!~!!)
  group_by(ID)%>%
  arrange(tau.date)%>%
  filter(row_number()==1)%>%
  ungroup()

filter_all <-filter_amy.tau%>%filter(between(amy.cdr.date.diff,-3,3)) # (EXCLUDES ANYONE WITHOUT CDR INFO!!~!!)
filter_all <-filter_all%>%filter(between(cdr.tau.date.diff,-3,3)) #n = 262
bl_filter_all <- filter_all%>% # n = 255
  group_by(ID)%>%
  arrange(tau.date)%>%
  filter(row_number()==1)%>%
  ungroup()


########### NEXT STEPS 
# write.csv(data,"./all_data_n824.csv",row.names = FALSE)
# write.csv(bl_data,"./baseline_data_n552.csv",row.names = FALSE)
# write.csv(filter_amy.tau,"./filter.amy.tau_n422.csv",row.names = FALSE)
write.csv(bl_filter_amy.tau,"./neuropsych_baseline.filter.amy.tau_n419.csv",row.names = FALSE)
# write.csv(filter_cdr.tau,"./filter.cdr.tau_n483.csv",row.names = FALSE)
# write.csv(bl_filter_cdr.tau,"./baseline.filter.cdr.tau_n360.csv",row.names = FALSE)
# write.csv(filter_amy.cdr,"./filter.amy.cdr_n422.csv",row.names = FALSE)
# write.csv(bl_filter_amy.cdr,"./baseline.filter.amy.cdr_n304.csv",row.names = FALSE)
# write.csv(filter_all,"./filter_all_n294.csv",row.names = FALSE)
# write.csv(bl_filter_all,"./baseline.filter_all_n256.csv",row.names = FALSE)
# 


