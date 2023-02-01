# Title: ADRC Tau Imaging Setup
# Author: Diana Hobbs
# Date: April 2022

#### Set Up ####
# Load packages 
pacman::p_load(tidyverse,jtools,naniar,lubridate,BiocManager,purrr)

###### Load data and select columns of interest

# tau !*!*!!*!*! First open in excel .xlsx & format PET date to 3/14/2012 and THEN save to .csv !*!*!!*!*!
tau <- read.csv("./Data/tau.csv")%>%
  mutate(tau_Date = as.Date(PET_Date, "%m/%d/%Y"))%>%
  mutate(Tau_Summary=av1451_fsuvr_rsf_tot_amygdala+
           av1451_fsuvr_rsf_tot_ctx_entorhi+
           av1451_fsuvr_rsf_tot_ctx_infertm+
           av1451_fsuvr_rsf_tot_ctx_latocc)%>%
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
  dplyr::select(ID,tau_Date,"bankssts","caudal anterior cingulate","caudal middle frontal",
       "cuneus","entorhinal","frontal pole","fusiform","inferior parietal","inferior temporal",
       "insula","isthmus cingulate","lateral occipital","lateral orbitofrontal","lingual",
       "medial orbitofrontal","middle temporal","paracentral","parahippocampal","pars opercularis",
       "pars orbitalis","pars triangularis","pericalcarine","postcentral","posterior cingulate",
       "precentral","precuneus","rostral anterior cingulate","rostral middle frontal","superior frontal",
       "superior parietal","superior temporal","supramarginal","temporal pole","transverse temporal",
       "amygdala","caudate","hippocampus","pallidum","putamen","thalamus proper","ventral DC",
       Tau_Summary,Tauopathy)%>%
  na.omit()

# MRI T3 !*!*!!*!*! First open in excel .xlsx & format MR date to 3/14/2012 and THEN save to .csv !*!*!!*!*!
mri <- read.csv("./Data/mri_t3.csv", header=TRUE)%>%
  mutate(MR_Date = as.Date(MR_Date, "%m/%d/%Y"))%>%
  mutate("bankssts"=MR_LV_SSTSBANK+MR_RV_SSTSBANK)%>%
  mutate("caudal anterior cingulate"=MR_LV_CAUDANTCNG+MR_RV_CAUDANTCNG)%>%
  mutate("caudal middle frontal"=MR_LV_CAUDMIDFRN+MR_RV_CAUDMIDFRN)%>%
  mutate("cuneus"=MR_LV_CUNEUS+MR_RV_CUNEUS)%>%
  mutate("entorhinal"=MR_LV_ENTORHINAL+MR_RV_ENTORHINAL)%>%
  mutate("frontal pole"=MR_LV_FRNPOLE+MR_RV_FRNPOLE)%>%
  mutate("fusiform"=MR_LV_FUSIFORM+MR_RV_FUSIFORM)%>%
  mutate("inferior parietal"=MR_LV_INFRPRTL+MR_RV_INFRPRTL)%>%
  mutate("inferior temporal"=MR_LV_INFRTMP+MR_RV_INFRTMP)%>%
  mutate("insula"=MR_LV_INSULA+MR_RV_INSULA)%>%
  mutate("isthmus cingulate"=MR_LV_ISTHMUSCNG+MR_RV_ISTHMUSCNG)%>%
  mutate("lateral occipital"=MR_LV_LATOCC+MR_RV_LATOCC)%>%
  mutate("lateral orbitofrontal"=MR_LV_LATORBFRN+MR_RV_LATORBFRN)%>%
  mutate("lingual"=MR_LV_LINGUAL+MR_RV_LINGUAL)%>%
  mutate("medial orbitofrontal"=MR_LV_MEDORBFRN+MR_RV_MEDORBFRN)%>%
  mutate("middle temporal"=MR_LV_MIDTMP+MR_RV_MIDTMP)%>%
  mutate("paracentral"=MR_LV_PARACNTRL+MR_RV_PARACNTRL)%>%
  mutate("parahippocampal"=MR_LV_PARAHPCMPL+MR_RV_PARAHPCMPL)%>%
  mutate("pars opercularis"=MR_LV_PARAOPRCLRS+MR_RV_PARAOPRCLRS)%>%
  mutate("pars orbitalis"=MR_LV_PARSORBLS+MR_RV_PARSORBLS)%>%
  mutate("pars triangularis"=MR_LV_PARSTRNGLRS+MR_RV_PARSTRNGLRS)%>%
  mutate("pericalcarine"=MR_LV_PERICLCRN+MR_RV_PERICLCRN)%>%
  mutate("postcentral"=MR_LV_POSTCNG+MR_RV_POSTCNG)%>%
  mutate("posterior cingulate"=MR_LV_POSTCNTRL+MR_RV_POSTCNTRL)%>%
  mutate("precentral"=MR_LV_PRECNTRL+MR_RV_PRECNTRL)%>%
  mutate("precuneus"=MR_LV_PRECUNEUS+MR_RV_PRECUNEUS)%>%
  mutate("rostral anterior cingulate"=MR_LV_ROSANTCNG+MR_RV_ROSANTCNG)%>%
  mutate("rostral middle frontal"=MR_LV_ROSMIDFRN+MR_RV_ROSMIDFRN)%>%
  mutate("superior frontal"=MR_LV_SUPERFRN+MR_RV_SUPERFRN)%>%
  mutate("superior parietal"=MR_LV_SUPERPRTL+MR_RV_SUPERPRTL)%>%
  mutate("superior temporal"=MR_LV_SUPERTMP+MR_RV_SUPERTMP)%>%
  mutate("supramarginal"=MR_LV_SUPRAMRGNL+MR_RV_SUPRAMRGNL)%>%
  mutate("temporal pole"=MR_LV_TMPPOLE+MR_RV_TMPPOLE)%>%
  mutate("transverse temporal"=MR_LV_TRANSTMP+MR_RV_TRANSTMP)%>%
  mutate("amygdala"=MR_LV_AMYGDALA+MR_RV_AMYGDALA)%>%
  mutate("caudate"=MR_LV_CAUD+MR_RV_CAUD)%>%
  mutate("hippocampus"=MR_LV_HIPPOCAMPUS+MR_RV_HIPPOCAMPUS)%>%
  mutate("pallidum"=MR_LV_PALLIDUM+MR_RV_PALLIDUM)%>%
  mutate("putamen"=MR_LV_PUTAMEN+MR_RV_PUTAMEN)%>%
  mutate("thalamus proper"=MR_LV_THALAMUS+MR_RV_THALAMUS)%>%
  mutate("ventral DC"=MR_LV_VENTRALDC+MR_RV_VENTRALDC)%>%
  rename(ICV=MR_TOTV_INTRACRANIAL, WMH=MR_TOTV_WMHYPOINTENSITIES)%>%
  dplyr::select(ID, MR_Date, Field_Strength,
         "bankssts","caudal anterior cingulate","caudal middle frontal","cuneus","entorhinal",
         "frontal pole","fusiform","inferior parietal","inferior temporal","insula","isthmus cingulate",
         "lateral occipital","lateral orbitofrontal","lingual","medial orbitofrontal","middle temporal",
         "paracentral","parahippocampal","pars opercularis","pars orbitalis","pars triangularis",
         "pericalcarine","postcentral","posterior cingulate","precentral","precuneus","rostral anterior cingulate",
         "rostral middle frontal","superior frontal","superior parietal","superior temporal","supramarginal",
         "temporal pole","transverse temporal","amygdala","caudate","hippocampus","pallidum","putamen","thalamus proper",
         "ventral DC",WMH,ICV)%>%
  na.omit()
  
# PiB !*!*!!*!*! First open in excel .xlsx & format PET date to 3/14/2012 and THEN save to .csv !*!*!!*!*!
PiB <- read.csv("./Data/pib.csv")%>%
  mutate(PET_Date = as.Date(PET_Date, "%m/%d/%Y"))%>%
  rename("bankssts"=pib_fsuvr_rsf_tot_ctx_sstsbank,
         "caudal anterior cingulate"=pib_fsuvr_rsf_tot_ctx_caudantcng,
         "caudal middle frontal"=pib_fsuvr_rsf_tot_ctx_caudmidfrn,
         "cuneus"=pib_fsuvr_rsf_tot_ctx_cuneus,
         "entorhinal"=pib_fsuvr_rsf_tot_ctx_entorhinal,
         "frontal pole"=pib_fsuvr_rsf_tot_ctx_frnpole,
         "fusiform"=pib_fsuvr_rsf_tot_ctx_fusiform,
         "inferior parietal"=pib_fsuvr_rsf_tot_ctx_inferprtl,
         "inferior temporal"=pib_fsuvr_rsf_tot_ctx_infertmp,
         "insula"=pib_fsuvr_rsf_tot_ctx_insula,
         "isthmus cingulate"=pib_fsuvr_rsf_tot_ctx_isthmuscng,
         "lateral occipital"=pib_fsuvr_rsf_tot_ctx_latocc,
         "lateral orbitofrontal"=pib_fsuvr_rsf_tot_ctx_latorbfrn,
         "lingual"=pib_fsuvr_rsf_tot_ctx_lingual,
         "medial orbitofrontal"=pib_fsuvr_rsf_tot_ctx_medorbfrn,
         "middle temporal"=pib_fsuvr_rsf_tot_ctx_midtmp,
         "paracentral"=pib_fsuvr_rsf_tot_ctx_paracntrl,
         "parahippocampal"=pib_fsuvr_rsf_tot_ctx_parahpcmpl,
         "pars opercularis"=pib_fsuvr_rsf_tot_ctx_parsopclrs,
         "pars orbitalis"=pib_fsuvr_rsf_tot_ctx_parsorbls,
         "pars triangularis"=pib_fsuvr_rsf_tot_ctx_parstrngls,
         "pericalcarine"=pib_fsuvr_rsf_tot_ctx_periclcrn,
         "postcentral"=pib_fsuvr_rsf_tot_ctx_postcntrl,
         "posterior cingulate"=pib_fsuvr_rsf_tot_ctx_postcng,
         "precentral"=pib_fsuvr_rsf_tot_ctx_precntrl,
         "precuneus"=pib_fsuvr_rsf_tot_ctx_precuneus,
         "rostral anterior cingulate"=pib_fsuvr_rsf_tot_ctx_rosantcng,
         "rostral middle frontal"=pib_fsuvr_rsf_tot_ctx_rosmidfrn,
         "superior frontal"=pib_fsuvr_rsf_tot_ctx_superfrn,
         "superior parietal"=pib_fsuvr_rsf_tot_ctx_superprtl,
         "superior temporal"=pib_fsuvr_rsf_tot_ctx_supertmp,
         "supramarginal"=pib_fsuvr_rsf_tot_ctx_supramrgnl,
         "temporal pole"=pib_fsuvr_rsf_tot_ctx_tmppole,
         "transverse temporal"=pib_fsuvr_rsf_tot_ctx_transtmp,
         "amygdala"=pib_fsuvr_rsf_tot_amygdala,
         "caudate"=pib_fsuvr_rsf_tot_caud,
         "hippocampus" =pib_fsuvr_rsf_tot_hippocampus,
         "pallidum"=pib_fsuvr_rsf_tot_pallidum,
         "putamen"=pib_fsuvr_rsf_tot_putamen,
         "thalamus proper"=pib_fsuvr_rsf_tot_thalamus_prpr,
         "ventral DC"=pib_fsuvr_rsf_tot_ventraldc)%>%
 dplyr::select(ID, PET_Date,"bankssts","caudal anterior cingulate","caudal middle frontal",
                "cuneus","entorhinal","frontal pole","fusiform","inferior parietal","inferior temporal",
                "insula","isthmus cingulate","lateral occipital","lateral orbitofrontal","lingual",
                "medial orbitofrontal","middle temporal","paracentral","parahippocampal","pars opercularis",
                "pars orbitalis","pars triangularis","pericalcarine","postcentral","posterior cingulate",
                "precentral","precuneus","rostral anterior cingulate","rostral middle frontal","superior frontal",
                "superior parietal","superior temporal","supramarginal","temporal pole","transverse temporal",
                "amygdala","caudate","hippocampus","pallidum","putamen","thalamus proper","ventral DC")%>%  
  na.omit()
