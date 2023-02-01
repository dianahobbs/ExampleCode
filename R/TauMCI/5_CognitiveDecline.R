pacman::p_load(tidyverse,jtools,naniar,lubridate,janitor,lme4,sfsmisc,foreign,MASS,scales,ggseg,cowplot,gtsummary,ggplot2,RColorBrewer)
pacman::p_load(tidyverse,jtools,naniar,lubridate,janitor,lme4,sfsmisc,foreign,MASS,scales,ggseg,cowplot,gtsummary,ggplot2,RColorBrewer,gtsummary,webshot)
webshot::install_phantomjs(force=T)

#Load data file
df <- read.csv("./Data/N151_Wide_5-25-22.csv")%>%
  select(ID,Tau_Date,AGE,SEX,EDUC,Converter_after_tauPET,Biomarker_group_metaROI,Biomarker_group_EC,Biomarker_group_IT,Biomarker_group_ANY)
df <- df%>%filter(Converter_after_tauPET==0)

YrBl_Coeff<-read.csv("./5_25_22/long2_coefficients_PACC_YrBl_6-14.csv")%>%
  rename(ID=X)

Decliners<-YrBl_Coeff%>%mutate(Decliner = case_when(
  YrBl <= (-0.0606043) ~ "1", 
  YrBl > (-0.0606043) ~ "0"))

table(Decliners$Decliner)

BaseCog_Wide_Decline<-merge(df,Decliners,by="ID")

split_Decliners <- BaseCog_Wide_Decline%>%
  select(Biomarker_group_metaROI, Decliner)%>%
  tbl_summary(by=Biomarker_group_metaROI,
              statistic=list(all_continuous()~ "{mean}({sd})",
                             all_categorical()~"{n}/{N}({p}%)"),
              digits = all_continuous()~2,
              label = c(Decliner = "Decliners"))%>%
  as_gt()%>%             
  gt::gtsave(filename = "./6_15/split_YrBl_Decliners_6-15.png")
split_Decliners


##### Grab only nonconverting A-T- individuals to get cutorr & Determine who is declining v stable in nonconverters 
# metaROI

# ANY
YrBl_Coeff_ANY<-read.csv("./5_25_22/long2_coefficients_PACC_YrBl_ANY_6-14.csv")%>%
  rename(ID=X)

Decliners_ANY<-YrBl_Coeff_ANY%>%mutate(Decliner = case_when(
  YrBl <= (-0.0606043) ~ "1", 
  YrBl > (-0.0606043) ~ "0"))


YrBl_Data_ANY<-merge(df, Decliners_ANY, by="ID")

split_YrBl_Decliners_ANY <- YrBl_Data_ANY%>%
  select(Biomarker_group_ANY, Decliner)%>%
  tbl_summary(by=Biomarker_group_ANY,
              statistic=list(all_continuous()~ "{mean}({sd})",
                             all_categorical()~"{n}/{N}({p}%)"),
              digits = all_continuous()~2,
              label = c(Decliner = "Decliners"))%>%
  as_gt()%>%             
  gt::gtsave(filename = "./6_15/split_YrBl_Decliners_ANY_6-15.png")
split_YrBl_Decliners_ANY

#EC
YrBl_Coeff_EC<-read.csv("./5_25_22/long2_coefficients_PACC_YrBl_EC_6-14.csv")%>%
  rename(ID=X)

Decliners_EC<-YrBl_Coeff_EC%>%mutate(Decliner = case_when(
  YrBl <= (-0.0606043) ~ "1", 
  YrBl > (-0.0606043) ~ "0"))
YrBl_Data_EC<-merge(df, Decliners_EC, by="ID")

split_YrBl_Decliners_EC <- YrBl_Data_EC%>%
  select(Biomarker_group_EC, Decliner)%>%
  tbl_summary(by=Biomarker_group_EC,
              statistic=list(all_continuous()~ "{mean}({sd})",
                             all_categorical()~"{n}/{N}({p}%)"),
              digits = all_continuous()~2,
              label = c(Decliner = "Decliners"))%>%
  as_gt()%>%             
  gt::gtsave(filename = "./6_15/split_YrBl_Decliners_EC_6-15.png")
split_YrBl_Decliners_EC


#IT
YrBl_Coeff_IT<-read.csv("./5_25_22/long2_coefficients_PACC_YrBl_IT_6-14.csv")%>%
  rename(ID=X)

Decliners_IT<-YrBl_Coeff_IT%>%mutate(Decliner = case_when(
  YrBl <= (-0.0606043) ~ "1", 
  YrBl > (-0.0606043) ~ "0"))

YrBl_Data_IT<-merge(df, Decliners_IT, by="ID")

split_YrBl_Decliners_IT <- YrBl_Data_IT%>%
  select(Biomarker_group_IT, Decliner)%>%
  tbl_summary(by=Biomarker_group_IT,
              statistic=list(all_continuous()~ "{mean}({sd})",
                             all_categorical()~"{n}/{N}({p}%)"),
              digits = all_continuous()~2,
              label = c(Decliner = "Decliners"))%>%
  as_gt()%>%             
  gt::gtsave(filename = "./6_15/split_YrBl_Decliners_IT_6-15.png")
split_YrBl_Decliners_IT

# ANY


