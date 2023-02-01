# Title: ADRC Cohort Tau/Age Regional Analyses for Brian
# Author: Diana Hobbs, edited code from Brian 
# Date: November 2022

#### Set Up ####
# Load packages 
pacman::p_load(tidyverse,readxl,writexl,ggplot2,scales,cowplot,ggseg,MASS,sfsmisc,robustbase,lme4)
#clear memory
#READ IN DATA from original ADRC_Setup.csv file "./baseline.filter.amy.tau_n420.csv"
sample <- read.csv("./neuropsych_baseline.filter.amy.tau_n419.csv")%>%
   mutate(cohort = case_when(is.na(cdr.date) ~ "DIAN",
                             !is.na(cdr.date) ~ "ADRC"))

tauage <- sample%>%filter(amy_pos == "A-") # n = 279

tauage_DIAN <- tauage%>%filter(cohort=="DIAN") # n = 34
tauage_ADRC <- tauage%>%filter(cohort=="ADRC") # n = 245
range(tauage_DIAN$age) # 18 to 66
range(tauage_ADRC$age) # 45 to 91

ggseg_list = read.csv("./ggseg_namelist.csv")

# reassign variable types
tauage$sex=as.factor(tauage$sex); tauage$Female=as.factor(tauage$Female); tauage$race.binary=as.factor(tauage$race.binary); 
tauage$sex <- relevel(tauage$sex, "Male") # set reference as we expect women to have greater levels


###### REGRESSION LOOP #####
roi <- tauage[,8:48]

# # Loop through regressions and keep output
region_fit.COF<-list() # output for regression
robust_fit.COF<-list() # output for robust regression
region_fit_int.COF<-list() 
robust_fit_int.COF<-list() 
model_robust<- matrix(0, 41, 5)
model_robust_int<- matrix(0, 41, 6)
outcount=1
# Create a 'for' loop 
for (i in 1:ncol(roi)){
  region_fit <- lm(roi[,i] ~ sex+centiloid+race.binary+age, data=tauage)
  robust_fit <- rlm(roi[,i] ~ sex+centiloid+race.binary+age, data=tauage)
  region_fit_int <- lm(roi[,i] ~ sex+centiloid+race.binary+age+age*sex, data=tauage)
  robust_fit_int <- rlm(roi[,i] ~ sex+centiloid+race.binary+age+age*sex, data=tauage)
  
  region_fit.sum <- summary(region_fit)
  robust_fit.sum <- summary(robust_fit)  
  region_fit_int.sum <- summary(region_fit_int)
  robust_fit_int.sum <- summary(robust_fit_int)
  
  region_fit.cof <- region_fit.sum$coefficients
  robust_fit.cof <- robust_fit.sum$coefficients
  region_fit_int.cof <- region_fit_int.sum$coefficients
  robust_fit_int.cof <- robust_fit_int.sum$coefficients
  
  region=colnames(roi)[i]
  model_robust[outcount,1]=region
  model_robust_int[outcount,1]=region
  colnames(model_robust) <- c("region", "p.adjust.sex", "p.adjust.centiloid", "p.adjust.race", "p.adjust.age")
  colnames(model_robust_int) <- c("region", "p.adjust.sex", "p.adjust.centiloid", "p.adjust.race", "p.adjust.age", "p.adjust.sex:age")
  
  p_robust=f.robftest(robust_fit, var = "sexFemale"); model_robust[outcount, 2]=p_robust$p.value 
  p_robust=f.robftest(robust_fit, var = "centiloid"); model_robust[outcount, 3]=p_robust$p.value
  p_robust=f.robftest(robust_fit, var = "race.binary1");model_robust[outcount, 4]=p_robust$p.value 
  p_robust=f.robftest(robust_fit, var = "age");model_robust[outcount, 5]=p_robust$p.value 
  
  p_robust_int=f.robftest(robust_fit_int, var = "sexFemale"); model_robust_int[outcount, 2]=p_robust_int$p.value 
  p_robust_int=f.robftest(robust_fit_int, var = "centiloid"); model_robust_int[outcount, 3]=p_robust_int$p.value
  p_robust_int=f.robftest(robust_fit_int, var = "race.binary1");model_robust_int[outcount, 4]=p_robust_int$p.value 
  p_robust_int=f.robftest(robust_fit_int, var = "age");model_robust_int[outcount, 5]=p_robust_int$p.value 
  p_robust_int=f.robftest(robust_fit_int, var = "sexFemale:age");model_robust_int[outcount, 6]=p_robust_int$p.value
  
  outcount=outcount+1
  region_fit.COF[[i]] <- region_fit.cof
  robust_fit.COF[[i]] <- robust_fit.cof
  region_fit_int.COF[[i]] <- region_fit_int.cof
  robust_fit_int.COF[[i]] <- robust_fit_int.cof
  
  names(region_fit.COF)[i] <- colnames(roi)[i]
  names(robust_fit.COF)[i] <- colnames(roi)[i]
  names(region_fit_int.COF)[i] <- colnames(roi)[i]
  names(robust_fit_int.COF)[i] <- colnames(roi)[i]
}
rm(region_fit,region_fit.cof,region_fit.sum,robust_fit,robust_fit.cof,robust_fit.sum,
   region_fit_int,region_fit_int.cof,region_fit_int.sum,robust_fit_int,robust_fit_int.cof,robust_fit_int.sum,
   i,roi,p_robust,p_robust_int,outcount,region)

region_fit.COF <- as.data.frame(region_fit.COF)
robust_fit.COF <- as.data.frame(robust_fit.COF)
region_fit_int.COF <- as.data.frame(region_fit_int.COF)
robust_fit_int.COF <- as.data.frame(robust_fit_int.COF)

##### Extract Estimates, T-value, P-value #####
# region_fit.COF
Estimate_region <- as.data.frame(t(region_fit.COF[,grepl(".Estimate", names(region_fit.COF))]))
TVal_region <- as.data.frame(t(region_fit.COF[,grepl(".t.value", names(region_fit.COF))]))
PVal_region <- as.data.frame(t(region_fit.COF[,grepl(".Pr...t..", names(region_fit.COF))]))
# robust_fit.COF
Estimate_robust <- as.data.frame(t(robust_fit.COF[,grepl(".Value", names(robust_fit.COF))]))
TVal_robust <- as.data.frame(t(robust_fit.COF[,grepl(".t.value", names(robust_fit.COF))]))
PVal_robust <- as.data.frame(model_robust)
colnames(PVal_robust) <- c("(Intercept)","sexFemale","centiloid","race.binary1","age")
# region_fit_int.COF
Estimate_region_int <- as.data.frame(t(region_fit_int.COF[,grepl(".Estimate", names(region_fit_int.COF))]))
TVal_region_int <- as.data.frame(t(region_fit_int.COF[,grepl(".t.value", names(region_fit_int.COF))]))
PVal_region_int <- as.data.frame(t(region_fit_int.COF[,grepl(".Pr...t..", names(region_fit_int.COF))]))
# robust_fit_int.COF
Estimate_robust_int <- as.data.frame(t(robust_fit_int.COF[,grepl(".Value", names(robust_fit_int.COF))]))
TVal_robust_int <- as.data.frame(t(robust_fit_int.COF[,grepl(".t.value", names(robust_fit_int.COF))]))
PVal_robust_int <- as.data.frame(model_robust_int)
colnames(PVal_robust_int) <- c("(Intercept)","sexFemale","centiloid","race.binary1","age","sexFemale:age")

##### Merge coefficient summaries & perform FDR correction using Benjamini and Hochberg method #####
# region
sex_region <- data.frame("B"=Estimate_region$sexFemale,"t"=TVal_region$sexFemale, "p"=PVal_region$sexFemale,"p_corr"=(p.adjust(PVal_region$sexFemale, method = "BH",n = length(PVal_region$sexFemale))))
centiloid_region <- data.frame("B"=Estimate_region$centiloid,"t"=TVal_region$centiloid,"p"=PVal_region$centiloid,"p_corr"=(p.adjust(PVal_region$centiloid, method = "BH",n = length(PVal_region$centiloid))))
race_region <- data.frame("B"=Estimate_region$race.binary1,"t"=TVal_region$race.binary1,"p"=PVal_region$race.binary1,"p_corr"=(p.adjust(PVal_region$race.binary1, method = "BH",n = length(PVal_region$race.binary1))))
age_region <-  data.frame("B"=Estimate_region$age,"t"=TVal_region$age,"p"=PVal_region$age,"p_corr"=(p.adjust(PVal_region$age, method = "BH",n = length(PVal_region$age))))
# robust
sex_robust <- data.frame("B"=Estimate_robust$sexFemale,"t"=TVal_robust$sexFemale, "p"=PVal_robust$sexFemale,"p_corr"=(p.adjust(PVal_robust$sexFemale, method = "BH",n = length(PVal_robust$sexFemale))))
centiloid_robust <- data.frame("B"=Estimate_robust$centiloid,"t"=TVal_robust$centiloid,"p"=PVal_robust$centiloid,"p_corr"=(p.adjust(PVal_robust$centiloid, method = "BH",n = length(PVal_robust$centiloid))))
race_robust <- data.frame("B"=Estimate_robust$race.binary1,"t"=TVal_robust$race.binary1,"p"=PVal_robust$race.binary1,"p_corr"=(p.adjust(PVal_robust$race.binary1, method = "BH",n = length(PVal_robust$race.binary1))))
age_robust <-  data.frame("B"=Estimate_robust$age,"t"=TVal_robust$age,"p"=PVal_robust$age,"p_corr"=(p.adjust(PVal_robust$age, method = "BH",n = length(PVal_robust$age))))
# region_int
sex_region_int <- data.frame("B"=Estimate_region_int$sexFemale,"t"=TVal_region_int$sexFemale, "p"=PVal_region_int$sexFemale,"p_corr"=(p.adjust(PVal_region_int$sexFemale, method = "BH",n = length(PVal_region_int$sexFemale))))
centiloid_region_int <- data.frame("B"=Estimate_region_int$centiloid,"t"=TVal_region_int$centiloid,"p"=PVal_region_int$centiloid,"p_corr"=(p.adjust(PVal_region_int$centiloid, method = "BH",n = length(PVal_region_int$centiloid))))
race_region_int <- data.frame("B"=Estimate_region_int$race.binary1,"t"=TVal_region_int$race.binary1,"p"=PVal_region_int$race.binary1,"p_corr"=(p.adjust(PVal_region_int$race.binary1, method = "BH",n = length(PVal_region_int$race.binary1))))
age_region_int <-  data.frame("B"=Estimate_region_int$age,"t"=TVal_region_int$age,"p"=PVal_region_int$age,"p_corr"=(p.adjust(PVal_region_int$age, method = "BH",n = length(PVal_region_int$age))))
sex_age_region_int <-  data.frame("B"=Estimate_region_int$"sexFemale:age","t"=TVal_region_int$"sexFemale:age","p"=PVal_region_int$"sexFemale:age","p_corr"=(p.adjust(PVal_region_int$"sexFemale:age", method = "BH",n = length(PVal_region_int$"sexFemale:age"))))
# robust_int 
sex_robust_int <- data.frame("B"=Estimate_robust_int$sexFemale,"t"=TVal_robust_int$sexFemale, "p"=PVal_robust_int$sexFemale,"p_corr"=(p.adjust(PVal_robust_int$sexFemale, method = "BH",n = length(PVal_robust_int$sexFemale))))
centiloid_robust_int <- data.frame("B"=Estimate_robust_int$centiloid,"t"=TVal_robust_int$centiloid,"p"=PVal_robust_int$centiloid,"p_corr"=(p.adjust(PVal_robust_int$centiloid, method = "BH",n = length(PVal_robust_int$centiloid))))
race_robust_int <- data.frame("B"=Estimate_robust_int$race.binary1,"t"=TVal_robust_int$race.binary1,"p"=PVal_robust_int$race.binary1,"p_corr"=(p.adjust(PVal_robust_int$race.binary1, method = "BH",n = length(PVal_robust_int$race.binary1))))
age_robust_int <-  data.frame("B"=Estimate_robust_int$age,"t"=TVal_robust_int$age,"p"=PVal_robust_int$age,"p_corr"=(p.adjust(PVal_robust_int$age, method = "BH",n = length(PVal_robust_int$age))))
sex_age_robust_int <-  data.frame("B"=Estimate_robust_int$"sexFemale:age","t"=TVal_robust_int$"sexFemale:age","p"=PVal_robust_int$"sexFemale:age","p_corr"=(p.adjust(PVal_robust_int$"sexFemale:age", method = "BH",n = length(PVal_robust_int$"sexFemale:age"))))
rm(Estimate_region,TVal_region,PVal_region,Estimate_robust,TVal_robust,PVal_robust,
   Estimate_region_int,TVal_region_int,PVal_region_int,Estimate_robust_int,TVal_robust_int,PVal_robust_int)

##### Extract cortical and subcortical regions for plotting in ggseg #####
cort_list <- read.csv("./ggseg_namelist.csv")
subcort_list <- cort_list$subcortical[1:7] #extract subcortical list
# Sex
CortSex_region <- data.frame(region=cort_list$cortical,filter(sex_region[1:34,])); SubcortSex_region <- data.frame(region=subcort_list,filter(sex_region[35:41,]))
CortSex_robust <- data.frame(region=cort_list$cortical,filter(sex_robust[1:34,])); SubcortSex_robust <- data.frame(region=subcort_list,filter(sex_robust[35:41,]))
CortSex_region_int <- data.frame(region=cort_list$cortical,filter(sex_region_int[1:34,])); SubcortSex_region_int <- data.frame(region=subcort_list,filter(sex_region_int[35:41,]))
CortSex_robust_int <- data.frame(region=cort_list$cortical,filter(sex_robust_int[1:34,])); SubcortSex_robust_int <- data.frame(region=subcort_list,filter(sex_robust_int[35:41,]))
# Centiloid 
CortCentiloid_region <- data.frame(region=cort_list$cortical,filter(centiloid_region[1:34,])); SubcortCentiloid_region <- data.frame(region=subcort_list,filter(centiloid_region[35:41,]))
CortCentiloid_robust <- data.frame(region=cort_list$cortical,filter(centiloid_robust[1:34,])); SubcortCentiloid_robust <- data.frame(region=subcort_list,filter(centiloid_robust[35:41,]))
CortCentiloid_region_int <- data.frame(region=cort_list$cortical,filter(centiloid_region_int[1:34,])); SubcortCentiloid_region_int <- data.frame(region=subcort_list,filter(centiloid_region_int[35:41,]))
CortCentiloid_robust_int <- data.frame(region=cort_list$cortical,filter(centiloid_robust_int[1:34,])); SubcortCentiloid_robust_int <- data.frame(region=subcort_list,filter(centiloid_robust_int[35:41,]))
# Race
CortRace_region <- data.frame(region=cort_list$cortical,filter(race_region[1:34,])); SubcortRace_region <- data.frame(region=subcort_list,filter(race_region[35:41,]))
CortRace_robust <- data.frame(region=cort_list$cortical,filter(race_robust[1:34,])); SubcortRace_robust <- data.frame(region=subcort_list,filter(race_robust[35:41,]))
CortRace_region_int <- data.frame(region=cort_list$cortical,filter(race_region_int[1:34,])); SubcortRace_region_int <- data.frame(region=subcort_list,filter(race_region_int[35:41,]))
CortRace_robust_int <- data.frame(region=cort_list$cortical,filter(race_robust_int[1:34,])); SubcortRace_robust_int <- data.frame(region=subcort_list,filter(race_robust_int[35:41,]))
# Age
CortAge_region <- data.frame(region=cort_list$cortical,filter(age_region[1:34,])); SubcortAge_region <- data.frame(region=subcort_list,filter(age_region[35:41,]))
CortAge_robust <- data.frame(region=cort_list$cortical,filter(age_robust[1:34,])); SubcortAge_robust <- data.frame(region=subcort_list,filter(age_robust[35:41,]))
CortAge_region_int <- data.frame(region=cort_list$cortical,filter(age_region_int[1:34,])); SubcortAge_region_int <- data.frame(region=subcort_list,filter(age_region_int[35:41,]))
CortAge_robust_int <- data.frame(region=cort_list$cortical,filter(age_robust_int[1:34,])); SubcortAge_robust_int <- data.frame(region=subcort_list,filter(age_robust_int[35:41,]))
# Sex*Age
CortSexAge_region_int <- data.frame(region=cort_list$cortical,filter(sex_age_region_int[1:34,])); SubcortSexAge_region_int <- data.frame(region=subcort_list,filter(sex_age_region_int[35:41,]))
CortSexAge_robust_int <- data.frame(region=cort_list$cortical,filter(sex_age_robust_int[1:34,])); SubcortSexAge_robust_int <- data.frame(region=subcort_list,filter(sex_age_robust_int[35:41,]))
# Remove unneeded variables 
rm(sex_region, centiloid_region, race_region, age_region,sex_robust, centiloid_robust, race_robust, age_robust,
   sex_region_int, centiloid_region_int, race_region_int, age_region_int, sex_age_region_int, sex_robust_int, centiloid_robust_int, race_robust_int, age_robust_int, sex_age_robust_int)

##### Only keep those that are significant 
# region
CortSex_region.Sig <- subset(CortSex_region, p_corr < 0.05); SubcortSex_region.Sig <- subset(SubcortSex_region, p_corr < 0.05)
CortCentiloid_region.Sig <- subset(CortCentiloid_region, p_corr < 0.05); SubcortCentiloid_region.Sig <- subset(SubcortCentiloid_region, p_corr < 0.05)
CortRace_region.Sig <- subset(CortRace_region, p_corr < 0.05); SubcortRace_region.Sig <- subset(SubcortRace_region, p_corr < 0.05)
CortAge_region.Sig <- subset(CortAge_region, p_corr < 0.05); SubcortAge_region.Sig <- subset(SubcortAge_region, p_corr < 0.05)
# robust
CortSex_robust.Sig <- subset(CortSex_robust, p_corr < 0.05); SubcortSex_robust.Sig <- subset(SubcortSex_robust, p_corr < 0.05)
CortCentiloid_robust.Sig <- subset(CortCentiloid_robust, p_corr < 0.05); SubcortCentiloid_robust.Sig <- subset(SubcortCentiloid_robust, p_corr < 0.05)
CortRace_robust.Sig <- subset(CortRace_robust, p_corr < 0.05); SubcortRace_robust.Sig <- subset(SubcortRace_robust, p_corr < 0.05)
CortAge_robust.Sig <- subset(CortAge_robust, p_corr < 0.05); SubcortAge_robust.Sig <- subset(SubcortAge_robust, p_corr < 0.05)
# region_int
CortSex_region_int.Sig <- subset(CortSex_region_int, p_corr < 0.05); SubcortSex_region_int.Sig <- subset(SubcortSex_region_int, p_corr < 0.05)
CortCentiloid_region_int.Sig <- subset(CortCentiloid_region_int, p_corr < 0.05); SubcortCentiloid_region_int.Sig <- subset(SubcortCentiloid_region_int, p_corr < 0.05)
CortRace_region_int.Sig <- subset(CortRace_region_int, p_corr < 0.05); SubcortRace_region_int.Sig <- subset(SubcortRace_region_int, p_corr < 0.05)
CortAge_region_int.Sig <- subset(CortAge_region_int, p_corr < 0.05); SubcortAge_region_int.Sig <- subset(SubcortAge_region_int, p_corr < 0.05)
CortSexAge_region_int.Sig <- subset(CortSexAge_region_int, p_corr < 0.05); SubcortSexAge_region_int.Sig <- subset(SubcortSexAge_region_int, p_corr < 0.05)
# robust_int
CortSex_robust_int.Sig <- subset(CortSex_robust_int, p_corr < 0.05); SubcortSex_robust_int.Sig <- subset(SubcortSex_robust_int, p_corr < 0.05)
CortCentiloid_robust_int.Sig <- subset(CortCentiloid_robust_int, p_corr < 0.05); SubcortCentiloid_robust_int.Sig <- subset(SubcortCentiloid_robust_int, p_corr < 0.05)
CortRace_robust_int.Sig <- subset(CortRace_robust_int, p_corr < 0.05); SubcortRace_robust_int.Sig <- subset(SubcortRace_robust_int, p_corr < 0.05)
CortAge_robust_int.Sig <- subset(CortAge_robust_int, p_corr < 0.05); SubcortAge_robust_int.Sig <- subset(SubcortAge_robust_int, p_corr < 0.05)
CortSexAge_robust_int.Sig <- subset(CortSexAge_robust_int, p_corr < 0.05); SubcortSexAge_robust_int.Sig <- subset(SubcortSexAge_robust_int, p_corr < 0.05)

##### Make brain images #####
#This theme specifies parameters for all of the plots. Change here to change any of them
custom_theme = list(
  scale_fill_gradientn(colours=c("dodgerblue4","light blue","white","firebrick","goldenrod"), limits=c(-4,4), oob=squish),
  #use scale_fill_gradient2 and the option midpoint to specify white as another option
  theme(text = element_text(colour = "white")),  theme(axis.text = element_blank(), axis.title = element_blank()), 
  theme(plot.background = element_rect(fill = "black",colour = NA)), theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), theme(panel.border = element_blank()), guides(fill="none", color="none"))

custom_theme_sub = list(
  scale_fill_gradientn(colours=c("dodgerblue4","light blue","white","firebrick","goldenrod"), limits=c(-4,4), oob=squish),
  theme(text = element_text(colour = "white")),  theme(axis.text = element_blank(), axis.title = element_blank()), 
  theme(plot.background = element_rect(fill = "black",colour = NA)), theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), theme(panel.border = element_blank()))

##### Plotting t-values on the brain atlas
### Sex
# region
CortSex_region.t=ggseg(.data=CortSex_region,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+labs(title="Sex")+custom_theme
SubcortSex_region.t=ggseg(.data=SubcortSex_region,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+custom_theme_sub
# robust
CortSex_robust.t=ggseg(.data=CortSex_robust,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+labs(title="Sex")+custom_theme
SubcortSex_robust.t=ggseg(.data=SubcortSex_robust,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+custom_theme_sub
# region_int
CortSex_region_int.t=ggseg(.data=CortSex_region_int,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+labs(title="Sex")+custom_theme
SubcortSex_region_int.t=ggseg(.data=SubcortSex_region_int,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+custom_theme_sub
# robust_int
CortSex_robust_int.t=ggseg(.data=CortSex_robust_int,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+labs(title="Sex")+custom_theme
SubcortSex_robust_int.t=ggseg(.data=SubcortSex_robust_int,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+custom_theme_sub

### Centiloid
# region
CortCentiloid_region.t=ggseg(.data=CortCentiloid_region,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+labs(title="Centiloid")+custom_theme
SubcortCentiloid_region.t=ggseg(.data=SubcortCentiloid_region,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+custom_theme_sub
# robust
CortCentiloid_robust.t=ggseg(.data=CortCentiloid_robust,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+labs(title="Centiloid")+custom_theme
SubcortCentiloid_robust.t=ggseg(.data=SubcortCentiloid_robust,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+custom_theme_sub
# region_int
CortCentiloid_region_int.t=ggseg(.data=CortCentiloid_region_int,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+labs(title="Centiloid")+custom_theme
SubcortCentiloid_region_int.t=ggseg(.data=SubcortCentiloid_region_int,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+custom_theme_sub
# robust_int
CortCentiloid_robust_int.t=ggseg(.data=CortCentiloid_robust_int,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+labs(title="Centiloid")+custom_theme
SubcortCentiloid_robust_int.t=ggseg(.data=SubcortCentiloid_robust_int,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+custom_theme_sub

### Race
# region
CortRace_region.t=ggseg(.data=CortRace_region,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+labs(title="Race")+custom_theme
SubcortRace_region.t=ggseg(.data=SubcortRace_region,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+custom_theme_sub
# robust
CortRace_robust.t=ggseg(.data=CortRace_robust,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+labs(title="Race")+custom_theme
SubcortRace_robust.t=ggseg(.data=SubcortRace_robust,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+custom_theme_sub
# region_int
CortRace_region_int.t=ggseg(.data=CortRace_region_int,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+labs(title="Race")+custom_theme
SubcortRace_region_int.t=ggseg(.data=SubcortRace_region_int,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+custom_theme_sub
# robust_int
CortRace_robust_int.t=ggseg(.data=CortRace_robust_int,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+labs(title="Race")+custom_theme
SubcortRace_robust_int.t=ggseg(.data=SubcortRace_robust_int,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+custom_theme_sub

### Age
# region
CortAge_region.t=ggseg(.data=CortAge_region,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+labs(title="Age")+custom_theme
SubcortAge_region.t=ggseg(.data=SubcortAge_region,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+custom_theme_sub
# robust
CortAge_robust.t=ggseg(.data=CortAge_robust,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+labs(title="Age")+custom_theme
SubcortAge_robust.t=ggseg(.data=SubcortAge_robust,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+custom_theme_sub
# region_int
CortAge_region_int.t=ggseg(.data=CortAge_region_int,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+labs(title="Age")+custom_theme
SubcortAge_region_int.t=ggseg(.data=SubcortAge_region_int,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+custom_theme_sub
# robust_int
CortAge_robust_int.t=ggseg(.data=CortAge_robust_int,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+labs(title="Age")+custom_theme
SubcortAge_robust_int.t=ggseg(.data=SubcortAge_robust_int,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+custom_theme_sub

### Sex*Age
# region_int
CortSexAge_region_int.t=ggseg(.data=CortSexAge_region_int,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+labs(title="Sex*Age")+custom_theme
SubcortSexAge_region_int.t=ggseg(.data=SubcortSexAge_region_int,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+custom_theme_sub
# robust_int
CortSexAge_robust_int.t=ggseg(.data=CortSexAge_robust_int,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+labs(title="Sex*Age")+custom_theme
SubcortSexAge_robust_int.t=ggseg(.data=SubcortSexAge_robust_int,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+custom_theme_sub

##### Plotting only significant ones on the brain atlas
### Sex
# region
CortSex_region.Sig.t=ggseg(.data=CortSex_region.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+labs(title="Sex")+custom_theme
SubcortSex_region.Sig.t=ggseg(.data=SubcortSex_region.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+custom_theme_sub
# robust
CortSex_robust.Sig.t=ggseg(.data=CortSex_robust.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+labs(title="Sex")+custom_theme
SubcortSex_robust.Sig.t=ggseg(.data=SubcortSex_robust.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+custom_theme_sub
# region_int.Sig
CortSex_region_int.Sig.t=ggseg(.data=CortSex_region_int.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+labs(title="Sex")+custom_theme
SubcortSex_region_int.Sig.t=ggseg(.data=SubcortSex_region_int.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+custom_theme_sub
# robust_int.Sig
CortSex_robust_int.Sig.t=ggseg(.data=CortSex_robust_int.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+labs(title="Sex")+custom_theme
SubcortSex_robust_int.Sig.t=ggseg(.data=SubcortSex_robust_int.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+custom_theme_sub

### Centiloid
# region
CortCentiloid_region.Sig.t=ggseg(.data=CortCentiloid_region.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+labs(title="Centiloid")+custom_theme
SubcortCentiloid_region.Sig.t=ggseg(.data=SubcortCentiloid_region.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+custom_theme_sub
# robust
CortCentiloid_robust.Sig.t=ggseg(.data=CortCentiloid_robust.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+labs(title="Centiloid")+custom_theme
SubcortCentiloid_robust.Sig.t=ggseg(.data=SubcortCentiloid_robust.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+custom_theme_sub
# region_int.Sig
CortCentiloid_region_int.Sig.t=ggseg(.data=CortCentiloid_region_int.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+labs(title="Centiloid")+custom_theme
SubcortCentiloid_region_int.Sig.t=ggseg(.data=SubcortCentiloid_region_int.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+custom_theme_sub
# robust_int.Sig
CortCentiloid_robust_int.Sig.t=ggseg(.data=CortCentiloid_robust_int.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+labs(title="Centiloid")+custom_theme
SubcortCentiloid_robust_int.Sig.t=ggseg(.data=SubcortCentiloid_robust_int.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+custom_theme_sub

### Race
# region
CortRace_region.Sig.t=ggseg(.data=CortRace_region.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+labs(title="Race")+custom_theme
SubcortRace_region.Sig.t=ggseg(.data=SubcortRace_region.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+custom_theme_sub
# robust
CortRace_robust.Sig.t=ggseg(.data=CortRace_robust.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+labs(title="Race")+custom_theme
SubcortRace_robust.Sig.t=ggseg(.data=SubcortRace_robust.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+custom_theme_sub
# region_int.Sig
CortRace_region_int.Sig.t=ggseg(.data=CortRace_region_int.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+labs(title="Race")+custom_theme
SubcortRace_region_int.Sig.t=ggseg(.data=SubcortRace_region_int.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+custom_theme_sub
# robust_int.Sig
CortRace_robust_int.Sig.t=ggseg(.data=CortRace_robust_int.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+labs(title="Race")+custom_theme
SubcortRace_robust_int.Sig.t=ggseg(.data=SubcortRace_robust_int.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+custom_theme_sub

### Age
# region
CortAge_region.Sig.t=ggseg(.data=CortAge_region.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+labs(title="Age")+custom_theme
SubcortAge_region.Sig.t=ggseg(.data=SubcortAge_region.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+custom_theme_sub
# robust
CortAge_robust.Sig.t=ggseg(.data=CortAge_robust.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+labs(title="Age")+custom_theme
SubcortAge_robust.Sig.t=ggseg(.data=SubcortAge_robust.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+custom_theme_sub
# region_int.Sig
CortAge_region_int.Sig.t=ggseg(.data=CortAge_region_int.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+labs(title="Age")+custom_theme
SubcortAge_region_int.Sig.t=ggseg(.data=SubcortAge_region_int.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+custom_theme_sub
# robust_int.Sig
CortAge_robust_int.Sig.t=ggseg(.data=CortAge_robust_int.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+labs(title="Age")+custom_theme
SubcortAge_robust_int.Sig.t=ggseg(.data=SubcortAge_robust_int.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+custom_theme_sub

### Sex*Age
# region_int.Sig
CortSexAge_region_int.Sig.t=ggseg(.data=CortSexAge_region_int.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+labs(title="Sex*Age")+custom_theme
SubcortSexAge_region_int.Sig.t=ggseg(.data=SubcortSexAge_region_int.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+custom_theme_sub
# robust_int.Sig
CortSexAge_robust_int.Sig.t=ggseg(.data=CortSexAge_robust_int.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+labs(title="Sex*Age")+custom_theme
SubcortSexAge_robust_int.Sig.t=ggseg(.data=SubcortSexAge_robust_int.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+custom_theme_sub

# remove unneeded variables 
# rm(CortSex_region,SubcortSex_region,CortSex_robust,SubcortSex_robust,CortSex_region_int,SubcortSex_region_int,CortSex_robust_int,SubcortSex_robust_int,
#    CortCentiloid_region,SubcortCentiloid_region,CortCentiloid_robust,SubcortCentiloid_robust,CortCentiloid_region_int,SubcortCentiloid_region_int,CortCentiloid_robust_int,SubcortCentiloid_robust_int,
#    CortRace_region,SubcortRace_region,CortRace_robust,SubcortRace_robust,CortRace_region_int,SubcortRace_region_int,CortRace_robust_int,SubcortRace_robust_int,
#    CortAge_region,SubcortAge_region,CortAge_robust,SubcortAge_robust,CortAge_region_int,SubcortAge_region_int,CortAge_robust_int,SubcortAge_robust_int,
#    CortSexAge_region_int,SubcortSexAge_region_int,CortSexAge_robust_int,SubcortSexAge_robust_int)


##### Merge brain panels into figures for the paper using cowplot #####
# first we are going to make mini figures for each cortical,subcortical merger
# region
sex_row=plot_grid(CortSex_region.t, SubcortSex_region.t,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
centiloid_row=plot_grid(CortCentiloid_region.t, SubcortCentiloid_region.t,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
race_row=plot_grid(CortRace_region.t, SubcortRace_region.t,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
age_row=plot_grid(CortAge_region.t, SubcortAge_region.t,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") ) 
merged_brain=plot_grid(sex_row, centiloid_row, race_row, age_row, ncol=1)+theme(plot.background = element_rect(fill = "black"), panel.border = element_blank())
# Save Figure
tiff('./Figures/Supplemental_Figure1_Unthresholded_region.tiff', units="in", width=7, height=9, res=200,type="cairo" )
print(merged_brain)
dev.off()
# robust
sex_row=plot_grid(CortSex_robust.t, SubcortSex_robust.t,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
centiloid_row=plot_grid(CortCentiloid_robust.t, SubcortCentiloid_robust.t,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
race_row=plot_grid(CortRace_robust.t, SubcortRace_robust.t,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
age_row=plot_grid(CortAge_robust.t, SubcortAge_robust.t,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") ) 
merged_brain=plot_grid(sex_row, centiloid_row, race_row, age_row, ncol=1)+theme(plot.background = element_rect(fill = "black"), panel.border = element_blank())
# Save Figure
tiff('./Figures/Supplemental_Figure1_Unthresholded_robust.tiff', units="in", width=7, height=9, res=200,type="cairo" )
print(merged_brain)
dev.off()
# region_int
sex_row=plot_grid(CortSex_region_int.t, SubcortSex_region_int.t,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
centiloid_row=plot_grid(CortCentiloid_region_int.t, SubcortCentiloid_region_int.t,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
race_row=plot_grid(CortRace_region_int.t, SubcortRace_region_int.t,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
age_row=plot_grid(CortAge_region_int.t, SubcortAge_region_int.t,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") ) 
sex_age_row=plot_grid(CortSexAge_region_int.t, SubcortSexAge_region_int.t,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
merged_brain=plot_grid(sex_row, centiloid_row, race_row, age_row, sex_age_row, ncol=1)+theme(plot.background = element_rect(fill = "black"), panel.border = element_blank())
# Save Figure
tiff('./Figures/Supplemental_Figure1_Unthresholded_region_int.tiff', units="in", width=7, height=9, res=200,type="cairo" )
print(merged_brain)
dev.off()
# robust_int
sex_row=plot_grid(CortSex_robust_int.t, SubcortSex_robust_int.t,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
centiloid_row=plot_grid(CortCentiloid_robust_int.t, SubcortCentiloid_robust_int.t,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
race_row=plot_grid(CortRace_robust_int.t, SubcortRace_robust_int.t,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
age_row=plot_grid(CortAge_robust_int.t, SubcortAge_robust_int.t,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") ) 
sex_age_row=plot_grid(CortSexAge_robust_int.t, SubcortSexAge_robust_int.t,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
merged_brain=plot_grid(sex_row, centiloid_row, race_row, age_row, sex_age_row, ncol=1)+theme(plot.background = element_rect(fill = "black"), panel.border = element_blank())
# Save Figure
tiff('./Figures/Supplemental_Figure1_Unthresholded_robust_int.tiff', units="in", width=7, height=9, res=200,type="cairo" )
print(merged_brain)
dev.off()

# now we are going to make mini figures for each significant cortical,subcortical merger
# region 
sex_row_sig=plot_grid(CortSex_region.Sig.t, SubcortSex_region.Sig.t,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
age_row_sig=plot_grid(CortAge_region.Sig.t, SubcortAge_region.Sig.t,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") ) 
merged_brain_sig=plot_grid(sex_row_sig, age_row_sig, ncol=1)+theme(plot.background = element_rect(fill = "black"), panel.border = element_blank())
# Save Figure
tiff('./Figures/Figure2_Thresholded_region.tiff', units="in", width=7, height=9, res=200,type="cairo" )
print(merged_brain_sig)
dev.off()
# robust
sex_row_sig=plot_grid(CortSex_robust.Sig.t, SubcortSex_robust.Sig.t,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
age_row_sig=plot_grid(CortAge_robust.Sig.t, SubcortAge_robust.Sig.t,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") ) 
merged_brain_sig=plot_grid(sex_row_sig, age_row_sig, ncol=1)+theme(plot.background = element_rect(fill = "black"), panel.border = element_blank())
# Save Figure
tiff('./Figures/Figure2_Thresholded_robust.tiff', units="in", width=7, height=9, res=200,type="cairo" )
print(merged_brain_sig)
dev.off()
# region_int
sex_row_sig=plot_grid(CortSex_region_int.Sig.t, SubcortSex_region_int.Sig.t,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
age_row_sig=plot_grid(CortAge_region_int.Sig.t, SubcortAge_region_int.Sig.t,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") ) 
sex_age_row_sig=plot_grid(CortSexAge_region_int.Sig.t, SubcortSexAge_region_int.Sig.t,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
merged_brain_sig=plot_grid(sex_row_sig, age_row_sig, sex_age_row_sig, ncol=1)+theme(plot.background = element_rect(fill = "black"), panel.border = element_blank())
# Save Figure
tiff('./Figures/Figure2_Thresholded_region_int.tiff', units="in", width=7, height=9, res=200,type="cairo" )
print(merged_brain_sig)
dev.off()
# robust_int
sex_row_sig=plot_grid(CortSex_robust_int.Sig.t, SubcortSex_robust_int.Sig.t,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
age_row_sig=plot_grid(CortAge_robust_int.Sig.t, SubcortAge_robust_int.Sig.t,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") ) 
sex_age_row_sig=plot_grid(CortSexAge_robust_int.Sig.t, SubcortSexAge_robust_int.Sig.t,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
merged_brain_sig=plot_grid(sex_row_sig, age_row_sig, sex_age_row_sig, ncol=1)+theme(plot.background = element_rect(fill = "black"), panel.border = element_blank())
# Save Figure
tiff('./Figures/Figure2_Thresholded_robust_int.tiff', units="in", width=7, height=9, res=200,type="cairo" )
print(merged_brain_sig)
dev.off()


##### Make Age Scatterplots #####
scatter_theme = list(
  geom_point(size=2),  geom_smooth(method=lm, se=FALSE, fullrange=TRUE, aes(colour = factor(sex))),  scale_color_manual(values=c('#3182BD','#CC6666')), #gold #E69F00
  theme_bw(), theme(text = element_text(size=16)),  theme(legend.position="none"), theme(legend.title=element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()), scale_y_continuous(
    labels = scales::number_format(accuracy = 0.1)))

# Regions with the largest age effects:
# region: cortical(16) - isthmus cingulate, lateral orbitofrontal,inferior temporal, pericalcarine, precuneus
#         subcortical(7) - putamen, pallidum, caudate, hippocampus, ventral DC
# robust: cortical(15) -  isthmus cingulate, lateral orbitofrontal, inferior temporal, bankssts, pericalcarine
#         subcortical(7)- putamen, pallidum, caudate, hippocampus, ventral DC

# Cortical
age_isthcing=ggplot(tauage, aes(x=age, y=isthmus.cingulate, linetype=sex )) +
  xlab("Age") + ylab("Isthmus Cingulate Tau")+  scatter_theme

age_latorb=ggplot(tauage, aes(x=age, y=lateral.orbitofrontal, linetype=sex )) +
  xlab("Age") + ylab("Lateral Orbital Frontal Tau")+  scatter_theme

age_inftemp=ggplot(tauage, aes(x=age, y=inferior.temporal, linetype=sex )) +
  xlab("Age") + ylab("Inferior Temporal Tau")+  scatter_theme

age_pericalc=ggplot(tauage, aes(x=age, y=pericalcarine, linetype=sex )) +
  xlab("Age") + ylab("Pericalcarine Tau")+  scatter_theme

age_precuneus=ggplot(tauage, aes(x=age, y=precuneus, linetype=sex )) +
  xlab("Age") + ylab("Precuneus Tau")+  scatter_theme

age_bank=ggplot(tauage, aes(x=age, y=bankssts, linetype=sex )) +
  xlab("Age") + ylab("Bankssts Tau")+  scatter_theme


# Subcortical: pallidum, putamen, amygdala
age_putamen=ggplot(tauage, aes(x=age, y=putamen, linetype=sex)) +
  xlab("Age") + ylab("Putamen Tau")+scatter_theme+ theme(axis.title.x = element_blank())

age_pallidum=ggplot(tauage, aes(x=age, y=pallidum, linetype=sex )) +
  xlab("Age") + ylab("Pallidum Tau")+  scatter_theme

age_caud=ggplot(tauage, aes(x=age, y=caudate, linetype=sex )) +
  xlab("Age") + ylab("Caudate Tau")+  scatter_theme + theme(axis.title.x = element_blank())

age_hipp=ggplot(tauage, aes(x=age, y=hippocampus, linetype=sex )) +
  xlab("Age") + ylab("Hippocampus Tau")+  scatter_theme + theme(axis.title.x = element_blank())

age_ventdc=ggplot(tauage, aes(x=age, y=ventral.DC, linetype=sex )) +
  xlab("Age") + ylab("Ventral DC Tau")+  scatter_theme + theme(axis.title.x = element_blank())

## Interaction of age & sex 
# region_int: cortical(0) - 0
#             subcortical(0)- 0 
# robust_int: cortical(0) - 0
#             subcortical(0) - 0

##### Make Sex Violin Plots #####
# Greatest sex effects are in rostral middle frontal, frontal pole, and middle temporal gyrus . Do RMF, then inferior temporal, and entorhinal 
violin_theme=list(
  geom_violin(trim=FALSE, alpha=1),  geom_boxplot(width=0.1, fill="white", outlier.shape = NA),scale_fill_manual(values=c('#3182BD','#CC6666'), labels=c("Men","Women")), #gold #E69F00
  theme_classic(), theme(text = element_text(size=16)),  theme(legend.position="none"), theme(legend.title=element_blank()), scale_x_discrete(breaks=c("male","female"), labels=c("Men","Women")), theme(axis.title.x = element_blank()), scale_y_continuous(
    labels = scales::number_format(accuracy = 0.1))) # to add dots instead of box plots geom_dotplot(binaxis='y', stackdir='center', dotsize=1)

# Regions with the largest sex effects:
# region: cortical(20) - rostral middle frontal, frontal pole, lateral occipital, inferior parietal, middle temporal, 
#         subcortical(0) - 0
# robust: cortical(17) - frontal pole, rostral middle frontal, lateral occipital, middle temporal, inferior parietal 
#         subcortical(0)- 0

# Cortical: frontal pole, lateral occipital, pars orbitalis, rostral middle frontal, postcentral
sex_frontpole=ggplot(tauage, aes(x=sex, y=frontal.pole, fill=sex))+  ylab("Frontal Pole Tau")+  violin_theme  #xlab("Group") + ylab("Temporal Pole Tau") #If you don't want a figure title
sex_latocc=ggplot(tauage, aes(x=sex, y=lateral.occipital, fill=sex))+ ylab("Lateral Occipital Tau")+  violin_theme  
sex_infpar=ggplot(tauage, aes(x=sex, y=inferior.parietal, fill=sex))+ ylab("Inferior Parietal Tau")+  violin_theme  
sex_rostmidfron=ggplot(tauage, aes(x=sex, y=rostral.middle.frontal, fill=sex))+ ylab("Rostral Middle Frontal Tau")+  violin_theme  
sex_midtemp=ggplot(tauage, aes(x=sex, y=middle.temporal, fill=sex))+ ylab("Middle Temporal Tau")+  violin_theme  


##### Merge panels into figures for the paper using cowplot ####
scatter_age1=plot_grid(age_latorb,age_inftemp,age_bank,age_isthcing,age_pericalc,age_precuneus, ncol=2, labels = c("A)", "B)", "C)","D)", "E)", "F)"), align = "v", axis = "bt")+theme(plot.background = element_rect(fill = "white"), panel.border = element_blank())
tiff('./Figures/Scatter_age_cort.tiff', units="in", width=8.5, height=8.5, res=200,type="cairo" )
print(scatter_age1)
dev.off()

scatter_age2=plot_grid(age_pallidum,age_putamen,age_caud,age_hipp,age_ventdc, ncol=2, labels = c("A)", "B)", "C)","D)", "E)"), align = "v", axis = "bt")+theme(plot.background = element_rect(fill = "white"), panel.border = element_blank())
tiff('./Figures/Scatter_age_subcort.tiff', units="in", width=8.5, height=8.5, res=200,type="cairo" )
print(scatter_age2)
dev.off()


violin_sex=plot_grid(sex_frontpole,sex_latocc,sex_infpar,sex_rostmidfron,sex_midtemp, ncol=2, labels = c("A)","B)","C)","D)","E)"), align = "v", axis = "bt")+theme(plot.background = element_rect(fill = NA), panel.border = element_blank())
tiff('./Figures/Violin_sex.tiff', units="in", width=8.5, height=8.5, res=200,type="cairo" )
print(violin_sex)
dev.off()

##### Write out significance tables as figures #####
#round tables to uniform length using something like round(sex.sig[, 2:5], digits=7)
# region
sex.region.all=rbind(CortSex_region, SubcortSex_region)
sex.region.all.sig=rbind(CortSex_region.Sig, SubcortSex_region.Sig)
centiloid.region.all=rbind(CortCentiloid_region, SubcortCentiloid_region)
race.region.all=rbind(CortRace_region, SubcortRace_region)
age.region.all=rbind(CortAge_region, SubcortAge_region)
age.region.all.sig=rbind(CortAge_region.Sig, SubcortAge_region.Sig)
# robust
sex.robust.all=rbind(CortSex_robust, SubcortSex_robust)
sex.robust.all.sig=rbind(CortSex_robust.Sig, SubcortSex_robust.Sig)
centiloid.robust.all=rbind(CortCentiloid_robust, SubcortCentiloid_robust)
race.robust.all=rbind(CortRace_robust, SubcortRace_robust)
age.robust.all=rbind(CortAge_robust, SubcortAge_robust)
age.robust.all.sig=rbind(CortAge_robust.Sig, SubcortAge_robust.Sig)
# region_int
sex.region_int.all=rbind(CortSex_region_int, SubcortSex_region_int)
centiloid.region_int.all=rbind(CortCentiloid_region_int, SubcortCentiloid_region_int)
race.region_int.all=rbind(CortRace_region_int, SubcortRace_region_int)
age.region_int.all=rbind(CortAge_region_int, SubcortAge_region_int)
sex_age.region_int.all=rbind(CortSexAge_region_int, SubcortSexAge_region_int)
sex_age.region_int.all.sig=rbind(CortSexAge_region_int.Sig, SubcortSexAge_region_int.Sig)
# robust_int
sex.robust_int.all=rbind(CortSex_robust_int, SubcortSex_robust_int)
centiloid.robust_int.all=rbind(CortCentiloid_robust_int, SubcortCentiloid_robust_int)
race.robust_int.all=rbind(CortRace_robust_int, SubcortRace_robust_int)
age.robust_int.all=rbind(CortAge_robust_int, SubcortAge_robust_int)
sex_age.robust_int.all=rbind(CortSexAge_robust_int, SubcortSexAge_robust_int)
sex_age.robust_int.all.sig=rbind(CortSexAge_robust_int.Sig, SubcortSexAge_robust_int.Sig)

# Sex
sex.round=cbind(sex.region.all,sex.robust.all[,2:5])
sex.round[,2]=round(sex.round[,2],digits=3);sex.round[,3]=round(sex.round[,3],digits=2);sex.round[,4]=round(sex.round[,4],digits=10);sex.round[,5]=round(sex.round[,5],digits=10);
sex.round[,6]=round(sex.round[,6],digits=3);sex.round[,7]=round(sex.round[,7],digits=2);sex.round[,8]=round(as.numeric(sex.round[,8]),digits=10);sex.round[,9]=round(sex.round[,9],digits=10)

sex.round.region.sig=sex.region.all.sig
sex.round.region.sig[,2]=round(sex.round.region.sig[,2],digits=3); sex.round.region.sig[,3]=round(sex.round.region.sig[,3],digits=2); sex.round.region.sig[,4]=round(sex.round.region.sig[,4],digits=10); sex.round.region.sig[,5]=round(sex.round.region.sig[,5],digits=10)

sex.round.robust.sig=sex.robust.all.sig
sex.round.robust.sig[,2]=round(sex.round.robust.sig[,2],digits=3); sex.round.robust.sig[,3]=round(sex.round.robust.sig[,3],digits=2); sex.round.robust.sig[,4]=round(as.numeric(sex.round.robust.sig[,4]),digits=10); sex.round.robust.sig[,5]=round(sex.round.robust.sig[,5],digits=10)

# Centiloid
centiloid.round=cbind(centiloid.region.all,centiloid.robust.all[,2:5]); 
centiloid.round[,2]=round(centiloid.round[,2],digits=4); centiloid.round[,3]=round(centiloid.round[,3],digits=2); centiloid.round[,4]=round(centiloid.round[,4],digits=2); centiloid.round[,5]=round(centiloid.round[,5],digits=2)
centiloid.round[,6]=round(centiloid.round[,6],digits=4); centiloid.round[,7]=round(centiloid.round[,7],digits=2); centiloid.round[,8]=round(as.numeric(centiloid.round[,8]),digits=2); centiloid.round[,9]=round(centiloid.round[,9],digits=2)

# Race
race.round=cbind(race.region.all,race.robust.all[,2:5]); 
race.round[,2]=round(race.round[,2],digits=4); race.round[,3]=round(race.round[,3],digits=2); race.round[,4]=round(race.round[,4],digits=2); race.round[,5]=round(race.round[,5],digits=2)
race.round[,6]=round(race.round[,6],digits=4); race.round[,7]=round(race.round[,7],digits=2); race.round[,8]=round(as.numeric(race.round[,8]),digits=2); race.round[,9]=round(race.round[,9],digits=2)

# Age
age.round=cbind(age.region.all,age.robust.all[,2:5])
age.round[,2]=round(age.round[,2],digits=3);age.round[,3]=round(age.round[,3],digits=2);age.round[,4]=round(age.round[,4],digits=10);age.round[,5]=round(age.round[,5],digits=10);
age.round[,6]=round(age.round[,6],digits=3);age.round[,7]=round(age.round[,7],digits=2);age.round[,8]=round(as.numeric(age.round[,8]),digits=10);age.round[,9]=round(age.round[,9],digits=10)

age.round.region.sig=age.region.all.sig
age.round.region.sig=age.region.all.sig; age.round.region.sig[,2]=round(age.round.region.sig[,2],digits=3); age.round.region.sig[,3]=round(age.round.region.sig[,3],digits=2); age.round.region.sig[,4]=round(age.round.region.sig[,4],digits=10); age.round.region.sig[,5]=round(age.round.region.sig[,5],digits=10)

age.round.robust.sig=age.robust.all.sig
age.round.robust.sig=age.robust.all.sig; age.round.robust.sig[,2]=round(age.round.robust.sig[,2],digits=3); age.round.robust.sig[,3]=round(age.round.robust.sig[,3],digits=2); age.round.robust.sig[,4]=round(as.numeric(age.round.robust.sig[,4]),digits=10); age.round.robust.sig[,5]=round(age.round.robust.sig[,5],digits=10)

# Sex*Age
sexage.round=cbind(sex_age.region_int.all,sex_age.robust_int.all[,2:5]); 
sexage.round[,2]=round(sexage.round[,2],digits=5); sexage.round[,3]=round(sexage.round[,3],digits=2); sexage.round[,4]=round(sexage.round[,4],digits=2); sexage.round[,5]=round(sexage.round[,5],digits=2)
sexage.round[,6]=round(sexage.round[,6],digits=5); sexage.round[,7]=round(sexage.round[,7],digits=2); sexage.round[,8]=round(as.numeric(sexage.round[,8]),digits=3); sexage.round[,9]=round(sexage.round[,9],digits=2) 

sexage.round.region.sig=sex_age.region_int.all.sig
sexage.round.region.sig[,2]=round(sexage.round.region.sig[,2],digits=3); sexage.round.region.sig[,3]=round(sexage.round.region.sig[,3],digits=2); sexage.round.region.sig[,4]=round(sexage.round.region.sig[,4],digits=10); sexage.round.region.sig[,5]=round(sexage.round.region.sig[,5],digits=10)

sexage.round.robust.sig=sex_age.robust_int.all.sig
sexage.round.robust.sig[,2]=round(sexage.round.robust.sig[,2],digits=3); sexage.round.robust.sig[,3]=round(sexage.round.robust.sig[,3],digits=2); sexage.round.robust.sig[,4]=round(as.numeric(sexage.round.robust.sig[,4]),digits=10); sexage.round.robust.sig[,5]=round(sexage.round.robust.sig[,5],digits=10)

#change headings
colnames(sex.round)[1]= "Region"; colnames(sex.round)[5]= "p (cor)";colnames(sex.round)[9]= "p (adj)"
colnames(sex.round.region.sig)[1]= "Region"; colnames(sex.round.region.sig)[5]= "p (adj)"
colnames(sex.round.robust.sig)[1]= "Region"; colnames(sex.round.robust.sig)[5]= "p (adj)"
colnames(centiloid.round)[1]= "Region"; colnames(centiloid.round)[5]= "p (cor)";colnames(centiloid.round)[9]= "p (adj)";
colnames(race.round)[1]= "Region"; colnames(race.round)[5]= "p (cor)";colnames(race.round)[9]= "p (adj)";
colnames(age.round)[1]= "Region"; colnames(age.round)[5]= "p (cor)";colnames(age.round)[9]= "p (adj)"
colnames(age.round.region.sig)[1]= "Region"; colnames(age.round.region.sig)[5]= "p (adj)"
colnames(age.round.robust.sig)[1]= "Region"; colnames(age.round.robust.sig)[5]= "p (adj)"
colnames(sexage.round)[1]= "Region"; colnames(sexage.round)[5]= "p (cor)";colnames(sexage.round)[9]= "p (adj)";
colnames(sexage.round.region.sig)[1]= "Region"; colnames(sexage.round.region.sig)[5]= "p (adj)"
colnames(sexage.round.robust.sig)[1]= "Region"; colnames(sexage.round.robust.sig)[5]= "p (adj)"

pacman::p_load(gridExtra, grid, gtable)

# Save Figure
#Just significant tables for the main paper
g <- tableGrob(age.round.region.sig, rows = NULL)
g <- gtable_add_grob(g, grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)), t = 2, b = nrow(g), l = 1, r = ncol(g))
g <- gtable_add_grob(g, grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)), t = 1, l = 1, r = ncol(g))
tiff('./Figures/age_region_table1.tiff', units="in", width=5.4, height=5.8, res=200,type="cairo" )
grid.draw(g) #grid.table(age.round_sig, rows=NULL)
dev.off()

g <- tableGrob(age.round.robust.sig, rows = NULL)
g <- gtable_add_grob(g, grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)), t = 2, b = nrow(g), l = 1, r = ncol(g))
g <- gtable_add_grob(g, grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)), t = 1, l = 1, r = ncol(g))
tiff('./Figures/age_robust_table1.tiff', units="in", width=6, height=5.8, res=200,type="cairo" )
grid.draw(g) #grid.table(age.round_sig, rows=NULL)
dev.off()

g <- tableGrob(sex.round.region.sig, rows = NULL)
g <- gtable_add_grob(g, grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)), t = 2, b = nrow(g), l = 1, r = ncol(g))
g <- gtable_add_grob(g, grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)), t = 1, l = 1, r = ncol(g))
tiff('./Figures/sex_region_table2.tiff', units="in", width=5.6, height=6.9, res=200,type="cairo" )
grid.draw(g) #grid.table(sex.round_sig, rows=NULL)
dev.off()

# no observations for sex.robust
# g <- tableGrob(sex.round.robust.sig, rows = NULL)
# g <- gtable_add_grob(g, grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)), t = 2, b = nrow(g), l = 1, r = ncol(g))
# g <- gtable_add_grob(g, grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)), t = 1, l = 1, r = ncol(g))
# tiff('./Figures/sex_robust_table2.tiff', units="in", width=5.6, height=6.9, res=200,type="cairo" )
# grid.draw(g) #grid.table(sex.round_sig, rows=NULL)
# dev.off()


##### Supplemental showing all regions #####

#Supplemental tables that show all regions as well as the robust regression values
g <- tableGrob(age.round, rows = NULL)
g <- gtable_add_grob(g, grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)), t = 2, b = nrow(g), l = 1, r = ncol(g))
g <- gtable_add_grob(g, grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)), t = 1, l = 1, r = ncol(g))
g <- gtable_add_grob(g, grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)), t = 1, b = nrow(g), l = 2, r = 5)
tiff('./Figures/Supplemental1_age.tiff', units="in", width=10, height=12, res=200,type="cairo" )
grid.draw(g);dev.off()

g <- tableGrob(sex.round, rows = NULL)
g <- gtable_add_grob(g, grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)), t = 2, b = nrow(g), l = 1, r = ncol(g))
g <- gtable_add_grob(g, grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)), t = 1, l = 1, r = ncol(g))
g <- gtable_add_grob(g, grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)), t = 1, b = nrow(g), l = 2, r = 5)
tiff('./Figures/Supplemental2_sex.tiff', units="in", width=10, height=12, res=200,type="cairo" )
grid.draw(g);dev.off()

g <- tableGrob(centiloid.round, rows = NULL)
g <- gtable_add_grob(g, grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)), t = 2, b = nrow(g), l = 1, r = ncol(g))
g <- gtable_add_grob(g, grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)), t = 1, l = 1, r = ncol(g))
g <- gtable_add_grob(g, grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)), t = 1, b = nrow(g), l = 2, r = 5)
tiff('./Figures/Supplemental3_centiloid.tiff', units="in", width=7.5, height=12, res=200,type="cairo" )
grid.draw(g);dev.off()

g <- tableGrob(race.round, rows = NULL)
g <- gtable_add_grob(g, grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)), t = 2, b = nrow(g), l = 1, r = ncol(g))
g <- gtable_add_grob(g, grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)), t = 1, l = 1, r = ncol(g))
g <- gtable_add_grob(g, grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)), t = 1, b = nrow(g), l = 2, r = 5)
tiff('./Figures/Supplemental4_race.tiff', units="in", width=7.5, height=12, res=200,type="cairo" )
grid.draw(g);dev.off()

g <- tableGrob(sexage.round, rows = NULL)
g <- gtable_add_grob(g, grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)), t = 2, b = nrow(g), l = 1, r = ncol(g))
g <- gtable_add_grob(g, grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)), t = 1, l = 1, r = ncol(g))
g <- gtable_add_grob(g, grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)), t = 1, b = nrow(g), l = 2, r = 5)
tiff('./Figures/Supplemental5_interaction.tiff', units="in", width=7.5, height=12, res=200,type="cairo" )
grid.draw(g);    dev.off()


tauage$cdr <- as.factor(tauage$cdr)
# CDR group and Age on inferior temporal and hippocampus '
scatter_theme_cdr = list(
  geom_point(size=2),  geom_smooth(method=lm, se=FALSE, fullrange=TRUE, aes(colour = factor(cdr))),  scale_color_manual(values=c('#3182BD','#CC6666')), #gold #E69F00
  theme_bw(), theme(text = element_text(size=16))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()), scale_y_continuous(
    labels = scales::number_format(accuracy = 0.1)))#,  theme(legend.position="none"), theme(legend.title=element_blank())

a <- ggplot(tauage, aes(x=age, y=inferior.temporal, shape=cdr, color=cdr)) +
  geom_point(size=2, alpha=0.5) +
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE, aes(colour = factor(cdr))) +
  scale_color_manual(values=c('#3182BD','#CC6666')) + #gold #E69F00
  labs(title="CDR v Age Inferior Temporal",
       x="Age", y = "Inferior Temporal")+
  theme_classic() 
tiff('./Figures/CDR_x_Age_InfTemp.tiff', units="in", width=8.5, height=8.5, res=200,type="cairo" )
dev.off()

b <- ggplot(tauage, aes(x=age, y=hippocampus, shape=cdr, color=cdr)) +
  geom_point(size=2, alpha=0.5) +
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE, aes(colour = factor(cdr))) +
  scale_color_manual(values=c('#3182BD','#CC6666')) + #gold #E69F00
  labs(title="CDR v Age Hippocampus",
       x="Age", y = "Hippocampus")+
  theme_classic() 
tiff('./Figures/CDR_x_Age_Hippo.tiff', units="in", width=8.5, height=8.5, res=200,type="cairo" )
dev.off()

cdr_x_age=plot_grid(a, b,  align="hv", axis="tb", rel_widths = c(1, 1))+theme(plot.background = element_rect(fill = "black") )  
tiff('./Figures/CDR_x_Age.tiff', units="in", width=8.5, height=8.5, res=200,type="cairo" )
dev.off()



##### Cognitive scores ##### 
cog <- tauage%>%dplyr::select(ID,age,cohort,amy_pos,cdr,ANIMALS,tmb)%>% #n = 255
                na.omit()

control_cog <- cog%>%filter(amy_pos=="A-" & cdr==0) # n = 241

# Baseline control cognitive summary score to get m/sd for subsequent analyses
mANIMALS <- mean(control_cog$ANIMALS, na.rm = T); sdANIMALS <- sd(control_cog$ANIMALS, na.rm = T)
mtmb <- mean(control_cog$tmb, na.rm = T); sdtmb <- sd(control_cog$tmb, na.rm = T)

cog$zANIMALS <- (cog$ANIMALS - mANIMALS)/sdANIMALS
cog$ztmb <- ((cog$tmb - mtmb)/sdtmb)*(-1)

cog$composite <- (cog$zANIMALS + cog$ztmb)/2 

## Regression tau age effect on cognition 
composite <- cog%>%dplyr::select(ID,composite)
data <- merge(tauage, composite, by="ID") # n = 255
# Significant region
data$region_avg = (data$isthmus.cingulate+ 
                     data$lateral.orbitofrontal+
                     data$inferior.temporal+ 
                     data$pericalcarine+ 
                     data$precuneus+ 
                     data$bankssts+
                     data$medial.orbitofrontal+
                     data$posterior.cingulate+
                     data$fusiform+ 
                     data$middle.temporal+
                     data$rostral.anterior.cingulate+
                     data$paracentral+
                     data$pars.triangularis+
                     data$supramarginal+ 
                     data$inferior.parietal+
                     data$frontal.pole)/16

data$robust_avg = (data$isthmus.cingulate+ 
                      data$lateral.orbitofrontal+
                      data$inferior.temporal+ 
                      data$pericalcarine+ 
                      data$precuneus+ 
                      data$bankssts+
                      data$medial.orbitofrontal+
                      data$posterior.cingulate+
                      data$fusiform+ 
                      data$middle.temporal+
                      data$rostral.anterior.cingulate+
                      data$paracentral+
                      data$pars.triangularis+
                      data$supramarginal+ 
                      data$transverse.temporal)/15


# regression cog ~ avg sig age tau regions
region <- lm(composite ~ region_avg, data=data)
region_cv <- lm(composite ~ age+sex+race+region_avg, data=data)
summary(region)
summary(region_cv)

robust <- lm(composite ~ robust_avg, data=data)
robust_cv <- lm(composite ~ age+sex+race+robust_avg, data=data)
summary(robust)
summary(robust_cv)

data$cdr <- as.factor(data$cdr)
aaa<- ggplot(data, aes(x=region_avg, y=composite, shape=cdr, color=cdr)) +
  geom_point(size=2, alpha=0.5) +
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE) +
  scale_color_manual(values=c('#3182BD','#CC6666')) + #gold #E69F00
  labs(title="Region: Cognition x Tau Separated by CDR",
       x="Average Tau", y = "Cognitive Composite Z-Score")+
  theme_classic() 
tiff('./Figures/Region_Cognition_Tau_byCDR.tiff', units="in", width=8.5, height=8.5, res=200,type="cairo" )
dev.off()

bbb<- ggplot(data, aes(x=region_avg, y=composite)) +
  geom_point(size=2, alpha=0.5) +
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE) +
  scale_color_manual(values=c('#3182BD','#CC6666')) + #gold #E69F00
  labs(title="Region: Cognition x Tau",
       x="Average Tau", y = "Cognitive Composite Z-Score")+
  theme_classic()
tiff('./Figures/Region_Cognition_Tau.tiff', units="in", width=8.5, height=8.5, res=200,type="cairo" )
dev.off()

ccc<- ggplot(data, aes(x=robust_avg, y=composite, shape=cdr, color=cdr)) +
  geom_point(size=2, alpha=0.5) +
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE) +
  scale_color_manual(values=c('#3182BD','#CC6666')) + #gold #E69F00
  labs(title="Robust: Cognition x Tau Separated by CDR",
       x="Average Tau", y = "Cognitive Composite Z-Score")+
  theme_classic() 
tiff('./Figures/Robust_Cognition_Tau_byCDR.tiff', units="in", width=8.5, height=8.5, res=200,type="cairo" )
dev.off()

ddd<- ggplot(data, aes(x=robust_avg, y=composite)) +
  geom_point(size=2, alpha=0.5) +
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE) +
  scale_color_manual(values=c('#3182BD','#CC6666')) + #gold #E69F00
  labs(title="Robust: Cognition x Tau",
       x="Average Tau", y = "Cognitive Composite Z-Score")+
  theme_classic()
tiff('./Figures/Robust_Cognition_Tau.tiff', units="in", width=8.5, height=8.5, res=200,type="cairo" )
dev.off()

# Plot the residuals 
cog_lm <- lm(composite ~ age + sex + race, data=data)
tau_lm <- lm(region_avg ~ age + sex + race, data=data)

cog_resid <- cog_lm$residuals
tau_resid <- tau_lm$residuals
ggplot(data, aes(x=tau_resid, y=cog_resid)) +
  geom_point(size=2, alpha=0.5) +
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE) +
  scale_color_manual(values=c('#3182BD','#CC6666')) + #gold #E69F00
  labs(title="Residual Comparison",
       x="Tau Residuals", y = "Cognitive Composite Residuals")+
  theme_classic()

cog_resid_std <- rstandard(cog_lm)
tau_resid_std <- rstandard(tau_lm)
ggplot(data, aes(x=tau_resid_std, y=cog_resid_std)) +
  geom_point(size=2, alpha=0.5) +
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE) +
  scale_color_manual(values=c('#3182BD','#CC6666')) + #gold #E69F00
  labs(title="Standard Residual Comparison",
       x="Tau Residuals", y = "Cognitive Composite Residuals")+
  theme_classic()
