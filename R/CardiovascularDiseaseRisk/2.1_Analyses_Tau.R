# Title: ADRC - Tau Regression Loop 
# Author: Diana Hobbs
# Date: April 2022

#Functions to set up environment
rm(list=ls()) #removes all variables from the current environment
options(scipen = 999) #forces R not to use scientific notation

#### Set Up ####
# Load packages 
pacman::p_load(tidyverse,jtools,naniar,lubridate,lme4,sfsmisc,foreign,MASS,scales,ggseg,cowplot,GeneNet)

###### LOAD & SETUP DATA 
tau <- read.csv("./Data/Rabin_Replication_TAU.csv")%>%
  mutate(CVD = str_replace_all(CVD_Colsp, c("<5%"="0", "5-10%"="1", "10-20%"="2",">20%"="3")))%>%
  na.omit()

# Reassign variable types & set reference 
tau$SEX <- as.factor(tau$SEX); tau$SEX <- relevel(tau$SEX, "Male") 
tau$APOE4 <- as.factor(tau$APOE4); tau$APOE4 <- relevel(tau$APOE4, "0")
tau$CVD <- as.factor(tau$CVD); tau$CVD <- relevel(tau$CVD, "0")

# Z transforming continuous variables 
#test<-tau%>%
#  mutate(BMI_z = scale(BMI))

###### REGRESSION LOOP
ROI <- tau[,20:60]

# Loop through regressions and keep output 
LM.COF<-list() # output for regression
RLM.COF<-list() # output for robust regression
model<- matrix(0, 41, 12)
outcount=1
# Create a 'for' loop 
for (i in 1:ncol(ROI)){
  lm <- lm(ROI[,i] ~ AGE + SEX + APOE4 + followup_time + CVD + Centiloid + CVD*Centiloid, data=tau)
  rlm <- rlm(ROI[,i] ~ AGE + SEX + APOE4 + followup_time + CVD + Centiloid + CVD*Centiloid, data=tau)
  lm.sum <- summary(lm)
  rlm.sum <- summary(rlm)
  lm.cof <- lm.sum$coefficients
  rlm.cof <- rlm.sum$coefficients
  region=colnames(ROI)[i]
  model[outcount,1]=region
  colnames(model) <- c("region", "p.adjust.age", "p.adjust.sex", "p.adjust.apoe", "p.adjust.followup",
                       "p.adjust.CVD1", "p.adjust.CVD2", "p.adjust.CVD3", "p.adjust.centiloid",
                       "p.adjust.CVD1*Cent", "p.adjust.CVD2*Cent", "p.adjust.CVD3*Cent")
  rob_p=f.robftest(rlm, var = "AGE"); model[outcount, 2]=rob_p$p.value 
  rob_p=f.robftest(rlm, var = "SEXFemale"); model[outcount, 3]=rob_p$p.value
  rob_p=f.robftest(rlm, var = "APOE41");model[outcount, 4]=rob_p$p.value 
  rob_p=f.robftest(rlm, var = "followup_time");model[outcount, 5]=rob_p$p.value 
  rob_p=f.robftest(rlm, var = "CVD1");model[outcount, 6]=rob_p$p.value 
  rob_p=f.robftest(rlm, var = "CVD2");model[outcount, 7]=rob_p$p.value 
  rob_p=f.robftest(rlm, var = "CVD3");model[outcount, 8]=rob_p$p.value 
  rob_p=f.robftest(rlm, var = "Centiloid");model[outcount, 9]=rob_p$p.value 
  rob_p=f.robftest(rlm, var = "CVD1:Centiloid");model[outcount, 10]=rob_p$p.value 
  rob_p=f.robftest(rlm, var = "CVD2:Centiloid");model[outcount, 11]=rob_p$p.value 
  rob_p=f.robftest(rlm, var = "CVD3:Centiloid");model[outcount, 12]=rob_p$p.value 

  outcount=outcount+1
  LM.COF[[i]] <- lm.cof
  RLM.COF[[i]] <- rlm.cof
  names(LM.COF)[i] <- colnames(ROI)[i]
  names(RLM.COF)[i] <- colnames(ROI)[i]
}
write.csv(LM.COF,"./Rabin_Replication_Output/Tau/Colsp_LM.COF_1.csv")
write.csv(RLM.COF,"./Rabin_Replication_Output/Tau/Colsp_RLM.COF_1.csv")
rm(lm,lm.cof,LM.COF,lm.sum,rlm,rlm.cof,RLM.COF,rlm.sum,i,ROI,rob_p,outcount,tau,region)

###### Organize loop output in usable format 
Data <- read.csv("./Rabin_Replication_Output/Tau/Colsp_LM.COF_1.csv")

# Extract Estimates 
Estimate <- Data[,grepl(".Estimate", names(Data))]
names(Estimate) = gsub(pattern = ".Estimate", replacement = "", x = names(Estimate))
Estimate <- t(Estimate)
colnames(Estimate) = c("(Intercept)","AGE","SEX","APOE4","followup_time", "CVD1","CVD2","CVD3","Centiloid","CVD1_Cent","CVD2_Cent","CVD3_Cent")
Estimate <- data.frame(Estimate)

TVal <- Data[,grepl(".t.value", names(Data))]
names(TVal) = gsub(pattern = ".t.value", replacement = "", x = names(TVal))
TVal <- t(TVal)
colnames(TVal) <- c("(Intercept)","AGE","SEX","APOE4","followup_time","CVD1","CVD2","CVD3","Centiloid","CVD1_Cent","CVD2_Cent","CVD3_Cent")
TVal <- data.frame(TVal)

PVal <- Data[,grepl(".Pr...t..", names(Data))]
names(PVal) = gsub(pattern = ".Pr...t..", replacement = "", x = names(PVal))
PVal <- t(PVal)
colnames(PVal) <- c("(Intercept)","AGE","SEX","APOE4","followup_time","CVD1","CVD2","CVD3","Centiloid","CVD1_Cent","CVD2_Cent","CVD3_Cent")
PVal <- data.frame(PVal)

# Merge coefficient summaries & Perform FDR correction using Benjamini and Hochberg method
Age <- data.frame("B"=Estimate$AGE,"t"=TVal$AGE,"p"=PVal$AGE,"p_corr"=(p.adjust(PVal$AGE, method = "BH",n = length(PVal$AGE))))
Sex <- data.frame("B"=Estimate$SEX,"t"=TVal$SEX,"p"=PVal$SEX,"p_corr"=(p.adjust(PVal$SEX, method = "BH",n = length(PVal$SEX))))
APOE4 <- data.frame("B"=Estimate$APOE4,"t"=TVal$APOE4,"p"=PVal$APOE4,"p_corr"=(p.adjust(PVal$APOE4, method = "BH",n = length(PVal$APOE4))))
followup_time <- data.frame("B"=Estimate$followup_time,"t"=TVal$followup_time,"p"=PVal$followup_time,"p_corr"=(p.adjust(PVal$followup_time, method = "BH",n = length(PVal$followup_time))))
CVD1 <- data.frame("B"=Estimate$CVD1,"t"=TVal$CVD1,"p"=PVal$CVD1,"p_corr"=(p.adjust(PVal$CVD1, method = "BH",n = length(PVal$CVD1))))
CVD2 <- data.frame("B"=Estimate$CVD2,"t"=TVal$CVD2,"p"=PVal$CVD2,"p_corr"=(p.adjust(PVal$CVD2, method = "BH",n = length(PVal$CVD2))))
CVD3 <- data.frame("B"=Estimate$CVD3,"t"=TVal$CVD3,"p"=PVal$CVD3,"p_corr"=(p.adjust(PVal$CVD3, method = "BH",n = length(PVal$CVD3))))
Centiloid <- data.frame("B"=Estimate$Centiloid,"t"=TVal$Centiloid,"p"=PVal$Centiloid,"p_corr"=(p.adjust(PVal$Centiloid, method = "BH",n = length(PVal$Centiloid))))
CVD1_Cent <- data.frame("B"=Estimate$CVD1_Cent,"t"=TVal$CVD1_Cent,"p"=PVal$CVD1_Cent,"p_corr"=(p.adjust(PVal$CVD1_Cent, method = "BH",n = length(PVal$CVD1_Cent))))
CVD2_Cent <- data.frame("B"=Estimate$CVD2_Cent,"t"=TVal$CVD2_Cent,"p"=PVal$CVD2_Cent,"p_corr"=(p.adjust(PVal$CVD2_Cent, method = "BH",n = length(PVal$CVD2_Cent))))
CVD3_Cent <- data.frame("B"=Estimate$CVD3_Cent,"t"=TVal$CVD3_Cent,"p"=PVal$CVD3_Cent,"p_corr"=(p.adjust(PVal$CVD3_Cent, method = "BH",n = length(PVal$CVD3_Cent))))
rm(Data,Estimate,TVal,PVal)

cort_list <- read.csv("./BrianFiles/ggseg_namelist.csv")
subcort_list <- cort_list$subcortical[1:7] #extract subcortical list

CortAge <- data.frame(region=cort_list$cortical,filter(Age[1:34,])); SubcortAge <- data.frame(region=subcort_list,filter(Age[35:41,]))
CortSex <- data.frame(region=cort_list$cortical,filter(Sex[1:34,])); SubcortSex <- data.frame(region=subcort_list,filter(Sex[35:41,]))
CortAPOE4 <- data.frame(region=cort_list$cortical,filter(APOE4[1:34,])); SubcortAPOE4 <- data.frame(region=subcort_list,filter(APOE4[35:41,]))
Cortfollowup_time <- data.frame(region=cort_list$cortical,filter(followup_time[1:34,])); Subcortfollowup_time <- data.frame(region=subcort_list,filter(followup_time[35:41,]))
CortCVD1 <- data.frame(region=cort_list$cortical,filter(CVD1[1:34,])); SubcortCVD1 <- data.frame(region=subcort_list,filter(CVD1[35:41,]))
CortCVD2 <- data.frame(region=cort_list$cortical,filter(CVD2[1:34,])); SubcortCVD2 <- data.frame(region=subcort_list,filter(CVD2[35:41,]))
CortCVD3 <- data.frame(region=cort_list$cortical,filter(CVD3[1:34,])); SubcortCVD3 <- data.frame(region=subcort_list,filter(CVD3[35:41,]))
CortCentiloid <- data.frame(region=cort_list$cortical,filter(Centiloid[1:34,])); SubcortCentiloid <- data.frame(region=subcort_list,filter(Centiloid[35:41,]))
CortCVD1_Cent <- data.frame(region=cort_list$cortical,filter(CVD1_Cent[1:34,])); SubcortCVD1_Cent <- data.frame(region=subcort_list,filter(CVD1_Cent[35:41,]))
CortCVD2_Cent <- data.frame(region=cort_list$cortical,filter(CVD2_Cent[1:34,])); SubcortCVD2_Cent <- data.frame(region=subcort_list,filter(CVD2_Cent[35:41,]))
CortCVD3_Cent <- data.frame(region=cort_list$cortical,filter(CVD3_Cent[1:34,])); SubcortCVD3_Cent <- data.frame(region=subcort_list,filter(CVD3_Cent[35:41,]))
rm(Age,Sex,APOE4,followup_time,CVD1,CVD2,CVD3,Centiloid,CVD1_Cent,CVD2_Cent,CVD3_Cent)

CortAge.Sig <- subset(CortAge, p_corr < 0.05); SubcortAge.Sig <- subset(SubcortAge, p_corr < 0.05)
CortSex.Sig <- subset(CortSex, p_corr < 0.05); SubcortSex.Sig <- subset(SubcortSex, p_corr < 0.05)
CortAPOE4.Sig <- subset(CortAPOE4, p_corr < 0.05); SubcortAPOE4.Sig <- subset(SubcortAPOE4, p_corr < 0.05)
Cortfollowup_time.Sig <- subset(Cortfollowup_time, p_corr < 0.05); Subcortfollowup_time.Sig <- subset(Subcortfollowup_time, p_corr < 0.05)
CortCVD1.Sig <- subset(CortCVD1, p_corr < 0.05); SubcortCVD1.Sig <- subset(SubcortCVD1, p_corr < 0.05)
CortCVD2.Sig <- subset(CortCVD2, p_corr < 0.05); SubcortCVD2.Sig <- subset(SubcortCVD2, p_corr < 0.05)
CortCVD3.Sig <- subset(CortCVD3, p_corr < 0.05); SubcortCVD3.Sig <- subset(SubcortCVD3, p_corr < 0.05)
CortCentiloid.Sig <- subset(CortCentiloid, p_corr < 0.05); SubcortCentiloid.Sig <- subset(SubcortCentiloid, p_corr < 0.05)
CortCVD1_Cent.Sig <- subset(CortCVD1_Cent, p_corr < 0.05); SubcortCVD1_Cent.Sig <- subset(SubcortCVD1_Cent, p_corr < 0.05)
CortCVD2_Cent.Sig <- subset(CortCVD2_Cent, p_corr < 0.05); SubcortCVD2_Cent.Sig <- subset(SubcortCVD2_Cent, p_corr < 0.05)
CortCVD3_Cent.Sig <- subset(CortCVD3_Cent, p_corr < 0.05); SubcortCVD3_Cent.Sig <- subset(SubcortCVD3_Cent, p_corr < 0.05)


###### Make images 
# This theme specifies parameters for all of the plots. Change here to change any of them
custom_theme = list(
  scale_fill_gradientn(colours=c("#046265","#59C3B0","#E6E2BF","#EDA153","#D94227"), limits=c(-4,4), oob=squish),
  #use scale_fill_gradient2 and the option midpoint to specify white as another option
  theme(text = element_text(colour = "white")),  theme(axis.text = element_blank(), axis.title = element_blank()), 
  theme(plot.background = element_rect(fill = "black",colour = NA)), theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), theme(panel.border = element_blank()), guides(fill="none", color="none"))

custom_theme_sub = list(
  scale_fill_gradientn(colours=c("#046265","#59C3B0","#E6E2BF","#EDA153","#D94227"), limits=c(-4,4), oob=squish),
  theme(text = element_text(colour = "white")),  theme(axis.text = element_blank(), axis.title = element_blank()), 
  theme(plot.background = element_rect(fill = "black",colour = NA)), theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), theme(panel.border = element_blank()))

# Plotting t-values on the brain atlas
CortAge.t=ggseg(.data=CortAge,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+
  labs(title="Age") +  custom_theme
print(CortAge.t)
SubcortAge.t=ggseg(.data=SubcortAge,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+
  #labs(title="Age") + 
  custom_theme_sub
print(SubcortAge.t)

CortSex.t=ggseg(.data=CortSex,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+
  labs(title="Sex") +  custom_theme
print(CortSex.t)
SubcortSex.t=ggseg(.data=SubcortSex,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+
  #labs(title="Sex") + 
  custom_theme_sub
print(SubcortSex.t)

CortAPOE4.t=ggseg(.data=CortAPOE4,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+
  labs(title="APOE4") +  custom_theme
print(CortAPOE4.t)
SubcortAPOE4.t=ggseg(.data=SubcortAPOE4,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+
  #labs(title="APOE4") + 
  custom_theme_sub
print(SubcortAPOE4.t)

Cortfollowup_time.t=ggseg(.data=Cortfollowup_time,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+
  labs(title="followup_time") +  custom_theme
print(Cortfollowup_time.t)
Subcortfollowup_time.t=ggseg(.data=Subcortfollowup_time,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+
  #labs(title="APOE4") + 
  custom_theme_sub
print(Subcortfollowup_time.t)

CortCVD1.t=ggseg(.data=CortCVD1,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+
  labs(title="CVD1") +  custom_theme
print(CortCVD1.t)
SubcortCVD1.t=ggseg(.data=SubcortCVD1,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+
  #labs(title="CVD1") + 
  custom_theme_sub
print(SubcortCVD1.t)

CortCVD2.t=ggseg(.data=CortCVD2,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+
  labs(title="CVD2") +  custom_theme
print(CortCVD2.t)
SubcortCVD2.t=ggseg(.data=SubcortCVD2,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+
  #labs(title="CVD2") + 
  custom_theme_sub
print(SubcortCVD2.t)

CortCVD3.t=ggseg(.data=CortCVD3,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+
  labs(title="CVD3") +  custom_theme
print(CortCVD3.t)
SubcortCVD3.t=ggseg(.data=SubcortCVD3,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+
  #labs(title="CVD3") + 
  custom_theme_sub
print(SubcortCVD3.t)

CortCentiloid.t=ggseg(.data=CortCentiloid,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+
  labs(title="Centiloid") +  custom_theme
print(CortCentiloid.t)
SubcortCentiloid.t=ggseg(.data=SubcortCentiloid,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+
  #labs(title="Centiloid") + 
  custom_theme_sub
print(SubcortCentiloid.t)

CortCVD1_Cent.t=ggseg(.data=CortCVD1_Cent,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+
  labs(title="CVD1_Cent") +  custom_theme
print(CortCVD1_Cent.t)
SubcortCVD1_Cent.t=ggseg(.data=SubcortCVD1_Cent,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+
  #labs(title="CVD1") + 
  custom_theme_sub
print(SubcortCVD1_Cent.t)

CortCVD2_Cent.t=ggseg(.data=CortCVD2_Cent,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+
  labs(title="CVD2_Cent") +  custom_theme
print(CortCVD2_Cent.t)
SubcortCVD2_Cent.t=ggseg(.data=SubcortCVD2_Cent,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+
  #labs(title="CVD2") + 
  custom_theme_sub
print(SubcortCVD2_Cent.t)

CortCVD3_Cent.t=ggseg(.data=CortCVD3_Cent,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+
  labs(title="CVD3_Cent") +  custom_theme
print(CortCVD3_Cent.t)
SubcortCVD3_Cent.t=ggseg(.data=SubcortCVD3_Cent,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+
  #labs(title="CVD3") + 
  custom_theme_sub
print(SubcortCVD3_Cent.t)


# only significant ones
CortAge.Sig.t=ggseg(.data=CortAge.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+
  labs(title="Age") +  custom_theme
print(CortAge.Sig.t)
SubcortAge.Sig.t=ggseg(.data=SubcortAge.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+
  #labs(title="Age") + 
  custom_theme_sub
print(SubcortAge.Sig.t)

CortSex.Sig.t=ggseg(.data=CortSex.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+
  labs(title="Sex") +  custom_theme
print(CortSex.Sig.t)
SubcortSex.Sig.t=ggseg(.data=SubcortSex.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+
  #labs(title="Sex") + 
  custom_theme_sub
print(SubcortSex.Sig.t)

CortAPOE4.Sig.t=ggseg(.data=CortAPOE4.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+
  labs(title="APOE4") +  custom_theme
print(CortAPOE4.Sig.t)
SubcortAPOE4.Sig.t=ggseg(.data=SubcortAPOE4.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+
  #labs(title="APOE4") + 
  custom_theme_sub
print(SubcortAPOE4.Sig.t)

Cortfollowup_time.Sig.t=ggseg(.data=Cortfollowup_time.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+
  labs(title="followup_time") +  custom_theme
print(Cortfollowup_time.Sig.t)
Subcortfollowup_time.Sig.t=ggseg(.data=Subcortfollowup_time.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+
  #labs(title="APOE4") + 
  custom_theme_sub
print(Subcortfollowup_time.Sig.t)

CortCVD1.Sig.t=ggseg(.data=CortCVD1.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+
  labs(title="CVD1") +  custom_theme
print(CortCVD1.Sig.t)
SubcortCVD1.Sig.t=ggseg(.data=SubcortCVD1.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+
  #labs(title="CVD1") + 
  custom_theme_sub
print(SubcortCVD1.Sig.t)

CortCVD2.Sig.t=ggseg(.data=CortCVD2.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+
  labs(title="CVD2") +  custom_theme
print(CortCVD2.Sig.t)
SubcortCVD2.Sig.t=ggseg(.data=SubcortCVD2.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+
  #labs(title="CVD2") + 
  custom_theme_sub
print(SubcortCVD2.Sig.t)

CortCVD3.Sig.t=ggseg(.data=CortCVD3.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+
  labs(title="CVD3") +  custom_theme
print(CortCVD3.Sig.t)
SubcortCVD3.Sig.t=ggseg(.data=SubcortCVD3.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+
  #labs(title="CVD3") + 
  custom_theme_sub
print(SubcortCVD3.Sig.t)

CortCentiloid.Sig.t=ggseg(.data=CortCentiloid.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+
  labs(title="Centiloid") +  custom_theme
print(CortCentiloid.Sig.t)
SubcortCentiloid.Sig.t=ggseg(.data=SubcortCentiloid.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+
  #labs(title="Centiloid") + 
  custom_theme_sub
print(SubcortCentiloid.Sig.t)

CortCVD1_Cent.Sig.t=ggseg(.data=CortCVD1_Cent.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+
  labs(title="CVD1_Cent") +  custom_theme
print(CortCVD1_Cent.Sig.t)
SubcortCVD1_Cent.Sig.t=ggseg(.data=SubcortCVD1_Cent.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+
  #labs(title="CVD1") + 
  custom_theme_sub
print(SubcortCVD1_Cent.Sig.t)

CortCVD2_Cent.Sig.t=ggseg(.data=CortCVD2_Cent.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+
  labs(title="CVD2_Cent") +  custom_theme
print(CortCVD2_Cent.Sig.t)
SubcortCVD2_Cent.Sig.t=ggseg(.data=SubcortCVD2_Cent.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+
  #labs(title="CVD2") + 
  custom_theme_sub
print(SubcortCVD2_Cent.Sig.t)

CortCVD3_Cent.Sig.t=ggseg(.data=CortCVD3_Cent.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+
  labs(title="CVD3_Cent") +  custom_theme
print(CortCVD3_Cent.Sig.t)
SubcortCVD3_Cent.Sig.t=ggseg(.data=SubcortCVD3_Cent.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+
  #labs(title="CVD3") + 
  custom_theme_sub
print(SubcortCVD3_Cent.Sig.t)

rm(CortAge,SubcortAge,CortAge.Sig,SubcortAge.Sig,CortSex,SubcortSex,CortSex.Sig,SubcortSex.Sig,CortAPOE4,SubcortAPOE4,CortAPOE4.Sig,SubcortAPOE4.Sig,
   Cortfollowup_time,Subcortfollowup_time,Cortfollowup_time.Sig,Subcortfollowup_time.Sig,
   CortCVD1,SubcortCVD1,CortCVD1.Sig,SubcortCVD1.Sig,CortCVD2,SubcortCVD2,CortCVD2.Sig,SubcortCVD2.Sig,CortCVD3,SubcortCVD3,CortCVD3.Sig,SubcortCVD3.Sig,
   CortCentiloid,SubcortCentiloid,CortCentiloid.Sig,SubcortCentiloid.Sig,CortCVD1_Cent,SubcortCVD1_Cent,CortCVD1_Cent.Sig,SubcortCVD1_Cent.Sig,
   CortCVD2_Cent,SubcortCVD2_Cent,CortCVD2_Cent.Sig,SubcortCVD2_Cent.Sig,CortCVD3_Cent,SubcortCVD3_Cent,CortCVD3_Cent.Sig,SubcortCVD3_Cent.Sig)

# combine plots to get figures
age_row=plot_grid(CortAge.t, SubcortAge.t,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") ) 
sex_row=plot_grid(CortSex.t, SubcortSex.t,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
apoe_row=plot_grid(CortAPOE4.t, SubcortAPOE4.t,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
followup_time_row=plot_grid(Cortfollowup_time.t, Subcortfollowup_time.t,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
cvd1_row=plot_grid(CortCVD1.t, SubcortCVD1.t,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
cvd2_row=plot_grid(CortCVD2.t, SubcortCVD2.t,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
cvd3_row=plot_grid(CortCVD3.t, SubcortCVD3.t,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
Centiloid_row=plot_grid(CortCentiloid.t, SubcortCentiloid.t,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
cvd1_cent_row=plot_grid(CortCVD1_Cent.t, SubcortCVD1_Cent.t,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
cvd2_cent_row=plot_grid(CortCVD2_Cent.t, SubcortCVD2_Cent.t,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
cvd3_cent_row=plot_grid(CortCVD3_Cent.t, SubcortCVD3_Cent.t,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
merged_brain=plot_grid(age_row, sex_row, apoe_row, followup_time_row, cvd1_row, cvd2_row, cvd3_row, Centiloid_row, cvd1_cent_row, cvd2_cent_row, cvd3_cent_row, ncol=1)+theme(plot.background = element_rect(fill = "black"), panel.border = element_blank())
# Save Figure
tiff('./Rabin_Replication_Output/Tau/Colsp_Unthresholded_1.tiff', units="in", width=7, height=9, res=200,type="cairo" )
print(merged_brain)
dev.off()

age_row=plot_grid(CortAge.Sig.t, SubcortAge.Sig.t,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") ) 
sex_row=plot_grid(CortSex.Sig.t, SubcortSex.Sig.t,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
apoe_row=plot_grid(CortAPOE4.Sig.t, SubcortAPOE4.Sig.t,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
followup_time_row=plot_grid(Cortfollowup_time.Sig.t, Subcortfollowup_time.Sig.t,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
cvd1_row=plot_grid(CortCVD1.Sig.t, SubcortCVD1.Sig.t,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
cvd2_row=plot_grid(CortCVD2.Sig.t, SubcortCVD2.Sig.t,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
cvd3_row=plot_grid(CortCVD3.Sig.t, SubcortCVD3.Sig.t,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
Centiloid_row=plot_grid(CortCentiloid.Sig.t, SubcortCentiloid.Sig.t,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
cvd1_cent_row=plot_grid(CortCVD1_Cent.Sig.t, SubcortCVD1_Cent.Sig.t,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
cvd2_cent_row=plot_grid(CortCVD2_Cent.Sig.t, SubcortCVD2_Cent.Sig.t,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
cvd3_cent_row=plot_grid(CortCVD3_Cent.Sig.t, SubcortCVD3_Cent.Sig.t,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
merged_brain.sig=plot_grid(age_row, sex_row, apoe_row, followup_time_row, cvd1_row, cvd2_row, cvd3_row, Centiloid_row, cvd1_cent_row, cvd2_cent_row, cvd3_cent_row, ncol=1)+theme(plot.background = element_rect(fill = "black"), panel.border = element_blank())
# Save Figure
tiff('./Rabin_Replication_Output/Tau/Colsp_Thresholded_1.tiff', units="in", width=7, height=9, res=200,type="cairo")
print(merged_brain.sig)
dev.off()
rm(age_row, sex_row, Centiloid_row, apoe_row, followup_time_row, cvd1_row, cvd2_row, cvd3_row, cvd1_cent_row, cvd2_cent_row, cvd3_cent_row, merged_brain,merged_brain.sig)

###### Robust regression ######
Data <- read.csv("./Rabin_Replication_Output/Tau/Colsp_RLM.COF_1.csv")
# Extract Estimates 
Estimate <- Data[,grepl(".Value", names(Data))]
names(Estimate) = gsub(pattern = ".Estimate", replacement = "", x = names(Estimate))
Estimate <- t(Estimate)
colnames(Estimate) = c("(Intercept)","AGE","SEX","APOE4","followup_time","CVD1","CVD2","CVD3","Centiloid","CVD1_Cent","CVD2_Cent","CVD3_Cent")
Estimate <- data.frame(Estimate)

TVal <- Data[,grepl(".t.value", names(Data))]
names(TVal) = gsub(pattern = ".t.value", replacement = "", x = names(TVal))
TVal <- t(TVal)
colnames(TVal) <- c("(Intercept)","AGE","SEX","APOE4","followup_time","CVD1","CVD2","CVD3","Centiloid","CVD1_Cent","CVD2_Cent","CVD3_Cent")
TVal <- data.frame(TVal)

PVal <- data.frame(model)
names(PVal) = gsub(pattern = "P.adjust.", replacement = "", x = names(PVal))
colnames(PVal) <- c("region","AGE","SEX","APOE4","followup_time","CVD1","CVD2","CVD3","Centiloid","CVD1_Cent","CVD2_Cent","CVD3_Cent")

# Merge coefficient summaries 
Age <- data.frame("B"=Estimate$AGE,"t"=TVal$AGE,"p"=PVal$AGE,"p_corr"=(p.adjust(PVal$AGE, method = "BH",n = length(PVal$AGE))))
Sex <- data.frame("B"=Estimate$SEX,"t"=TVal$SEX,"p"=PVal$SEX,"p_corr"=(p.adjust(PVal$SEX, method = "BH",n = length(PVal$SEX))))
APOE4 <- data.frame("B"=Estimate$APOE4,"t"=TVal$APOE4,"p"=PVal$APOE4,"p_corr"=(p.adjust(PVal$APOE4, method = "BH",n = length(PVal$APOE4))))
followup_time <- data.frame("B"=Estimate$followup_time,"t"=TVal$followup_time,"p"=PVal$followup_time,"p_corr"=(p.adjust(PVal$followup_time, method = "BH",n = length(PVal$followup_time))))
CVD1 <- data.frame("B"=Estimate$CVD1,"t"=TVal$CVD1,"p"=PVal$CVD1,"p_corr"=(p.adjust(PVal$CVD1, method = "BH",n = length(PVal$CVD1))))
CVD2 <- data.frame("B"=Estimate$CVD2,"t"=TVal$CVD2,"p"=PVal$CVD2,"p_corr"=(p.adjust(PVal$CVD2, method = "BH",n = length(PVal$CVD2))))
CVD3 <- data.frame("B"=Estimate$CVD3,"t"=TVal$CVD3,"p"=PVal$CVD3,"p_corr"=(p.adjust(PVal$CVD3, method = "BH",n = length(PVal$CVD3))))
Centiloid <- data.frame("B"=Estimate$Centiloid,"t"=TVal$Centiloid,"p"=PVal$Centiloid,"p_corr"=(p.adjust(PVal$Centiloid, method = "BH",n = length(PVal$Centiloid))))
CVD1_Cent <- data.frame("B"=Estimate$CVD1_Cent,"t"=TVal$CVD1_Cent,"p"=PVal$CVD1_Cent,"p_corr"=(p.adjust(PVal$CVD1_Cent, method = "BH",n = length(PVal$CVD1_Cent))))
CVD2_Cent <- data.frame("B"=Estimate$CVD2_Cent,"t"=TVal$CVD2_Cent,"p"=PVal$CVD2_Cent,"p_corr"=(p.adjust(PVal$CVD2_Cent, method = "BH",n = length(PVal$CVD2_Cent))))
CVD3_Cent <- data.frame("B"=Estimate$CVD3_Cent,"t"=TVal$CVD3_Cent,"p"=PVal$CVD3_Cent,"p_corr"=(p.adjust(PVal$CVD3_Cent, method = "BH",n = length(PVal$CVD3_Cent))))
rm(Data,Estimate,TVal,PVal)

CortAge <- data.frame(region=cort_list$cortical,filter(Age[1:34,])); SubcortAge <- data.frame(region=subcort_list,filter(Age[35:41,]))
CortSex <- data.frame(region=cort_list$cortical,filter(Sex[1:34,])); SubcortSex <- data.frame(region=subcort_list,filter(Sex[35:41,]))
CortAPOE4 <- data.frame(region=cort_list$cortical,filter(APOE4[1:34,])); SubcortAPOE4 <- data.frame(region=subcort_list,filter(APOE4[35:41,]))
Cortfollowup_time <- data.frame(region=cort_list$cortical,filter(followup_time[1:34,])); Subcortfollowup_time <- data.frame(region=subcort_list,filter(followup_time[35:41,]))
CortCVD1 <- data.frame(region=cort_list$cortical,filter(CVD1[1:34,])); SubcortCVD1 <- data.frame(region=subcort_list,filter(CVD1[35:41,]))
CortCVD2 <- data.frame(region=cort_list$cortical,filter(CVD2[1:34,])); SubcortCVD2 <- data.frame(region=subcort_list,filter(CVD2[35:41,]))
CortCVD3 <- data.frame(region=cort_list$cortical,filter(CVD3[1:34,])); SubcortCVD3 <- data.frame(region=subcort_list,filter(CVD3[35:41,]))
CortCentiloid <- data.frame(region=cort_list$cortical,filter(Centiloid[1:34,])); SubcortCentiloid <- data.frame(region=subcort_list,filter(Centiloid[35:41,]))
CortCVD1_Cent <- data.frame(region=cort_list$cortical,filter(CVD1_Cent[1:34,])); SubcortCVD1_Cent <- data.frame(region=subcort_list,filter(CVD1_Cent[35:41,]))
CortCVD2_Cent <- data.frame(region=cort_list$cortical,filter(CVD2_Cent[1:34,])); SubcortCVD2_Cent <- data.frame(region=subcort_list,filter(CVD2_Cent[35:41,]))
CortCVD3_Cent <- data.frame(region=cort_list$cortical,filter(CVD3_Cent[1:34,])); SubcortCVD3_Cent <- data.frame(region=subcort_list,filter(CVD3_Cent[35:41,]))
rm(Age,Sex,APOE4,CVD1,CVD2,CVD3,Centiloid,CVD1_Cent,CVD2_Cent,CVD3_Cent)

CortAge.Sig <- subset(CortAge, p_corr < 0.05); SubcortAge.Sig <- subset(SubcortAge, p_corr < 0.05)
CortSex.Sig <- subset(CortSex, p_corr < 0.05); SubcortSex.Sig <- subset(SubcortSex, p_corr < 0.05)
CortAPOE4.Sig <- subset(CortAPOE4, p_corr < 0.05); SubcortAPOE4.Sig <- subset(SubcortAPOE4, p_corr < 0.05)
Cortfollowup_time.Sig <- subset(Cortfollowup_time, p_corr < 0.05); Subcortfollowup_time.Sig <- subset(Subcortfollowup_time, p_corr < 0.05)
CortCVD1.Sig <- subset(CortCVD1, p_corr < 0.05); SubcortCVD1.Sig <- subset(SubcortCVD1, p_corr < 0.05)
CortCVD2.Sig <- subset(CortCVD2, p_corr < 0.05); SubcortCVD2.Sig <- subset(SubcortCVD2, p_corr < 0.05)
CortCVD3.Sig <- subset(CortCVD3, p_corr < 0.05); SubcortCVD3.Sig <- subset(SubcortCVD3, p_corr < 0.05)
CortCentiloid.Sig <- subset(CortCentiloid, p_corr < 0.05); SubcortCentiloid.Sig <- subset(SubcortCentiloid, p_corr < 0.05)
CortCVD1_Cent.Sig <- subset(CortCVD1_Cent, p_corr < 0.05); SubcortCVD1_Cent.Sig <- subset(SubcortCVD1_Cent, p_corr < 0.05)
CortCVD2_Cent.Sig <- subset(CortCVD2_Cent, p_corr < 0.05); SubcortCVD2_Cent.Sig <- subset(SubcortCVD2_Cent, p_corr < 0.05)
CortCVD3_Cent.Sig <- subset(CortCVD3_Cent, p_corr < 0.05); SubcortCVD3_Cent.Sig <- subset(SubcortCVD3_Cent, p_corr < 0.05)

###### Make images 
# Plotting t-values on the brain atlas
CortAge.t_rob=ggseg(.data=CortAge,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+
  labs(title="Age") +  custom_theme
print(CortAge.t_rob)
SubcortAge.t_rob=ggseg(.data=SubcortAge,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+
  #labs(title="Age") + 
  custom_theme_sub
print(SubcortAge.t_rob)

CortSex.t_rob=ggseg(.data=CortSex,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+
  labs(title="Sex") +  custom_theme
print(CortSex.t_rob)
SubcortSex.t_rob=ggseg(.data=SubcortSex,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+
  #labs(title="Sex") + 
  custom_theme_sub
print(SubcortSex.t_rob)

CortAPOE4.t_rob=ggseg(.data=CortAPOE4,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+
  labs(title="APOE4") +  custom_theme
print(CortAPOE4.t_rob)
SubcortAPOE4.t_rob=ggseg(.data=SubcortAPOE4,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+
  #labs(title="APOE4") + 
  custom_theme_sub
print(SubcortAPOE4.t_rob)

Cortfollowup_time.t_rob=ggseg(.data=Cortfollowup_time,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+
  labs(title="followup_time") +  custom_theme
print(Cortfollowup_time.t_rob)
Subcortfollowup_time.t_rob=ggseg(.data=Subcortfollowup_time,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+
  #labs(title="APOE4") + 
  custom_theme_sub
print(Subcortfollowup_time.t_rob)

CortCVD1.t_rob=ggseg(.data=CortCVD1,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+
  labs(title="CVD1") +  custom_theme
print(CortCVD1.t_rob)
SubcortCVD1.t_rob=ggseg(.data=SubcortCVD1,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+
  #labs(title="CVD1") + 
  custom_theme_sub
print(SubcortCVD1.t_rob)

CortCVD2.t_rob=ggseg(.data=CortCVD2,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+
  labs(title="CVD2") +  custom_theme
print(CortCVD2.t_rob)
SubcortCVD2.t_rob=ggseg(.data=SubcortCVD2,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+
  #labs(title="CVD2") + 
  custom_theme_sub
print(SubcortCVD2.t_rob)

CortCVD3.t_rob=ggseg(.data=CortCVD3,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+
  labs(title="CVD3") +  custom_theme
print(CortCVD3.t_rob)
SubcortCVD3.t_rob=ggseg(.data=SubcortCVD3,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+
  #labs(title="CVD3") + 
  custom_theme_sub
print(SubcortCVD3.t_rob)

CortCentiloid.t_rob=ggseg(.data=CortCentiloid,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+
  labs(title="Centiloid") +  custom_theme
print(CortCentiloid.t_rob)
SubcortCentiloid.t_rob=ggseg(.data=SubcortCentiloid,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+
  #labs(title="Centiloid") + 
  custom_theme_sub
print(SubcortCentiloid.t_rob)

CortCVD1_Cent.t_rob=ggseg(.data=CortCVD1_Cent,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+
  labs(title="CVD1_Cent") +  custom_theme
print(CortCVD1_Cent.t_rob)
SubcortCVD1_Cent.t_rob=ggseg(.data=SubcortCVD1_Cent,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+
  #labs(title="CVD1") + 
  custom_theme_sub
print(SubcortCVD1_Cent.t_rob)

CortCVD2_Cent.t_rob=ggseg(.data=CortCVD2_Cent,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+
  labs(title="CVD2_Cent") +  custom_theme
print(CortCVD2_Cent.t_rob)
SubcortCVD2_Cent.t_rob=ggseg(.data=SubcortCVD2_Cent,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+
  #labs(title="CVD2") + 
  custom_theme_sub
print(SubcortCVD2_Cent.t_rob)

CortCVD3_Cent.t_rob=ggseg(.data=CortCVD3_Cent,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+
  labs(title="CVD3_Cent") +  custom_theme
print(CortCVD3_Cent.t_rob)
SubcortCVD3_Cent.t_rob=ggseg(.data=SubcortCVD3_Cent,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+
  #labs(title="CVD3") + 
  custom_theme_sub
print(SubcortCVD3_Cent.t_rob)

# only significant ones
CortAge.Sig.t_rob=ggseg(.data=CortAge.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+
  labs(title="Age") +  custom_theme
print(CortAge.Sig.t_rob)
SubcortAge.Sig.t_rob=ggseg(.data=SubcortAge.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+
  #labs(title="Age") + 
  custom_theme_sub
print(SubcortAge.Sig.t_rob)

CortSex.Sig.t_rob=ggseg(.data=CortSex.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+
  labs(title="Sex") +  custom_theme
print(CortSex.Sig.t_rob)
SubcortSex.Sig.t_rob=ggseg(.data=SubcortSex.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+
  #labs(title="Sex") + 
  custom_theme_sub
print(SubcortSex.Sig.t_rob)

CortAPOE4.Sig.t_rob=ggseg(.data=CortAPOE4.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+
  labs(title="APOE4") +  custom_theme
print(CortAPOE4.Sig.t_rob)
SubcortAPOE4.Sig.t_rob=ggseg(.data=SubcortAPOE4.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+
  #labs(title="APOE4") + 
  custom_theme_sub
print(SubcortAPOE4.Sig.t_rob)

Cortfollowup_time.Sig.t_rob=ggseg(.data=Cortfollowup_time.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+
  labs(title="followup_time") +  custom_theme
print(Cortfollowup_time.Sig.t_rob)
Subcortfollowup_time.Sig.t_rob=ggseg(.data=Subcortfollowup_time.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+
  #labs(title="APOE4") + 
  custom_theme_sub
print(Subcortfollowup_time.Sig.t_rob)

CortCVD1.Sig.t_rob=ggseg(.data=CortCVD1.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+
  labs(title="CVD1") +  custom_theme
print(CortCVD1.Sig.t_rob)
SubcortCVD1.Sig.t_rob=ggseg(.data=SubcortCVD1.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+
  #labs(title="CVD1") + 
  custom_theme_sub
print(SubcortCVD1.Sig.t_rob)

CortCVD2.Sig.t_rob=ggseg(.data=CortCVD2.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+
  labs(title="CVD2") +  custom_theme
print(CortCVD2.Sig.t_rob)
SubcortCVD2.Sig.t_rob=ggseg(.data=SubcortCVD2.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+
  #labs(title="CVD2") + 
  custom_theme_sub
print(SubcortCVD2.Sig.t_rob)

CortCVD3.Sig.t_rob=ggseg(.data=CortCVD3.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+
  labs(title="CVD3") +  custom_theme
print(CortCVD3.Sig.t_rob)
SubcortCVD3.Sig.t_rob=ggseg(.data=SubcortCVD3.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+
  #labs(title="CVD3") + 
  custom_theme_sub
print(SubcortCVD3.Sig.t_rob)

CortCentiloid.Sig.t_rob=ggseg(.data=CortCentiloid.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+
  labs(title="Centiloid") +  custom_theme
print(CortCentiloid.Sig.t_rob)
SubcortCentiloid.Sig.t_rob=ggseg(.data=SubcortCentiloid.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+
  #labs(title="Centiloid") + 
  custom_theme_sub
print(SubcortCentiloid.Sig.t_rob)

CortCVD1_Cent.Sig.t_rob=ggseg(.data=CortCVD1_Cent.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+
  labs(title="CVD1_Cent") +  custom_theme
print(CortCVD1_Cent.Sig.t_rob)
SubcortCVD1_Cent.Sig.t_rob=ggseg(.data=SubcortCVD1_Cent.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+
  #labs(title="CVD1") + 
  custom_theme_sub
print(SubcortCVD1_Cent.Sig.t_rob)

CortCVD2_Cent.Sig.t_rob=ggseg(.data=CortCVD2_Cent.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+
  labs(title="CVD2_Cent") +  custom_theme
print(CortCVD2_Cent.Sig.t_rob)
SubcortCVD2_Cent.Sig.t_rob=ggseg(.data=SubcortCVD2_Cent.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+
  #labs(title="CVD2") + 
  custom_theme_sub
print(SubcortCVD2_Cent.Sig.t_rob)

CortCVD3_Cent.Sig.t_rob=ggseg(.data=CortCVD3_Cent.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left")+
  labs(title="CVD3_Cent") +  custom_theme
print(CortCVD3_Cent.Sig.t_rob)
SubcortCVD3_Cent.Sig.t_rob=ggseg(.data=SubcortCVD3_Cent.Sig,mapping=aes(fill=t),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+
  #labs(title="CVD3") + 
  custom_theme_sub
print(SubcortCVD3_Cent.Sig.t_rob)

rm(CortAge,SubcortAge,CortAge.Sig,SubcortAge.Sig,CortSex,SubcortSex,CortSex.Sig,SubcortSex.Sig,CortCentiloid,SubcortCentiloid,CortCentiloid.Sig,SubcortCentiloid.Sig,
    CortAPOE4,SubcortAPOE4,CortAPOE4.Sig,SubcortAPOE4.Sig,CortCVD1,SubcortCVD1,CortCVD1.Sig,SubcortCVD1.Sig,CortCVD2,SubcortCVD2,CortCVD2.Sig,SubcortCVD2.Sig,
   CortCVD3,SubcortCVD3,CortCVD3.Sig,SubcortCVD3.Sig,Cortfollowup_time,Subcortfollowup_time,Cortfollowup_time.Sig,Subcortfollowup_time.Sig,
   CortCVD1_Cent,SubcortCVD1_Cent,CortCVD1_Cent.Sig,SubcortCVD1_Cent.Sig,CortCVD2_Cent,SubcortCVD2_Cent,CortCVD2_Cent.Sig,SubcortCVD2_Cent.Sig,
   CortCVD3_Cent,SubcortCVD3_Cent,CortCVD3_Cent.Sig,SubcortCVD3_Cent.Sig)

# combine plots to get figures
age_row=plot_grid(CortAge.t_rob, SubcortAge.t_rob,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") ) 
sex_row=plot_grid(CortSex.t_rob, SubcortSex.t_rob,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
apoe_row=plot_grid(CortAPOE4.t_rob, SubcortAPOE4.t_rob,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
followup_time_row=plot_grid(Cortfollowup_time.t_rob, Subcortfollowup_time.t_rob,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
cvd1_row=plot_grid(CortCVD1.t_rob, SubcortCVD1.t_rob,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
cvd2_row=plot_grid(CortCVD2.t_rob, SubcortCVD2.t_rob,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
cvd3_row=plot_grid(CortCVD3.t_rob, SubcortCVD3.t_rob,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
Centiloid_row=plot_grid(CortCentiloid.t_rob, SubcortCentiloid.t_rob,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
cvd1_cent_row=plot_grid(CortCVD1_Cent.t_rob, SubcortCVD1_Cent.t_rob,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
cvd2_cent_row=plot_grid(CortCVD2_Cent.t_rob, SubcortCVD2_Cent.t_rob,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
cvd3_cent_row=plot_grid(CortCVD3_Cent.t_rob, SubcortCVD3_Cent.t_rob,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
merged_brain_rob=plot_grid(age_row, sex_row, apoe_row, followup_time_row, cvd1_row, cvd2_row, cvd3_row, Centiloid_row, cvd1_cent_row, cvd2_cent_row, cvd3_cent_row, ncol=1)+theme(plot.background = element_rect(fill = "black"), panel.border = element_blank())
# Save Figure
tiff('./Rabin_Replication_Output/Tau/Colsp_Unthresholded_rob_1.tiff', units="in", width=7, height=9, res=200,type="cairo" )
print(merged_brain_rob)
dev.off()

age_row=plot_grid(CortAge.Sig.t_rob, SubcortAge.Sig.t_rob,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") ) 
sex_row=plot_grid(CortSex.Sig.t_rob, SubcortSex.Sig.t_rob,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
apoe_row=plot_grid(CortAPOE4.Sig.t_rob, SubcortAPOE4.Sig.t_rob,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
followup_time_row=plot_grid(Cortfollowup_time.Sig.t_rob, Subcortfollowup_time.Sig.t_rob,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
cvd1_row=plot_grid(CortCVD1.Sig.t_rob, SubcortCVD1.Sig.t_rob,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
cvd2_row=plot_grid(CortCVD2.Sig.t_rob, SubcortCVD2.Sig.t_rob,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
cvd3_row=plot_grid(CortCVD3.Sig.t_rob, SubcortCVD3.Sig.t_rob,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
Centiloid_row=plot_grid(CortCentiloid.Sig.t_rob, SubcortCentiloid.Sig.t_rob,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
cvd1_cent_row=plot_grid(CortCVD1_Cent.Sig.t_rob, SubcortCVD1_Cent.Sig.t_rob,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
cvd2_cent_row=plot_grid(CortCVD2_Cent.Sig.t_rob, SubcortCVD2_Cent.Sig.t_rob,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
cvd3_cent_row=plot_grid(CortCVD3_Cent.Sig.t_rob, SubcortCVD3_Cent.Sig.t_rob,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
merged_brain.sig_rob=plot_grid(age_row, sex_row, apoe_row, followup_time_row, cvd1_row, cvd2_row, cvd3_row, Centiloid_row, cvd1_cent_row, cvd2_cent_row, cvd3_cent_row, ncol=1)+theme(plot.background = element_rect(fill = "black"), panel.border = element_blank())
# Save Figure
tiff('./Rabin_Replication_Output/Tau/Colsp_Thresholded_rob_1.tiff', units="in", width=7, height=9, res=200,type="cairo")
print(merged_brain.sig_rob)
dev.off()

rm(age_row, sex_row, Centiloid_row, apoe_row, cvd1_row, cvd2_row, cvd3_row, 
   CortAge.t,SubcortAge.t,CortAge.Sig.t,SubcortAge.Sig.t,CortSex.t,SubcortSex.t,CortSex.Sig.t,SubcortSex.Sig.t,CortCentiloid.t,SubcortCentiloid.t,CortCentiloid.Sig.t,SubcortCentiloid.Sig.t,
   CortAPOE4.t,SubcortAPOE4.t,CortAPOE4.Sig.t,SubcortAPOE4.Sig.t, CortCVD1.t,SubcortCVD1.t,CortCVD1.Sig.t,SubcortCVD1.Sig.t,CortCVD2.t,SubcortCVD2.t,CortCVD2.Sig.t,SubcortCVD2.Sig.t,
   CortCVD3.t,SubcortCVD3.t,CortCVD3.Sig.t,SubcortCVD3.Sig.t,CortCVD1_Cent.t,SubcortCVD1_Cent.t,CortCVD1_Cent.Sig.t,SubcortCVD1_Cent.Sig.t,CortCVD2_Cent.t,SubcortCVD2_Cent.t,CortCVD2_Cent.Sig.t,SubcortCVD2_Cent.Sig.t,
   CortCVD3_Cent.t,SubcortCVD3_Cent.t,CortCVD3_Cent.Sig.t,SubcortCVD3_Cent.Sig.t,CortAge.t_rob,SubcortAge.t_rob,CortAge.Sig.t_rob,SubcortAge.Sig.t_rob,CortSex.t_rob,SubcortSex.t_rob,CortSex.Sig.t_rob,SubcortSex.Sig.t_rob,
   CortCentiloid.t_rob,SubcortCentiloid.t_rob,CortCentiloid.Sig.t_rob,SubcortCentiloid.Sig.t_rob,CortAPOE4.t_rob,SubcortAPOE4.t_rob,CortAPOE4.Sig.t_rob,SubcortAPOE4.Sig.t_rob,
   CortCVD1.t_rob,SubcortCVD1.t_rob,CortCVD1.Sig.t_rob,SubcortCVD1.Sig.t_rob,CortCVD2.t_rob,SubcortCVD2.t_rob,CortCVD2.Sig.t_rob,SubcortCVD2.Sig.t_rob,
   CortCVD3.t_rob,SubcortCVD3.t_rob,CortCVD3.Sig.t_rob,SubcortCVD3.Sig.t_rob,merged_brain_rob,merged_brain.sig_rob,model,
   CortCVD1_Cent.t_rob,SubcortCVD1_Cent.t_rob,CortCVD1_Cent.Sig.t_rob,SubcortCVD1_Cent.Sig.t_rob,CortCVD2_Cent.t_rob,SubcortCVD2_Cent.t_rob,CortCVD2_Cent.Sig.t_rob,SubcortCVD2_Cent.Sig.t_rob,
   CortCVD3_Cent.t_rob,SubcortCVD3_Cent.t_rob,CortCVD3_Cent.Sig.t_rob,SubcortCVD3_Cent.Sig.t_rob,merged_brain_rob,merged_brain.sig_rob,model,
   Cortfollowup_time.t,Subcortfollowup_time.t,Cortfollowup_time.Sig.t,Subcortfollowup_time.Sig.t,)

