# Title: ADRC - PiB Regression Loop 
# Author: Diana Hobbs
# Date: April 2022

#Functions to set up environment
rm(list=ls()) #removes all variables from the current environment
options(scipen = 999) #forces R not to use scientific notation

#### Set Up ####
# Load packages 
pacman::p_load(tidyverse,jtools,naniar,lubridate,lme4,sfsmisc,foreign,MASS,scales,ggseg,cowplot,GeneNet)

###### LOAD & SETUP DATA 
PiB <- read.csv("./Data/Rabin_Replication_PIB.csv")%>%
  mutate(CVD = str_replace_all(CVD_Colsp, c("<5%"="0", "5-10%"="1", "10-20%"="2",">20%"="3")))


# Reassign variable types & set reference 
PiB$SEX <- as.factor(PiB$SEX); PiB$SEX <- relevel(PiB$SEX, "Male") 
PiB$APOE4 <- as.factor(PiB$APOE4); PiB$APOE4 <- relevel(PiB$APOE4, "0")
PiB$CVD <- as.factor(PiB$CVD); PiB$CVD <- relevel(PiB$CVD, "0")

# Z transforming continuous variables 
#test<-PiB%>%
#  mutate(BMI_z = scale(BMI))

###### REGRESSION LOOP
ROI <- PiB[,20:60]

# Loop through regressions and keep output 
LM.COF<-list() # output for regression
RLM.COF<-list() # output for robust regression
model<- matrix(0, 41, 7)
outcount=1
# Create a 'for' loop 
for (i in 1:ncol(ROI)){
  lm <- lm(ROI[,i] ~ AGE + SEX + APOE4 + CVD, data=PiB)
  rlm <- rlm(ROI[,i] ~ AGE + SEX + APOE4 + CVD, data=PiB)
  lm.sum <- summary(lm)
  rlm.sum <- summary(rlm)
  lm.cof <- lm.sum$coefficients
  rlm.cof <- rlm.sum$coefficients
  region=colnames(ROI)[i]
  model[outcount,1]=region
  colnames(model) <- c("region", "p.adjust.age", "p.adjust.sex", "p.adjust.apoe","p.adjust.CVD1", "p.adjust.CVD2", "p.adjust.CVD3")
  rob_p=f.robftest(rlm, var = "AGE"); model[outcount, 2]=rob_p$p.value 
  rob_p=f.robftest(rlm, var = "SEXFemale"); model[outcount, 3]=rob_p$p.value
  rob_p=f.robftest(rlm, var = "APOE41");model[outcount, 4]=rob_p$p.value 
  rob_p=f.robftest(rlm, var = "CVD1");model[outcount, 5]=rob_p$p.value 
  rob_p=f.robftest(rlm, var = "CVD2");model[outcount, 6]=rob_p$p.value 
  rob_p=f.robftest(rlm, var = "CVD3");model[outcount, 7]=rob_p$p.value 
  
  outcount=outcount+1
  LM.COF[[i]] <- lm.cof
  RLM.COF[[i]] <- rlm.cof
  names(LM.COF)[i] <- colnames(ROI)[i]
  names(RLM.COF)[i] <- colnames(ROI)[i]
}
write.csv(LM.COF,"./Rabin_Replication_Output/PiB/Colsp_LM.COF_1.csv")
write.csv(RLM.COF,"./Rabin_Replication_Output/PiB/Colsp_RLM.COF_1.csv")
rm(lm,lm.cof,LM.COF,lm.sum,rlm,rlm.cof,RLM.COF,rlm.sum,i,ROI,rob_p,outcount,PiB,region)

###### Organize loop output in usable format 
Data <- read.csv("./Rabin_Replication_Output/PiB/Colsp_LM.COF_1.csv")

# Extract Estimates 
Estimate <- Data[,grepl(".Estimate", names(Data))]
names(Estimate) = gsub(pattern = ".Estimate", replacement = "", x = names(Estimate))
Estimate <- t(Estimate)
colnames(Estimate) = c("(Intercept)","AGE","SEX","APOE4","CVD1","CVD2","CVD3")
Estimate <- data.frame(Estimate)

TVal <- Data[,grepl(".t.value", names(Data))]
names(TVal) = gsub(pattern = ".t.value", replacement = "", x = names(TVal))
TVal <- t(TVal)
colnames(TVal) <- c("(Intercept)","AGE","SEX","APOE4","CVD1","CVD2","CVD3")
TVal <- data.frame(TVal)

PVal <- Data[,grepl(".Pr...t..", names(Data))]
names(PVal) = gsub(pattern = ".Pr...t..", replacement = "", x = names(PVal))
PVal <- t(PVal)
colnames(PVal) <- c("(Intercept)","AGE","SEX","APOE4","CVD1","CVD2","CVD3")
PVal <- data.frame(PVal)

# Merge coefficient summaries & Perform FDR correction using Benjamini and Hochberg method
Age <- data.frame("B"=Estimate$AGE,"t"=TVal$AGE,"p"=PVal$AGE,"p_corr"=(p.adjust(PVal$AGE, method = "BH",n = length(PVal$AGE))))
Sex <- data.frame("B"=Estimate$SEX,"t"=TVal$SEX,"p"=PVal$SEX,"p_corr"=(p.adjust(PVal$SEX, method = "BH",n = length(PVal$SEX))))
APOE4 <- data.frame("B"=Estimate$APOE4,"t"=TVal$APOE4,"p"=PVal$APOE4,"p_corr"=(p.adjust(PVal$APOE4, method = "BH",n = length(PVal$APOE4))))
CVD1 <- data.frame("B"=Estimate$CVD1,"t"=TVal$CVD1,"p"=PVal$CVD1,"p_corr"=(p.adjust(PVal$CVD1, method = "BH",n = length(PVal$CVD1))))
CVD2 <- data.frame("B"=Estimate$CVD2,"t"=TVal$CVD2,"p"=PVal$CVD2,"p_corr"=(p.adjust(PVal$CVD2, method = "BH",n = length(PVal$CVD2))))
CVD3 <- data.frame("B"=Estimate$CVD3,"t"=TVal$CVD3,"p"=PVal$CVD3,"p_corr"=(p.adjust(PVal$CVD3, method = "BH",n = length(PVal$CVD3))))
rm(Data,Estimate,TVal,PVal)

cort_list <- read.csv("./BrianFiles/ggseg_namelist.csv")
subcort_list <- cort_list$subcortical[1:7] #extract subcortical list

CortAge <- data.frame(region=cort_list$cortical,filter(Age[1:34,])); SubcortAge <- data.frame(region=subcort_list,filter(Age[35:41,]))
CortSex <- data.frame(region=cort_list$cortical,filter(Sex[1:34,])); SubcortSex <- data.frame(region=subcort_list,filter(Sex[35:41,]))
CortAPOE4 <- data.frame(region=cort_list$cortical,filter(APOE4[1:34,])); SubcortAPOE4 <- data.frame(region=subcort_list,filter(APOE4[35:41,]))
CortCVD1 <- data.frame(region=cort_list$cortical,filter(CVD1[1:34,])); SubcortCVD1 <- data.frame(region=subcort_list,filter(CVD1[35:41,]))
CortCVD2 <- data.frame(region=cort_list$cortical,filter(CVD2[1:34,])); SubcortCVD2 <- data.frame(region=subcort_list,filter(CVD2[35:41,]))
CortCVD3 <- data.frame(region=cort_list$cortical,filter(CVD3[1:34,])); SubcortCVD3 <- data.frame(region=subcort_list,filter(CVD3[35:41,]))
rm(Age,Sex,APOE4,CVD1,CVD2,CVD3)

CortAge.Sig <- subset(CortAge, p_corr < 0.05); SubcortAge.Sig <- subset(SubcortAge, p_corr < 0.05)
CortSex.Sig <- subset(CortSex, p_corr < 0.05); SubcortSex.Sig <- subset(SubcortSex, p_corr < 0.05)
CortAPOE4.Sig <- subset(CortAPOE4, p_corr < 0.05); SubcortAPOE4.Sig <- subset(SubcortAPOE4, p_corr < 0.05)
CortCVD1.Sig <- subset(CortCVD1, p_corr < 0.05); SubcortCVD1.Sig <- subset(SubcortCVD1, p_corr < 0.05)
CortCVD2.Sig <- subset(CortCVD2, p_corr < 0.05); SubcortCVD2.Sig <- subset(SubcortCVD2, p_corr < 0.05)
CortCVD3.Sig <- subset(CortCVD3, p_corr < 0.05); SubcortCVD3.Sig <- subset(SubcortCVD3, p_corr < 0.05)

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

rm(CortAge,SubcortAge,CortAge.Sig,SubcortAge.Sig,CortSex,SubcortSex,CortSex.Sig,SubcortSex.Sig,CortAPOE4,SubcortAPOE4,CortAPOE4.Sig,SubcortAPOE4.Sig,
   CortCVD1,SubcortCVD1,CortCVD1.Sig,SubcortCVD1.Sig,CortCVD2,SubcortCVD2,CortCVD2.Sig,SubcortCVD2.Sig,CortCVD3,SubcortCVD3,CortCVD3.Sig,SubcortCVD3.Sig)

# combine plots to get figures
age_row=plot_grid(CortAge.t, SubcortAge.t,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") ) 
sex_row=plot_grid(CortSex.t, SubcortSex.t,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
apoe_row=plot_grid(CortAPOE4.t, SubcortAPOE4.t,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
cvd1_row=plot_grid(CortCVD1.t, SubcortCVD1.t,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
cvd2_row=plot_grid(CortCVD2.t, SubcortCVD2.t,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
cvd3_row=plot_grid(CortCVD3.t, SubcortCVD3.t,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
merged_brain=plot_grid(age_row, sex_row, apoe_row, cvd1_row, cvd2_row, cvd3_row, ncol=1)+theme(plot.background = element_rect(fill = "black"), panel.border = element_blank())
# Save Figure
tiff('./Rabin_Replication_Output/PiB/Colsp_Unthresholded_1.tiff', units="in", width=7, height=9, res=200,type="cairo" )
print(merged_brain)
dev.off()

age_row=plot_grid(CortAge.Sig.t, SubcortAge.Sig.t,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") ) 
sex_row=plot_grid(CortSex.Sig.t, SubcortSex.Sig.t,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
apoe_row=plot_grid(CortAPOE4.Sig.t, SubcortAPOE4.Sig.t,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
cvd1_row=plot_grid(CortCVD1.Sig.t, SubcortCVD1.Sig.t,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
cvd2_row=plot_grid(CortCVD2.Sig.t, SubcortCVD2.Sig.t,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
cvd3_row=plot_grid(CortCVD3.Sig.t, SubcortCVD3.Sig.t,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
merged_brain.sig=plot_grid(age_row, sex_row, apoe_row, cvd1_row, cvd2_row, cvd3_row, ncol=1)+theme(plot.background = element_rect(fill = "black"), panel.border = element_blank())
# Save Figure
tiff('./Rabin_Replication_Output/PiB/Colsp_Thresholded_1.tiff', units="in", width=7, height=9, res=200,type="cairo")
print(merged_brain.sig)
dev.off()
rm(age_row, sex_row, apoe_row, cvd1_row, cvd2_row, cvd3_row, merged_brain,merged_brain.sig)

###### Robust regression ######
Data <- read.csv("./Rabin_Replication_Output/PiB/Colsp_RLM.COF_1.csv")
# Extract Estimates 
Estimate <- Data[,grepl(".Value", names(Data))]
names(Estimate) = gsub(pattern = ".Estimate", replacement = "", x = names(Estimate))
Estimate <- t(Estimate)
colnames(Estimate) = c("(Intercept)","AGE","SEX","APOE4","CVD1","CVD2","CVD3")
Estimate <- data.frame(Estimate)

TVal <- Data[,grepl(".t.value", names(Data))]
names(TVal) = gsub(pattern = ".t.value", replacement = "", x = names(TVal))
TVal <- t(TVal)
colnames(TVal) <- c("(Intercept)","AGE","SEX","APOE4","CVD1","CVD2","CVD3")
TVal <- data.frame(TVal)

PVal <- data.frame(model)
names(PVal) = gsub(pattern = "P.adjust.", replacement = "", x = names(PVal))
colnames(PVal) <- c("region","AGE","SEX","APOE4","CVD1","CVD2","CVD3")

# Merge coefficient summaries 
Age <- data.frame("B"=Estimate$AGE,"t"=TVal$AGE,"p"=PVal$AGE,"p_corr"=(p.adjust(PVal$AGE, method = "BH",n = length(PVal$AGE))))
Sex <- data.frame("B"=Estimate$SEX,"t"=TVal$SEX,"p"=PVal$SEX,"p_corr"=(p.adjust(PVal$SEX, method = "BH",n = length(PVal$SEX))))
APOE4 <- data.frame("B"=Estimate$APOE4,"t"=TVal$APOE4,"p"=PVal$APOE4,"p_corr"=(p.adjust(PVal$APOE4, method = "BH",n = length(PVal$APOE4))))
CVD1 <- data.frame("B"=Estimate$CVD1,"t"=TVal$CVD1,"p"=PVal$CVD1,"p_corr"=(p.adjust(PVal$CVD1, method = "BH",n = length(PVal$CVD1))))
CVD2 <- data.frame("B"=Estimate$CVD2,"t"=TVal$CVD2,"p"=PVal$CVD2,"p_corr"=(p.adjust(PVal$CVD2, method = "BH",n = length(PVal$CVD2))))
CVD3 <- data.frame("B"=Estimate$CVD3,"t"=TVal$CVD3,"p"=PVal$CVD3,"p_corr"=(p.adjust(PVal$CVD3, method = "BH",n = length(PVal$CVD3))))
rm(Data,Estimate,TVal,PVal)

CortAge <- data.frame(region=cort_list$cortical,filter(Age[1:34,])); SubcortAge <- data.frame(region=subcort_list,filter(Age[35:41,]))
CortSex <- data.frame(region=cort_list$cortical,filter(Sex[1:34,])); SubcortSex <- data.frame(region=subcort_list,filter(Sex[35:41,]))
CortAPOE4 <- data.frame(region=cort_list$cortical,filter(APOE4[1:34,])); SubcortAPOE4 <- data.frame(region=subcort_list,filter(APOE4[35:41,]))
CortCVD1 <- data.frame(region=cort_list$cortical,filter(CVD1[1:34,])); SubcortCVD1 <- data.frame(region=subcort_list,filter(CVD1[35:41,]))
CortCVD2 <- data.frame(region=cort_list$cortical,filter(CVD2[1:34,])); SubcortCVD2 <- data.frame(region=subcort_list,filter(CVD2[35:41,]))
CortCVD3 <- data.frame(region=cort_list$cortical,filter(CVD3[1:34,])); SubcortCVD3 <- data.frame(region=subcort_list,filter(CVD3[35:41,]))
rm(Age,Sex,APOE4,CVD1,CVD2,CVD3)

CortAge.Sig <- subset(CortAge, p_corr < 0.05); SubcortAge.Sig <- subset(SubcortAge, p_corr < 0.05)
CortSex.Sig <- subset(CortSex, p_corr < 0.05); SubcortSex.Sig <- subset(SubcortSex, p_corr < 0.05)
CortAPOE4.Sig <- subset(CortAPOE4, p_corr < 0.05); SubcortAPOE4.Sig <- subset(SubcortAPOE4, p_corr < 0.05)
CortCVD1.Sig <- subset(CortCVD1, p_corr < 0.05); SubcortCVD1.Sig <- subset(SubcortCVD1, p_corr < 0.05)
CortCVD2.Sig <- subset(CortCVD2, p_corr < 0.05); SubcortCVD2.Sig <- subset(SubcortCVD2, p_corr < 0.05)
CortCVD3.Sig <- subset(CortCVD3, p_corr < 0.05); SubcortCVD3.Sig <- subset(SubcortCVD3, p_corr < 0.05)

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

rm(CortAge,SubcortAge,CortAge.Sig,SubcortAge.Sig,CortSex,SubcortSex,CortSex.Sig,SubcortSex.Sig,CortAPOE4,SubcortAPOE4,CortAPOE4.Sig,SubcortAPOE4.Sig,
   CortCVD1,SubcortCVD1,CortCVD1.Sig,SubcortCVD1.Sig,CortCVD2,SubcortCVD2,CortCVD2.Sig,SubcortCVD2.Sig,   CortCVD3,SubcortCVD3,CortCVD3.Sig,SubcortCVD3.Sig)

# combine plots to get figures
age_row=plot_grid(CortAge.t_rob, SubcortAge.t_rob,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") ) 
sex_row=plot_grid(CortSex.t_rob, SubcortSex.t_rob,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
apoe_row=plot_grid(CortAPOE4.t_rob, SubcortAPOE4.t_rob,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
cvd1_row=plot_grid(CortCVD1.t_rob, SubcortCVD1.t_rob,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
cvd2_row=plot_grid(CortCVD2.t_rob, SubcortCVD2.t_rob,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
cvd3_row=plot_grid(CortCVD3.t_rob, SubcortCVD3.t_rob,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
merged_brain_rob=plot_grid(age_row, sex_row, apoe_row,cvd1_row, cvd2_row, cvd3_row,ncol=1)+theme(plot.background = element_rect(fill = "black"), panel.border = element_blank())
# Save Figure
tiff('./Rabin_Replication_Output/PiB/Colsp_Unthresholded_rob_1.tiff', units="in", width=7, height=9, res=200,type="cairo" )
print(merged_brain_rob)
dev.off()

age_row=plot_grid(CortAge.Sig.t_rob, SubcortAge.Sig.t_rob,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") ) 
sex_row=plot_grid(CortSex.Sig.t_rob, SubcortSex.Sig.t_rob,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
apoe_row=plot_grid(CortAPOE4.Sig.t_rob, SubcortAPOE4.Sig.t_rob,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
cvd1_row=plot_grid(CortCVD1.Sig.t_rob, SubcortCVD1.Sig.t_rob,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
cvd2_row=plot_grid(CortCVD2.Sig.t_rob, SubcortCVD2.Sig.t_rob,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
cvd3_row=plot_grid(CortCVD3.Sig.t_rob, SubcortCVD3.Sig.t_rob,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
merged_brain.sig_rob=plot_grid(age_row, sex_row, apoe_row, cvd1_row, cvd2_row, cvd3_row, ncol=1)+theme(plot.background = element_rect(fill = "black"), panel.border = element_blank())
# Save Figure
tiff('./Rabin_Replication_Output/PiB/Colsp_Thresholded_rob_1.tiff', units="in", width=7, height=9, res=200,type="cairo")
print(merged_brain.sig_rob)
dev.off()

rm(age_row, sex_row, Centiloid_row, apoe_row, cvd1_row, cvd2_row, cvd3_row,CortAge.t,SubcortAge.t,CortAge.Sig.t,SubcortAge.Sig.t,CortSex.t,SubcortSex.t,CortSex.Sig.t,SubcortSex.Sig.t,
   CortAPOE4.t,SubcortAPOE4.t,CortAPOE4.Sig.t,SubcortAPOE4.Sig.t, CortCVD1.t,SubcortCVD1.t,CortCVD1.Sig.t,SubcortCVD1.Sig.t,CortCVD2.t,SubcortCVD2.t,CortCVD2.Sig.t,SubcortCVD2.Sig.t,
   CortCVD3.t,SubcortCVD3.t,CortCVD3.Sig.t,SubcortCVD3.Sig.t,CortAge.t_rob,SubcortAge.t_rob,CortAge.Sig.t_rob,SubcortAge.Sig.t_rob,CortSex.t_rob,SubcortSex.t_rob,CortSex.Sig.t_rob,SubcortSex.Sig.t_rob,
   CortAPOE4.t_rob,SubcortAPOE4.t_rob,CortAPOE4.Sig.t_rob,SubcortAPOE4.Sig.t_rob,CortCVD1.t_rob,SubcortCVD1.t_rob,CortCVD1.Sig.t_rob,SubcortCVD1.Sig.t_rob,CortCVD2.t_rob,SubcortCVD2.t_rob,CortCVD2.Sig.t_rob,SubcortCVD2.Sig.t_rob,
   CortCVD3.t_rob,SubcortCVD3.t_rob,CortCVD3.Sig.t_rob,SubcortCVD3.Sig.t_rob,merged_brain_rob,merged_brain.sig_rob,model)

