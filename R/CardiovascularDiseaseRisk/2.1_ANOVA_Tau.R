# Title: ADRC - Tau Regression Loop ANOVA
# Author: Diana Hobbs
# Date: April 2022

#### Set Up ####
# Load packages 
pacman::p_load(tidyverse,jtools,naniar,lubridate,lme4,sfsmisc,foreign,MASS,scales,ggseg,cowplot)

###### LOAD & SETUP DATA 
tau <- read.csv("./Data/Rabin_Replication_TAU.csv")%>%
  mutate(CVD = str_replace_all(CVD_Colsp, c("<5%"="0", "5-10%"="1", "10-20%"="2",">20%"="3")))

# Reassign variable types & set reference 
tau$SEX <- as.factor(tau$SEX); tau$SEX <- relevel(tau$SEX, "Male") 
tau$APOE4 <- as.factor(tau$APOE4); tau$APOE4 <- relevel(tau$APOE4, "0")
tau$CVD <- as.factor(tau$CVD); tau$CVD <- relevel(tau$CVD, "0")

###### REGRESSION LOOP
ROI <- tau[,20:60]

# Loop through regressions and keep output 
MODEL <- list() # output for ANOVA 

# Create a 'for' loop 
for (i in 1:ncol(ROI)){
  model <- anova(lm(ROI[,i] ~ AGE + SEX + APOE4 + followup_time + CVD + Centiloid + CVD*Centiloid, data=tau))
  MODEL[[i]] <- model
  names(MODEL)[i] <- colnames(ROI)[i]
}
write.csv(MODEL,"./Rabin_Replication_Output/Tau/Colsp_Anova_1.csv")
rm(model,MODEL,i,tau,ROI)

###### Organize loop output in usable format 
Data <- read.csv("./Rabin_Replication_Output/Tau/Colsp_Anova_1.csv")

# Extract F and p-values 
Fvalue <- Data[ ,grepl("F.value", names(Data))]
names(Fvalue) = gsub(pattern = ".F.value", replacement = "", x = names(Fvalue))
Fvalue <- t(Fvalue)
colnames(Fvalue) = c("AGE","SEX","APOE4","followup_time", "CVD","Centiloid","CVD_Cent","Residuals")
Fvalue <- as.data.frame(Fvalue)

Pvalue <- Data[ ,grepl("Pr..F.", names(Data))]
names(Pvalue) = gsub(pattern = "Pr..F.", replacement = "", x = names(Pvalue))
Pvalue <- t(Pvalue)
colnames(Pvalue) = c("AGE","SEX","APOE4","followup_time","CVD","Centiloid","CVD_Cent","Residuals")
Pvalue <- as.data.frame(Pvalue)

# Merge coefficient summaries & Perform FDR correction using Benjamini and Hochberg method
Age <- data.frame("f"=Fvalue$AGE,"p"=Pvalue$AGE,"p_corr"=(p.adjust(Pvalue$AGE, method = "BH",n = length(Pvalue$AGE))))
Sex <- data.frame("f"=Fvalue$SEX,"p"=Pvalue$SEX,"p_corr"=(p.adjust(Pvalue$SEX, method = "BH",n = length(Pvalue$SEX))))
APOE4 <- data.frame("f"=Fvalue$APOE4,"p"=Pvalue$APOE4,"p_corr"=(p.adjust(Pvalue$APOE4, method = "BH",n = length(Pvalue$APOE4))))
followup_time <- data.frame("f"=Fvalue$followup_time,"p"=Pvalue$followup_time,"p_corr"=(p.adjust(Pvalue$followup_time, method = "BH",n = length(Pvalue$followup_time))))
Centiloid <- data.frame("f"=Fvalue$Centiloid,"p"=Pvalue$Centiloid,"p_corr"=(p.adjust(Pvalue$Centiloid, method = "BH",n = length(Pvalue$Centiloid))))
CVD <- data.frame("f"=Fvalue$CVD,"p"=Pvalue$CVD,"p_corr"=(p.adjust(Pvalue$CVD, method = "BH",n = length(Pvalue$CVD))))
CVD_Cent <- data.frame("f"=Fvalue$CVD_Cent,"p"=Pvalue$CVD_Cent,"p_corr"=(p.adjust(Pvalue$CVD_Cent, method = "BH",n = length(Pvalue$CVD_Cent))))
rm(Data,Fvalue,Pvalue)

cort_list <- read.csv("./BrianFiles/ggseg_namelist.csv")
subcort_list <- cort_list$subcortical[1:7] #extract subcortical list

CortAge <- data.frame(region=cort_list$cortical,filter(Age[1:34,])); SubcortAge <- data.frame(region=subcort_list,filter(Age[35:41,]))
CortSex <- data.frame(region=cort_list$cortical,filter(Sex[1:34,])); SubcortSex <- data.frame(region=subcort_list,filter(Sex[35:41,]))
CortAPOE4 <- data.frame(region=cort_list$cortical,filter(APOE4[1:34,])); SubcortAPOE4 <- data.frame(region=subcort_list,filter(APOE4[35:41,]))
Cortfollowup_time <- data.frame(region=cort_list$cortical,filter(followup_time[1:34,])); Subcortfollowup_time <- data.frame(region=subcort_list,filter(followup_time[35:41,]))
CortCentiloid <- data.frame(region=cort_list$cortical,filter(Centiloid[1:34,])); SubCortCentiloid <- data.frame(region=subcort_list,filter(Centiloid[35:41,]))
CortCVD <- data.frame(region=cort_list$cortical,filter(CVD[1:34,])); SubcortCVD <- data.frame(region=subcort_list,filter(CVD[35:41,]))
CortCVD_Cent <- data.frame(region=cort_list$cortical,filter(CVD_Cent[1:34,])); SubcortCVD_Cent <- data.frame(region=subcort_list,filter(CVD_Cent[35:41,]))
rm(Age,Sex,APOE4,followup_time,CVD,CVD_Cent)

CortAge.Sig <- subset(CortAge, p_corr < 0.05); SubcortAge.Sig <- subset(SubcortAge, p_corr < 0.05)
CortSex.Sig <- subset(CortSex, p_corr < 0.05); SubcortSex.Sig <- subset(SubcortSex, p_corr < 0.05)
CortAPOE4.Sig <- subset(CortAPOE4, p_corr < 0.05); SubcortAPOE4.Sig <- subset(SubcortAPOE4, p_corr < 0.05)
Cortfollowup_time.Sig <- subset(Cortfollowup_time, p_corr < 0.05); Subcortfollowup_time.Sig <- subset(Subcortfollowup_time, p_corr < 0.05)
CortCentiloid.Sig <- subset(CortCentiloid, p_corr < 0.05); SubCortCentiloid.Sig <- subset(SubCortCentiloid, p_corr < 0.05)
CortCVD.Sig <- subset(CortCVD, p_corr < 0.05); SubcortCVD.Sig <- subset(SubcortCVD, p_corr < 0.05)
CortCVD_Cent.Sig <- subset(CortCVD_Cent, p_corr < 0.05); SubcortCVD_Cent.Sig <- subset(SubcortCVD_Cent, p_corr < 0.05)

###### Make images 
# This theme specifies parameters for all of the plots. Change here to change any of them
custom_theme = list(
  scale_fill_gradientn(colours=c("#0E5052","#4EB99E","#DFDDAE","#E5903E","#CE2E1C"), limits=c(0,10), oob=squish),
  #use scale_fill_gradient2 and the option midpoint to specify white as another option
  theme(text = element_text(colour = "white")),  theme(axis.text = element_blank(), axis.title = element_blank()), 
  theme(plot.background = element_rect(fill = "black",colour = NA)), theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), theme(panel.border = element_blank()), guides(fill="none", color="none"))

custom_theme_sub = list(
  scale_fill_gradientn(colours=c("#0E5052","#4EB99E","#DFDDAE","#E5903E","#CE2E1C"), limits=c(0,10), oob=squish),
  theme(text = element_text(colour = "white")),  theme(axis.text = element_blank(), axis.title = element_blank()), 
  theme(plot.background = element_rect(fill = "black",colour = NA)), theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), theme(panel.border = element_blank()))

# Plotting f-values on the brain atlas
CortAge.f=ggseg(.data=CortAge,mapping=aes(fill=f),position = "stacked", colour="black",size=.5, hemisphere="left")+
  labs(title="Age") +  custom_theme
print(CortAge.f)
SubcortAge.f=ggseg(.data=SubcortAge,mapping=aes(fill=f),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+
  #labs(title="Age") + 
  custom_theme_sub
print(SubcortAge.f)

CortSex.f=ggseg(.data=CortSex,mapping=aes(fill=f),position = "stacked", colour="black",size=.5, hemisphere="left")+
  labs(title="Sex") +  custom_theme
print(CortSex.f)
SubcortSex.f=ggseg(.data=SubcortSex,mapping=aes(fill=f),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+
  #labs(title="Sex") + 
  custom_theme_sub
print(SubcortSex.f)

CortAPOE4.f=ggseg(.data=CortAPOE4,mapping=aes(fill=f),position = "stacked", colour="black",size=.5, hemisphere="left")+
  labs(title="APOE4") +  custom_theme
print(CortAPOE4.f)
SubcortAPOE4.f=ggseg(.data=SubcortAPOE4,mapping=aes(fill=f),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+
  #labs(title="APOE4") + 
  custom_theme_sub
print(SubcortAPOE4.f)

Cortfollowup_time.f=ggseg(.data=Cortfollowup_time,mapping=aes(fill=f),position = "stacked", colour="black",size=.5, hemisphere="left")+
  labs(title="followup_time") +  custom_theme
print(Cortfollowup_time.f)
Subcortfollowup_time.f=ggseg(.data=Subcortfollowup_time,mapping=aes(fill=f),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+
  #labs(title="APOE4") + 
  custom_theme_sub
print(Subcortfollowup_time.f)

CortCentiloid.f=ggseg(.data=CortCentiloid,mapping=aes(fill=f),position = "stacked", colour="black",size=.5, hemisphere="left")+
  labs(title="Centiloid") +  custom_theme
print(CortCentiloid.f)
SubCortCentiloid.f=ggseg(.data=SubCortCentiloid,mapping=aes(fill=f),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+
  #labs(title="APOE4") + 
  custom_theme_sub
print(SubCortCentiloid.f)

CortCVD.f=ggseg(.data=CortCVD,mapping=aes(fill=f),position = "stacked", colour="black",size=.5, hemisphere="left")+
  labs(title="CVD") +  custom_theme
print(CortCVD.f)
SubcortCVD.f=ggseg(.data=SubcortCVD,mapping=aes(fill=f),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+
  #labs(title="CVD1") + 
  custom_theme_sub
print(SubcortCVD.f)

CortCVD_Cent.f=ggseg(.data=CortCVD_Cent,mapping=aes(fill=f),position = "stacked", colour="black",size=.5, hemisphere="left")+
  labs(title="CVD_Cent") +  custom_theme
print(CortCVD_Cent.f)
SubcortCVD_Cent.f=ggseg(.data=SubcortCVD_Cent,mapping=aes(fill=f),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+
  #labs(title="CVD1") + 
  custom_theme_sub
print(SubcortCVD_Cent.f)

# only significant ones
CortAge.Sig.f=ggseg(.data=CortAge.Sig,mapping=aes(fill=f),position = "stacked", colour="black",size=.5, hemisphere="left")+
  labs(title="Age") +  custom_theme
print(CortAge.Sig.f)
SubcortAge.Sig.f=ggseg(.data=SubcortAge.Sig,mapping=aes(fill=f),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+
  #labs(title="Age") + 
  custom_theme_sub
print(SubcortAge.Sig.f)

CortSex.Sig.f=ggseg(.data=CortSex.Sig,mapping=aes(fill=f),position = "stacked", colour="black",size=.5, hemisphere="left")+
  labs(title="Sex") +  custom_theme
print(CortSex.Sig.f)
SubcortSex.Sig.f=ggseg(.data=SubcortSex.Sig,mapping=aes(fill=f),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+
  #labs(title="Sex") + 
  custom_theme_sub
print(SubcortSex.Sig.f)

CortAPOE4.Sig.f=ggseg(.data=CortAPOE4.Sig,mapping=aes(fill=f),position = "stacked", colour="black",size=.5, hemisphere="left")+
  labs(title="APOE4") +  custom_theme
print(CortAPOE4.Sig.f)
SubcortAPOE4.Sig.f=ggseg(.data=SubcortAPOE4.Sig,mapping=aes(fill=f),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+
  #labs(title="APOE4") + 
  custom_theme_sub
print(SubcortAPOE4.Sig.f)

Cortfollowup_time.Sig.f=ggseg(.data=Cortfollowup_time.Sig,mapping=aes(fill=f),position = "stacked", colour="black",size=.5, hemisphere="left")+
  labs(title="followup_time") +  custom_theme
print(Cortfollowup_time.Sig.f)
Subcortfollowup_time.Sig.f=ggseg(.data=Subcortfollowup_time.Sig,mapping=aes(fill=f),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+
  #labs(title="APOE4") + 
  custom_theme_sub
print(Subcortfollowup_time.Sig.f)

CortCentiloid.Sig.f=ggseg(.data=CortCentiloid.Sig,mapping=aes(fill=f),position = "stacked", colour="black",size=.5, hemisphere="left")+
  labs(title="Centiloid") +  custom_theme
print(CortCentiloid.Sig.f)
SubCortCentiloid.Sig.f=ggseg(.data=SubCortCentiloid.Sig,mapping=aes(fill=f),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+
  #labs(title="APOE4") + 
  custom_theme_sub
print(SubCortCentiloid.Sig.f)

CortCVD.Sig.f=ggseg(.data=CortCVD.Sig,mapping=aes(fill=f),position = "stacked", colour="black",size=.5, hemisphere="left")+
  labs(title="CVD") +  custom_theme
print(CortCVD.Sig.f)
SubcortCVD.Sig.f=ggseg(.data=SubcortCVD.Sig,mapping=aes(fill=f),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+
  #labs(title="CVD1") + 
  custom_theme_sub
print(SubcortCVD.Sig.f)

CortCVD_Cent.Sig.f=ggseg(.data=CortCVD_Cent.Sig,mapping=aes(fill=f),position = "stacked", colour="black",size=.5, hemisphere="left")+
  labs(title="CVD_Cent") +  custom_theme
print(CortCVD_Cent.Sig.f)
SubcortCVD_Cent.Sig.f=ggseg(.data=SubcortCVD_Cent.Sig,mapping=aes(fill=f),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+
  #labs(title="CVD1") + 
  custom_theme_sub
print(SubcortCVD_Cent.Sig.f)

# for combined figure 
Tau_CortCVD.f=ggseg(.data=CortCVD,mapping=aes(fill=f),position = "stacked", colour="black",size=.5, hemisphere="left")+
  labs(title="tau") +  custom_theme
Tau_CortCVD.Sig.f=ggseg(.data=CortCVD.Sig,mapping=aes(fill=f),position = "stacked", colour="black",size=.5, hemisphere="left")+custom_theme
Tau_SubCortCVD.f=ggseg(.data=SubcortCVD,mapping=aes(fill=f),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+
  labs(title="tau") +  custom_theme_sub
Tau_SubCortCVD.Sig.f=ggseg(.data=SubcortCVD.Sig,mapping=aes(fill=f),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+custom_theme_sub

rm(CortAge,SubcortAge,CortAge.Sig,SubcortAge.Sig,CortSex,SubcortSex,CortSex.Sig,SubcortSex.Sig,CortAPOE4,SubcortAPOE4,CortAPOE4.Sig,SubcortAPOE4.Sig,
   Cortfollowup_time,Subcortfollowup_time,Cortfollowup_time.Sig,Subcortfollowup_time.Sig,CortCVD,SubcortCVD,CortCVD.Sig,SubcortCVD.Sig,CortCVD_Cent,SubcortCVD_Cent,CortCVD_Cent.Sig,SubcortCVD_Cent.Sig)

# combine plots to get figures
age_row=plot_grid(CortAge.f, SubcortAge.f,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") ) 
sex_row=plot_grid(CortSex.f, SubcortSex.f,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
apoe_row=plot_grid(CortAPOE4.f, SubcortAPOE4.f,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
followup_time_row=plot_grid(Cortfollowup_time.f, Subcortfollowup_time.f,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
cent_row=plot_grid(CortCentiloid.f, SubCortCentiloid.f,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
cvd_row=plot_grid(CortCVD.f, SubcortCVD.f,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
cvd_cent_row=plot_grid(CortCVD_Cent.f, SubcortCVD_Cent.f,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
merged_brain=plot_grid(age_row, sex_row, apoe_row, followup_time_row, cent_row, cvd_row, cvd_cent_row,ncol=1)+theme(plot.background = element_rect(fill = "black"), panel.border = element_blank())
# Save Figure
tiff('./Rabin_Replication_Output/Tau/Colsp_Anova_Unthresholded_1.tiff', units="in", width=7, height=9, res=200,type="cairo" )
print(merged_brain)
# For combined figure 
Colsp_TauUnthresh_cvd_row=plot_grid(Tau_CortCVD.f, align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
Colsp_Tau_row=plot_grid(Tau_CortCVD.f, Tau_CortCVD.Sig.f)  
dev.off()
Colsp_TauUnthresh_cvd_row_SubCort=plot_grid(Tau_SubCortCVD.f, align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
Colsp_Tau_row_SubCort=plot_grid(Tau_SubCortCVD.f, Tau_SubCortCVD.Sig.f)  
dev.off()

age_row=plot_grid(CortAge.Sig.f, SubcortAge.Sig.f,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") ) 
sex_row=plot_grid(CortSex.Sig.f, SubcortSex.Sig.f,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
apoe_row=plot_grid(CortAPOE4.Sig.f, SubcortAPOE4.Sig.f,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
followup_time_row=plot_grid(Cortfollowup_time.Sig.f, Subcortfollowup_time.Sig.f,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
cent_row=plot_grid(CortCentiloid.Sig.f, SubCortCentiloid.Sig.f,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
cvd_row=plot_grid(CortCVD.Sig.f, SubcortCVD.Sig.f,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
cvd_cent_row=plot_grid(CortCVD_Cent.Sig.f, SubcortCVD_Cent.Sig.f,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
merged_brain.sig=plot_grid(age_row, sex_row, apoe_row, followup_time_row, cent_row, cvd_row, cvd_cent_row, ncol=1)+theme(plot.background = element_rect(fill = "black"), panel.border = element_blank())
# Save Figure
tiff('./Rabin_Replication_Output/Tau/Colsp_Anova_Thresholded_1.tiff', units="in", width=7, height=9, res=200,type="cairo")
print(merged_brain.sig)
# For combined figure
Colsp_TauThresh_cvd_row=plot_grid(Tau_CortCVD.Sig.f, align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
dev.off()
Colsp_TauThresh_cvd_row_SubCort=plot_grid(Tau_SubCortCVD.Sig.f, align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
dev.off()
rm(age_row, sex_row, apoe_row, cvd_row, cvd_cent_row, merged_brain,merged_brain.sig,CortAge.f,SubcortAge.f,CortAge.Sig.f,SubcortAge.Sig.f,
   CortSex.f,SubcortSex.f,CortSex.Sig.f,SubcortSex.Sig.f,CortAPOE4.f,SubcortAPOE4.f,CortAPOE4.Sig.f,SubcortAPOE4.Sig.f,
   Cortfollowup_time.f,Subcortfollowup_time.f,Cortfollowup_time.Sig.f,Subcortfollowup_time.Sig.f,
   CortCVD.f,SubcortCVD.f,CortCVD.Sig.f,SubcortCVD.Sig.f,CortCVD_Cent.f,SubcortCVD_Cent.f,CortCVD_Cent.Sig.f,SubcortCVD_Cent.Sig.f)

