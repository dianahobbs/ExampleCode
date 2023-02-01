# Title: ADRC - PiB Regression Loop 
# Author: Diana Hobbs
# Date: April 2022

#### Set Up ####
# Load packages 
pacman::p_load(tidyverse,jtools,naniar,lubridate,lme4,sfsmisc,foreign,MASS,scales,ggseg,cowplot)

###### LOAD & SETUP DATA 
PiB <- read.csv("./Data/Rabin_Replication_PIB.csv")%>%
  mutate(CVD = str_replace_all(CVD_Colsp, c("<5%"="0", "5-10%"="1", "10-20%"="2",">20%"="3")))

# Reassign variable types & set reference 
PiB$SEX <- as.factor(PiB$SEX); PiB$SEX <- relevel(PiB$SEX, "Male") 
PiB$APOE4 <- as.factor(PiB$APOE4); PiB$APOE4 <- relevel(PiB$APOE4, "0")
PiB$CVD <- as.factor(PiB$CVD); PiB$CVD <- relevel(PiB$CVD, "0")

###### REGRESSION LOOP
ROI <- PiB[,20:60]

# Loop through regressions and keep output 
MODEL <- list() # output for ANOVA 

# Create a 'for' loop 
for (i in 1:ncol(ROI)){
  model <- anova(lm(ROI[,i] ~ AGE + SEX + APOE4 + CVD, data=PiB))
  MODEL[[i]] <- model
  names(MODEL)[i] <- colnames(ROI)[i]
}
write.csv(MODEL,"./Rabin_Replication_Output/PiB/Colsp_Anova_1.csv")
rm(model,MODEL,i,PiB,ROI)

###### Organize loop output in usable format 
Data <- read.csv("./Rabin_Replication_Output/PiB/Colsp_Anova_1.csv")

# Extract F and p-values 
Fvalue <- Data[ ,grepl("F.value", names(Data))]
names(Fvalue) = gsub(pattern = ".F.value", replacement = "", x = names(Fvalue))
Fvalue <- t(Fvalue)
colnames(Fvalue) = c("AGE","SEX","APOE4","CVD", "Residuals")
Fvalue <- as.data.frame(Fvalue)

Pvalue <- Data[ ,grepl("Pr..F.", names(Data))]
names(Pvalue) = gsub(pattern = "Pr..F.", replacement = "", x = names(Pvalue))
Pvalue <- t(Pvalue)
colnames(Pvalue) = c("AGE","SEX","APOE4","CVD","Residuals")
Pvalue <- as.data.frame(Pvalue)

# Merge coefficient summaries & Perform FDR correction using Benjamini and Hochberg method
Age <- data.frame("f"=Fvalue$AGE,"p"=Pvalue$AGE,"p_corr"=(p.adjust(Pvalue$AGE, method = "BH",n = length(Pvalue$AGE))))
Sex <- data.frame("f"=Fvalue$SEX,"p"=Pvalue$SEX,"p_corr"=(p.adjust(Pvalue$SEX, method = "BH",n = length(Pvalue$SEX))))
APOE4 <- data.frame("f"=Fvalue$APOE4,"p"=Pvalue$APOE4,"p_corr"=(p.adjust(Pvalue$APOE4, method = "BH",n = length(Pvalue$APOE4))))
CVD <- data.frame("f"=Fvalue$CVD,"p"=Pvalue$CVD,"p_corr"=(p.adjust(Pvalue$CVD, method = "BH",n = length(Pvalue$CVD))))
rm(Data,Fvalue,Pvalue)

cort_list <- read.csv("./BrianFiles/ggseg_namelist.csv")
subcort_list <- cort_list$subcortical[1:7] #extract subcortical list

CortAge <- data.frame(region=cort_list$cortical,filter(Age[1:34,])); SubcortAge <- data.frame(region=subcort_list,filter(Age[35:41,]))
CortSex <- data.frame(region=cort_list$cortical,filter(Sex[1:34,])); SubcortSex <- data.frame(region=subcort_list,filter(Sex[35:41,]))
CortAPOE4 <- data.frame(region=cort_list$cortical,filter(APOE4[1:34,])); SubcortAPOE4 <- data.frame(region=subcort_list,filter(APOE4[35:41,]))
CortCVD <- data.frame(region=cort_list$cortical,filter(CVD[1:34,])); SubcortCVD <- data.frame(region=subcort_list,filter(CVD[35:41,]))
rm(Age,Sex,APOE4,CVD)

CortAge.Sig <- subset(CortAge, p_corr < 0.05); SubcortAge.Sig <- subset(SubcortAge, p_corr < 0.05)
CortSex.Sig <- subset(CortSex, p_corr < 0.05); SubcortSex.Sig <- subset(SubcortSex, p_corr < 0.05)
CortAPOE4.Sig <- subset(CortAPOE4, p_corr < 0.05); SubcortAPOE4.Sig <- subset(SubcortAPOE4, p_corr < 0.05)
CortCVD.Sig <- subset(CortCVD, p_corr < 0.05); SubcortCVD.Sig <- subset(SubcortCVD, p_corr < 0.05)

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

CortCVD.f=ggseg(.data=CortCVD,mapping=aes(fill=f),position = "stacked", colour="black",size=.5, hemisphere="left")+
  labs(title="CVD") +  custom_theme
print(CortCVD.f)
SubcortCVD.f=ggseg(.data=SubcortCVD,mapping=aes(fill=f),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+
  #labs(title="CVD1") + 
  custom_theme_sub
print(SubcortCVD.f)


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

CortCVD.Sig.f=ggseg(.data=CortCVD.Sig,mapping=aes(fill=f),position = "stacked", colour="black",size=.5, hemisphere="left")+
  labs(title="CVD") +  custom_theme
print(CortCVD.Sig.f)
SubcortCVD.Sig.f=ggseg(.data=SubcortCVD.Sig,mapping=aes(fill=f),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+
  #labs(title="CVD1") + 
  custom_theme_sub
print(SubcortCVD.Sig.f)

# for combined figure 
PiB_CortCVD.f=ggseg(.data=CortCVD,mapping=aes(fill=f),position = "stacked", colour="black",size=.5, hemisphere="left")+
  labs(title="PiB") +  custom_theme
PiB_CortCVD.Sig.f=ggseg(.data=CortCVD.Sig,mapping=aes(fill=f),position = "stacked", colour="black",size=.5, hemisphere="left")+custom_theme
PiB_SubCortCVD.f=ggseg(.data=SubcortCVD,mapping=aes(fill=f),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+
  labs(title="PiB") +  custom_theme_sub
PiB_SubCortCVD.Sig.f=ggseg(.data=SubcortCVD.Sig,mapping=aes(fill=f),position = "stacked", colour="black",size=.5, hemisphere="left", atlas=aseg)+custom_theme_sub

rm(CortAge,SubcortAge,CortAge.Sig,SubcortAge.Sig,CortSex,SubcortSex,CortSex.Sig,SubcortSex.Sig,CortAPOE4,SubcortAPOE4,CortAPOE4.Sig,SubcortAPOE4.Sig,
   CortCVD,SubcortCVD,CortCVD.Sig,SubcortCVD.Sig)

# combine plots to get figures
age_row=plot_grid(CortAge.f, SubcortAge.f,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") ) 
sex_row=plot_grid(CortSex.f, SubcortSex.f,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
apoe_row=plot_grid(CortAPOE4.f, SubcortAPOE4.f,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
cvd_row=plot_grid(CortCVD.f, SubcortCVD.f,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
merged_brain=plot_grid(age_row, sex_row, apoe_row,cvd_row,ncol=1)+theme(plot.background = element_rect(fill = "black"), panel.border = element_blank())
# Save Figure
tiff('./Rabin_Replication_Output/PiB/Colsp_Anova_Unthresholded_1.tiff', units="in", width=7, height=9, res=200,type="cairo" )
print(merged_brain)
# For combined figure 
Colsp_PiBUnthresh_cvd_row=plot_grid(PiB_CortCVD.f, align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
Colsp_PiB_row=plot_grid(PiB_CortCVD.f, PiB_CortCVD.Sig.f)  
dev.off()
Colsp_PiBUnthresh_cvd_row_SubCort=plot_grid(PiB_SubCortCVD.f, align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
Colsp_PiB_row_SubCort=plot_grid(PiB_SubCortCVD.f, PiB_SubCortCVD.Sig.f)  
dev.off()

age_row=plot_grid(CortAge.Sig.f, SubcortAge.Sig.f,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") ) 
sex_row=plot_grid(CortSex.Sig.f, SubcortSex.Sig.f,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
apoe_row=plot_grid(CortAPOE4.Sig.f, SubcortAPOE4.Sig.f,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
cvd_row=plot_grid(CortCVD.Sig.f, SubcortCVD.Sig.f,  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
merged_brain.sig=plot_grid(age_row, sex_row, apoe_row, cvd_row, ncol=1)+theme(plot.background = element_rect(fill = "black"), panel.border = element_blank())
# Save Figure
tiff('./Rabin_Replication_Output/PiB/Colsp_Anova_Thresholded_1.tiff', units="in", width=7, height=9, res=200,type="cairo")
print(merged_brain.sig)
# For combined figure
Colsp_PiBThresh_cvd_row=plot_grid(PiB_CortCVD.Sig.f, align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
dev.off()
Colsp_PiBThresh_cvd_row_SubCort=plot_grid(PiB_SubCortCVD.Sig.f, align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
dev.off()
rm(age_row, sex_row, apoe_row, cvd_row, merged_brain,merged_brain.sig,CortAge.f,SubcortAge.f,CortAge.Sig.f,SubcortAge.Sig.f,
   CortSex.f,SubcortSex.f,CortSex.Sig.f,SubcortSex.Sig.f,CortAPOE4.f,SubcortAPOE4.f,CortAPOE4.Sig.f,SubcortAPOE4.Sig.f,
   CortCVD.f,SubcortCVD.f,CortCVD.Sig.f,SubcortCVD.Sig.f)

