#Functions to set up environment
rm(list=ls()) #removes all variables from the current environment
options(scipen = 999) #forces R not to use scientific notation

#Import relevant libraries
library(tibble) 
library(tidyr) 
library(lme4) 
install.packages("lmerTest")
library(lmerTest)
library(ggplot2)

#Load data file
cogn_long<-read.csv("./Data/N151_Long_5-25-22.csv", header=TRUE)
head(cogn_long)

#Change binary variables to factors
cogn_long$Converter_after_tauPET <- as.factor(cogn_long$Converter_after_tauPET)
cogn_long$Biomarker_group_ANY <- as.factor(cogn_long$Biomarker_group_ANY)
cogn_long$Biomarker_group_ANY <- as.factor(cogn_long$Biomarker_group_ANY)
cogn_long$Biomarker_group_ANY <- as.factor(cogn_long$Biomarker_group_ANY)
cogn_long$Biomarker_group_ANY <- as.factor(cogn_long$Biomarker_group_ANY)


#Run lme for full sample, to save individual coefficients for 'decliner' analysis
#(ID is participant identifier. YearsfromBL is cognitive follow up visit in years)
Longcogn_PACC <- lmer (PACC ~  (1 + YrBl|ID), data=cogn_long, contrasts=c("contr.sum", "contr.poly"), REML=FALSE)
summary(Longcogn_PACC)
write.csv(coef(Longcogn_PACC)$ID, './5_25_22/long2_coefficients_PACC_YrBl_ANY_6-14.csv', row.names = TRUE)

#Run lme between pathology groups (regardless of MCI status)
data_group124 <- cogn_long[which(cogn_long$Biomarker_group_ANY!="3"),]

Longcogn_PACC_bygroup_ANY <- lmer (PACC ~ YrBl*Biomarker_group_ANY +(1 + YrBl|ID) + AGE + EDUC + SEX, data=data_group124, contrasts=c("contr.sum", "contr.poly"), REML=FALSE)
summary(Longcogn_PACC_bygroup_ANY)


#Supplementary figure: Separate plot by each metaROI biomarker group, separating by MCI status
#Note for supplementary plots using other tau regions to define tau positivity, just substitute biomarker group variable
#select just those in group 1 (i.e., A+T+)

cogn_group1<- cogn_long[which(cogn_long$Biomarker_group_ANY=="1"),]

#Run lme for A+T+ by MCI status for plotting

Longcogn_PACC_group1_bystatus <- lmer (PACC ~ YrBl*Converter_after_tauPET +(1 + YrBl|ID) + AGE + EDUC + SEX, data=cogn_group1, contrasts=c("contr.sum", "contr.poly"), REML=FALSE)
summary(Longcogn_PACC_group1_bystatus)

#plot PACC
fit_mod <- tibble(fitted(Longcogn_PACC_group1_bystatus))
fit_mod <- cbind(fit_mod, ID = cogn_group1$ID) 
fit_mod <- cbind(fit_mod, BaseCog_Visit = cogn_group1$BaseCog_Visit)
Longcogn_fit = merge(fit_mod, cogn_group1, by=c("ID","BaseCog_Visit"))

Logncogn_PACC_LME_spaghetti_group1_bystatus <- Longcogn_fit %>%
  ggplot(aes(x=YrBl, y=`fitted(Longcogn_PACC_group1_bystatus)`, color=Converter_after_tauPET)) +
  geom_smooth(method=lm, se=TRUE, aes(fill=Converter_after_tauPET), show.legend =FALSE) +
  geom_smooth(method=lm, se=FALSE, aes(fill=Converter_after_tauPET)) +
  geom_rug(aes(x = YrBl), sides="b") +
  scale_color_manual(values = c("skyblue2", "tomato"), labels = c("CU", "MCI")) + scale_fill_manual(values = c("skyblue2", "tomato")) + 
  geom_line(aes(x=YrBl, y=PACC, group = ID, color = Converter_after_tauPET) ,size = 0.50, alpha = .5) +
  labs(x="Visit (Years)", y="PACC (z-score)", title = 'Knight ADRC: A+T+')+
  theme(panel.background = element_blank(), panel.grid.major = element_line(colour = "gray94", size = 0.25), axis.line = element_line(colour = "black"),
        text = element_text(family = "sans", size = 36), plot.title = element_text(hjust = 0.5, size = 48),
        axis.text = element_text(size = 36), legend.text = element_text(size = 36), 
        legend.title=element_blank(), legend.key = element_rect(fill = NA, color = NA),
        plot.margin = unit(c(1, 0, 0, 2), "cm"))+
  guides(fill=FALSE)

print(Logncogn_PACC_LME_spaghetti_group1_bystatus) 
ggsave("./5_25_22/6-14_long_cogn_ADRC_spaghetti_group1_bystatus_ANY.png", plot = Logncogn_PACC_LME_spaghetti_group1_bystatus, device = "png", width = 80, height = 20, units = "cm", dpi = 500)

#select just those in group 2 (i.e., A+T-)

cogn_group2<- cogn_long[which(cogn_long$Biomarker_group_ANY=="2"),]


#Run lme for A+T- by MCI status for plotting

Longcogn_PACC_group2_bystatus <- lmer (PACC ~ YrBl*Converter_after_tauPET +(1 + YrBl|ID) + AGE + EDUC + SEX, data=cogn_group2, contrasts=c("contr.sum", "contr.poly"), REML=FALSE)
summary(Longcogn_PACC_group2_bystatus)

#plot PACC
fit_mod <- tibble(fitted(Longcogn_PACC_group2_bystatus))
fit_mod <- cbind(fit_mod, ID = cogn_group2$ID) 
fit_mod <- cbind(fit_mod, BaseCog_Visit = cogn_group2$BaseCog_Visit)
Longcogn_fit = merge(fit_mod, cogn_group2, by=c("ID","BaseCog_Visit"))


Logncogn_PACC_LME_spaghetti_group2_bystatus <- Longcogn_fit %>%
  ggplot(aes(x=YrBl, y=`fitted(Longcogn_PACC_group2_bystatus)`, color=Converter_after_tauPET)) +
  geom_smooth(method=lm, se=TRUE, aes(fill=Converter_after_tauPET), show.legend =FALSE) +
  geom_smooth(method=lm, se=FALSE, aes(fill=Converter_after_tauPET)) +
  geom_rug(aes(x = YrBl), sides="b") +
  scale_color_manual(values = c("skyblue2", "tomato"), labels = c("CU", "MCI")) + scale_fill_manual(values = c("skyblue2", "tomato")) + 
  geom_line(aes(x=YrBl, y=PACC, group = ID, color = Converter_after_tauPET) ,size = 0.50, alpha = .5) +
  labs(x="Visit (Years)", y="PACC (z-score)", title = 'Knight ADRC: A+T-')+
  theme(panel.background = element_blank(), panel.grid.major = element_line(colour = "gray94", size = 0.25), axis.line = element_line(colour = "black"),
        text = element_text(family = "sans", size = 36), plot.title = element_text(hjust = 0.5, size = 48),
        axis.text = element_text(size = 36), legend.text = element_text(size = 36), 
        legend.title=element_blank(), legend.key = element_rect(fill = NA, color = NA),
        plot.margin = unit(c(1, 0, 0, 2), "cm"))+
  guides(fill=FALSE)

print(Logncogn_PACC_LME_spaghetti_group2_bystatus) 
ggsave("./5_25_22/6-14_long_cogn_ADRC_spaghetti_group2_bystatus_ANY.png", plot = Logncogn_PACC_LME_spaghetti_group2_bystatus, device = "png", width = 80, height = 20, units = "cm", dpi = 500)

##For A-T+, two options: if just one subject, can do just spaghetti plot. 
##if more than one, can do lme plot plus spaghetti

#select just those in group 3 (i.e., A-T+)

cogn_group3<- cogn_long[which(cogn_long$Biomarker_group_ANY=="3"),]


#plot spaghetti plot


#Plot just spaghetti (no LME)

Logncogn_PACC_group3_line <- ggplot(data = cogn_group3) +
  geom_line(aes(x=YrBl, y=PACC, group = ID, colour = "skyblue2"), size = 0.50, alpha = 1.0)+
  labs(x="Visit (Years)", y="PACC (z-score)", title = 'Knight ADRC: A-T+')+
  scale_color_manual(name = "none", values = c("skyblue2"), labels = c("CU"))+
  theme(panel.background = element_blank(), panel.grid.major = element_line(colour = "gray94", size = 0.25), axis.line = element_line(colour = "black"),
        text = element_text(family = "sans", size = 36), plot.title = element_text(hjust = 0.5, size = 48),
        axis.text = element_text(size = 36), legend.text = element_text(size = 36), 
        legend.title=element_blank(), legend.key = element_rect(fill = NA, color = NA),
        plot.margin = unit(c(1, 0, 0, 2), "cm"))+
  ggplot2::scale_x_continuous(limits=c(0, 15), breaks = c(0, 5, 10, 15))+
  guides(fill=FALSE)

print(Logncogn_PACC_group3_line)
ggsave("./5_25_22/6-14_long_cogn_ADRC_group3_line_ANY.png", plot = Logncogn_PACC_group3_line, device = "png", width = 80, height = 20, units = "cm", dpi = 500)

#Run lme for A-T+ by MCI status for plotting

Longcogn_PACC_group3 <- lmer (PACC ~ (1 + YrBl|ID) + AGE + SEX + EDUC, data=cogn_group3, contrasts=c("contr.sum", "contr.poly"), REML=FALSE)
summary(Longcogn_PACC_group3)

#plot PACC
fit_mod <- tibble(fitted(Longcogn_PACC_group3))
fit_mod <- cbind(fit_mod, ID = cogn_group3$ID) 
fit_mod <- cbind(fit_mod, BaseCog_Visit = cogn_group3$BaseCog_Visit)
Longcogn_fit = merge(fit_mod, cogn_group3, by=c("ID","BaseCog_Visit"))

Logncogn_PACC_LME_spaghetti_group3 <- Longcogn_fit %>%
  ggplot(aes(x=YrBl, y=`fitted(Longcogn_PACC_group3)`, color=Converter_after_tauPET)) +
  geom_smooth(method=lm, se=TRUE, aes(fill=Converter_after_tauPET), show.legend =FALSE) +
  geom_smooth(method=lm, se=FALSE, aes(fill=Converter_after_tauPET)) +
  geom_rug(aes(x = YrBl), sides="b") +
  scale_color_manual(values = c("skyblue2", "tomato"), labels = c("CU", "MCI")) + scale_fill_manual(values = c("skyblue2", "tomato")) + 
  geom_line(aes(x=YrBl, y=PACC, group = ID, color = Converter_after_tauPET) ,size = 0.50, alpha = .5) +
  labs(x="Visit (Years)", y="PACC (z-score)", title = 'Knight ADRC: A-T+')+
  theme(panel.background = element_blank(), panel.grid.major = element_line(colour = "gray94", size = 0.25), axis.line = element_line(colour = "black"),
        text = element_text(family = "sans", size = 36), plot.title = element_text(hjust = 0.5, size = 48),
        axis.text = element_text(size = 36), legend.text = element_text(size = 36), 
        legend.title=element_blank(), legend.key = element_rect(fill = NA, color = NA),
        plot.margin = unit(c(1, 0, 0, 2), "cm"))+
  ggplot2::scale_x_continuous(limits=c(0, 13), breaks = c(0, 5, 10))+
  guides(fill=FALSE)

print(Logncogn_PACC_LME_spaghetti_group3) 
ggsave("./5_25_22/6-14_long_cogn_ADRC_spaghetti_group3_ANY.png", plot = Logncogn_PACC_LME_spaghetti_group3, device = "png", width = 80, height = 20, units = "cm", dpi = 500)


#select just those in group 4 (i.e., A-T-)

cogn_group4<- cogn_long[which(cogn_long$Biomarker_group_ANY=="4"),]


#Run lme for A-T- by MCI status for plotting

Longcogn_PACC_group4_bystatus <- lmer (PACC ~ YrBl*Converter_after_tauPET +(1 + YrBl|ID) + AGE + SEX + EDUC, data=cogn_group4, contrasts=c("contr.sum", "contr.poly"), REML=FALSE)
summary(Longcogn_PACC_group4_bystatus)


#plot PACC
fit_mod <- tibble(fitted(Longcogn_PACC_group4_bystatus))
fit_mod <- cbind(fit_mod, ID = cogn_group4$ID) 
fit_mod <- cbind(fit_mod, BaseCog_Visit = cogn_group4$BaseCog_Visit)
Longcogn_fit = merge(fit_mod, cogn_group4, by=c("ID","BaseCog_Visit"))

Logncogn_PACC_LME_spaghetti_group4_bystatus <- Longcogn_fit %>%
  ggplot(aes(x=YrBl, y=`fitted(Longcogn_PACC_group4_bystatus)`, color=Converter_after_tauPET)) +
  geom_smooth(method=lm, se=TRUE, aes(fill=Converter_after_tauPET), show.legend =FALSE) +
  geom_smooth(method=lm, se=FALSE, aes(fill=Converter_after_tauPET)) +
  geom_rug(aes(x = YrBl), sides="b") +
  scale_color_manual(values = c("skyblue2", "tomato"), labels = c("CU", "MCI")) + scale_fill_manual(values = c("skyblue2", "tomato")) + 
  geom_line(aes(x=YrBl, y=PACC, group = ID, color = Converter_after_tauPET) ,size = 0.50, alpha = .5) +
  labs(x="Visit (Years)", y="PACC (z-score)", title = 'Knight ADRC: A-T-')+
  theme(panel.background = element_blank(), panel.grid.major = element_line(colour = "gray94", size = 0.25), axis.line = element_line(colour = "black"),
        text = element_text(family = "sans", size = 36), plot.title = element_text(hjust = 0.5, size = 48),
        axis.text = element_text(size = 36), legend.text = element_text(size = 36), 
        legend.title=element_blank(), legend.key = element_rect(fill = NA, color = NA),
        plot.margin = unit(c(1, 0, 0, 2), "cm"))+
  guides(fill=FALSE)

print(Logncogn_PACC_LME_spaghetti_group4_bystatus) 
ggsave("./5_25_22/6-14_long_cogn_ADRC_spaghetti_group4_bystatus_ANY.png", plot = Logncogn_PACC_LME_spaghetti_group4_bystatus, device = "png", width = 80, height = 20, units = "cm", dpi = 500)

#for combining figures
library(patchwork)


longcogn_ADRC <-
  
  (Logncogn_PACC_LME_spaghetti_group1_bystatus  |
     
     Logncogn_PACC_LME_spaghetti_group2_bystatus  |
     
     Logncogn_PACC_LME_spaghetti_group3 |
     
     Logncogn_PACC_LME_spaghetti_group4_bystatus )

print(longcogn_ADRC)

ggsave("./5_25_22/Long_Biomarker_group_ANY_6-14.png", plot = longcogn_ADRC, device = "png", width = 80, height = 20, units = "cm", dpi = 500)

