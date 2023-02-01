# Title: Cox Regression Analyses for ADRC Cohort for Tau/MCI Replication Study (Sylvia/Cherie)
# Written by Cherie, edited by Diana
# Author: Diana Hobbs
# Date: May 2022

#Functions to set up environment
rm(list=ls()) #removes all variables from the current environment
options(scipen = 999) #forces R not to use scientific notation

#Install & load packages
#install.packages(c("survival", "survminer", "ggfortify"))
library("survival")
library("survminer")
library("ggfortify")
library("ggplot2")
library("RColorBrewer")

#Lload data file
data<-read.csv("./Data/N151_Wide_5-25-22.csv", header=TRUE)

head(data)

#Notes on data: 
#Biomarker_group variable should be in the format: 1: A+T+, 2: A+T-, 3:A-T+, 4: A-T-
#Converter_after_tauPET variable: 1: converter to MCI, 0: didn't convert
#'Cox_followuptime' variable: For noncoverters: total years of follow-up after tauPET; for
#converters: time to MCI conversion

#Ensure variables are in correct format
data$Cox_followuptime <- as.numeric(data$Cox_followuptime)
data$Biomarker_group_metaROI <- as.factor(data$Biomarker_group_metaROI)
data$Biomarker_group_EC <- as.factor(data$Biomarker_group_EC)
data$Biomarker_group_IT <- as.factor(data$Biomarker_group_IT)
data$Biomarker_group_ANY <- as.factor(data$Biomarker_group_ANY)
data<-data%>%
  mutate(SEX = str_replace_all(SEX, c("Male"="0", "Female"="1")),
         SEX = as.factor(SEX))
data$APOE4 <-as.factor(data$APOE4)
#data$Converter_after_tauPET <- as.factor(data$Converter_after_tauPET)

##COX REGRESSION MODELS FOR REPORTING

#First remove A-T+ participant/s
data_group124 <- data[which(data$Biomarker_group_metaROI!="3"),]

#Run main Cox model 
res.cox_multi <- coxph(Surv(Cox_followuptime, Converter_after_tauPET) ~ Biomarker_group_metaROI + AGE + SEX + EDUC + APOE4, data = data_group124)
summary(res.cox_multi)
res.cox_multi_sum<-summary(res.cox_multi); res.cox_multi_cof<-res.cox_multi_sum$coefficients
res.cox_multi_COF<-list(); res.cox_multi_COF<-res.cox_multi_cof
write.csv(res.cox_multi_COF,"./5_25_22/res.cox_multi_COF_N151_NoMMSE.csv")

#res.cox_multi_test<- matrix(0, 3, 3)
#colnames(res.cox_multi_test)<-c("test", "df", "pvalue"); rownames(res.cox_multi_test)<-c("logtest", "sctest", "waldtest")
#res.cox_multi_test[1,]=res.cox_multi_sum$logtest; res.cox_multi_test[2,]=res.cox_multi_sum$sctest; res.cox_multi_test[3,]=res.cox_multi_sum$waldtest
#write.csv(res.cox_multi_test,"./5_25_22/res.cox_multi_test_N151.csv")

#Cox model with just demog/clinical variables
res.cox_demog<- coxph(Surv(Cox_followuptime, Converter_after_tauPET) ~ AGE + SEX + EDUC + APOE4 + MMSE_TauCDR1, data = data_group124)
summary(res.cox_demog)
res.cox_demog_sum<-summary(res.cox_demog); res.cox_demog_cof<-res.cox_demog_sum$coefficients
res.cox_demog_COF<-list(); res.cox_demog_COF<-res.cox_demog_cof
write.csv(res.cox_demog_COF,"./5_25_22/res.cox_demog_COF_N151_YesMMSE.csv")

#res.cox_demog_test<- matrix(0, 3, 3)
#colnames(res.cox_demog_test)<-c("test", "df", "pvalue"); rownames(res.cox_demog_test)<-c("logtest", "sctest", "waldtest")
#res.cox_demog_test[1,]=res.cox_demog_sum$logtest; res.cox_demog_test[2,]=res.cox_demog_sum$sctest; res.cox_demog_test[3,]=res.cox_demog_sum$waldtest
#write.csv(res.cox_demog_test,"./5_25_22/res.cox_demog_test_N151.csv")

#Main Cox model plus continuous measures of N
res.cox_multi_vol <- coxph(Surv(Cox_followuptime, Converter_after_tauPET) ~ Biomarker_group_metaROI + AGE + SEX + EDUC + APOE4 + PercentHIPP, data = data_group124)
summary(res.cox_multi_vol)
res.cox_multi_vol_sum<-summary(res.cox_multi_vol); res.cox_multi_vol_cof<-res.cox_multi_vol_sum$coefficients
res.cox_multi_vol_COF<-list(); res.cox_multi_vol_COF<-res.cox_multi_vol_cof
write.csv(res.cox_multi_vol_COF,"./5_25_22/res.cox_multi_vol_COF_N151_NoMMSE.csv")

res.cox_multi_thickness <- coxph(Surv(Cox_followuptime, Converter_after_tauPET) ~ Biomarker_group_metaROI + AGE + SEX + EDUC + APOE4 + Cort_Thick, data = data_group124)
summary(res.cox_multi_thickness)
res.cox_multi_thickness_sum<-summary(res.cox_multi_thickness); res.cox_multi_thickness_cof<-res.cox_multi_thickness_sum$coefficients
res.cox_multi_thickness_COF<-list(); res.cox_multi_thickness_COF<-res.cox_multi_thickness_cof
write.csv(res.cox_multi_thickness_COF,"./5_25_22/res.cox_multi_thickness_COF_N151_NoMMSE.csv")
res.cox_multi_thickness_sum$concordance
#res.cox_multi_vol_test<- matrix(0, 3, 3)
#colnames(res.cox_multi_vol_test)<-c("test", "df", "pvalue"); rownames(res.cox_multi_vol_test)<-c("logtest", "sctest", "waldtest")
#res.cox_multi_vol_test[1,]=res.cox_multi_vol_sum$logtest; res.cox_multi_vol_test[2,]=res.cox_multi_vol_sum$sctest; res.cox_multi_vol_test[3,]=res.cox_multi_vol_sum$waldtest
#write.csv(res.cox_multi_vol_test,"./5_25_22/res.cox_multi_vol_test_N151.csv")

#res.cox_multi_thickness <- coxph(Surv(Cox_followuptime, Converter_after_tauPET) ~ Biomarker_group_metaROI + AGE + SEX + EDUC + APOE4 + Cort_Thick, data = data_group124)
#summary(res.cox_multi_thickness)
#res.cox_multi_thickness_sum<-summary(res.cox_multi_thickness); res.cox_multi_thickness_cof<-res.cox_multi_thickness_sum$coefficients
#res.cox_multi_thickness_COF<-list(); res.cox_multi_thickness_COF<-res.cox_multi_thickness_cof
#write.csv(res.cox_multi_thickness_COF,"./5_25_22/res.cox_multi_thickness_COF_N151.csv")

#res.cox_multi_thickness_test<- matrix(0, 3, 3)
#colnames(res.cox_multi_thickness_test)<-c("test", "df", "pvalue"); rownames(res.cox_multi_thickness_test)<-c("logtest", "sctest", "waldtest")
#res.cox_multi_thickness_test[1,]=res.cox_multi_thickness_sum$logtest; res.cox_multi_thickness_test[2,]=res.cox_multi_thickness_sum$sctest; res.cox_multi_thickness_test[3,]=res.cox_multi_thickness_sum$waldtest
#write.csv(res.cox_multi_thickness_test,"./5_25_22/res.cox_multi_thickness_test_N151.csv")

#Cox model with different tau regions
data_group124 <- data[which(data$Biomarker_group_EC!="3"),]

res.cox_ec <- coxph(Surv(Cox_followuptime, Converter_after_tauPET) ~ Biomarker_group_EC+ AGE + SEX + EDUC + APOE4, data = data_group124)
summary(res.cox_ec)
res.cox_ec_sum<-summary(res.cox_ec); res.cox_ec_cof<-res.cox_ec_sum$coefficients
res.cox_ec_COF<-list(); res.cox_ec_COF<-res.cox_ec_cof
write.csv(res.cox_ec_COF,"./5_25_22/res.cox_ec_COF_N151_NoMMSE.csv")

#res.cox_ec_test<- matrix(0, 3, 3)
#colnames(res.cox_ec_test)<-c("test", "df", "pvalue"); rownames(res.cox_ec_test)<-c("logtest", "sctest", "waldtest")
#res.cox_ec_test[1,]=res.cox_ec_sum$logtest; res.cox_ec_test[2,]=res.cox_ec_sum$sctest; res.cox_ec_test[3,]=res.cox_ec_sum$waldtest
#write.csv(res.cox_ec_test,"./5_25_22/res.cox_ec_test_N151.csv")

data_group124 <- data[which(data$Biomarker_group_IT!="3"),]

res.cox_it <- coxph(Surv(Cox_followuptime, Converter_after_tauPET) ~ Biomarker_group_IT + AGE + SEX + EDUC + APOE4, data = data_group124)
summary(res.cox_it)
res.cox_it_sum<-summary(res.cox_it); res.cox_it_cof<-res.cox_it_sum$coefficients
res.cox_it_COF<-list(); res.cox_it_COF<-res.cox_it_cof
write.csv(res.cox_it_COF,"./5_25_22/res.cox_it_COF_N151_NoMMSE.csv")

#res.cox_it_test<- matrix(0, 3, 3)
#colnames(res.cox_it_test)<-c("test", "df", "pvalue"); rownames(res.cox_it_test)<-c("logtest", "sctest", "waldtest")
#res.cox_it_test[1,]=res.cox_it_sum$logtest; res.cox_it_test[2,]=res.cox_it_sum$sctest; res.cox_it_test[3,]=res.cox_it_sum$waldtest
#write.csv(res.cox_it_test,"./5_25_22/res.cox_it_test_N151.csv")

data_group124 <- data[which(data$Biomarker_group_ANY!="3"),]

res.cox_any <- coxph(Surv(Cox_followuptime, Converter_after_tauPET) ~ Biomarker_group_ANY + AGE + SEX + EDUC + APOE4, data = data_group124)
summary(res.cox_any)
res.cox_any_sum<-summary(res.cox_any); res.cox_any_cof<-res.cox_any_sum$coefficients
res.cox_any_COF<-list(); res.cox_any_COF<-res.cox_any_cof
write.csv(res.cox_any_COF,"./5_25_22/res.cox_any_COF_N151_NoMMSE.csv")

#res.cox_any_test<- matrix(0, 3, 3)
#colnames(res.cox_any_test)<-c("test", "df", "pvalue"); rownames(res.cox_any_test)<-c("logtest", "sctest", "waldtest")
#res.cox_any_test[1,]=res.cox_any_sum$logtest; res.cox_any_test[2,]=res.cox_any_sum$sctest; res.cox_any_test[3,]=res.cox_any_sum$waldtest
#write.csv(res.cox_any_test,"./5_25_22/res.cox_any_test_N151.csv")


#CSV files for log test and concordance
res.cox_concord<- matrix(0, 7, 5)
colnames(res.cox_concord)<-c("C", "se(C)", "likelihood ratio", "df", "p")
rownames(res.cox_concord)<-c("multi","demog","multi_vol", "cortthick", "EC(124)", "IT(124)", "ANY(124)")
res.cox_concord[1,1:2]=res.cox_multi_sum$concordance
res.cox_concord[2,1:2]=res.cox_demog_sum$concordance
res.cox_concord[3,1:2]=res.cox_multi_vol_sum$concordance
res.cox_concord[4,1:2]=res.cox_multi_thickness_sum$concordance
res.cox_concord[5,1:2]=res.cox_ec_sum$concordance
res.cox_concord[6,1:2]=res.cox_it_sum$concordance
res.cox_concord[7,1:2]=res.cox_any_sum$concordance
#Must open this file to put in likelihood ratio tests 
write.csv(res.cox_concord,"./5_25_22/res.cox_concord_N151_NoMMSEexceptDemog.csv")

###In plots, use all data for visualisation purposes (i.e., do not exclude A-T+)
#PLOTTING WITH TEMPORAL METAROI

#Plot survival curves by biomarker group
fit <- survfit(Surv(Cox_followuptime, Converter_after_tauPET) ~ Biomarker_group_metaROI, data = data)
plotadrc <- autoplot(fit, conf.int = FALSE, censor = FALSE, surv.size = 2) +
  ggplot2::theme(panel.background = element_blank(), panel.grid.major = element_line(colour = "gray94", size = 0.25), axis.line = element_line(colour = "black"),
                                                                                            text = element_text(family = "sans", size = 24), plot.title = element_text(hjust = 0.5, size = 28),axis.text = element_text(size = 24), 
                                                                                            legend.position = "none", legend.text = element_text(size = 24), legend.key = element_rect(fill = NA, colour = NA)) +
  ggplot2::labs(x="Time (Years)", y="Survival without MCI", title = 'Knight ADRC', 
                colour = "Biomarker group") + ggplot2::scale_color_manual(values = c("#FC8D62", "#E78AC3", "#66C2A5", "#FFD92F"))+#, labels = c("A+T+", "A+T-", "A-T+", "A-T-"))+
  ggplot2::scale_x_continuous(limits=c(0, 6), breaks = c(0, 2, 4, 6))+
  ggplot2::scale_y_continuous(limits=c(.58, 1), breaks = c(.60, .80, 1), labels=c("60%", "80%", "100%"))


print(plotadrc)

#Save as image
ggsave("./5_25_22/20220525fig_adrc_survival_N151_NoMMSE.png", plot = plotadrc, device = "png", width = 20.41, height = 15.72, units = "cm", dpi = 500)


#PLOTTING WITH ENTORHINAL CORTEX

#Plot survival curves by biomarker group
fit <- survfit(Surv(Cox_followuptime, Converter_after_tauPET) ~ Biomarker_group_EC, data = data)
plotadrc_ec <- autoplot(fit, conf.int = FALSE, censor = FALSE, surv.size = 2) + 
  ggplot2::theme(panel.background = element_blank(), 
                 panel.grid.major = element_line(colour = "gray94", size = 0.25), 
                 axis.line = element_line(colour = "black"), text = element_text(family = "sans", size = 24), plot.title = element_text(hjust = 0.5, size = 28), 
                 axis.text = element_text(size = 24), legend.position = "none", legend.text = element_text(size = 24), legend.key = element_rect(fill = NA, colour = NA))+
  ggplot2::labs(x="Time (Years)", y="Survival without MCI", title = 'Knight ADRC', 
                colour = "Biomarker group (EC)") + ggplot2::scale_color_manual(values = c("#FC8D62", "#E78AC3", "#66C2A5", "#FFD92F"), labels = c("A+T+", "A+T-", "A-T+", "A-T-"))+
  ggplot2::scale_x_continuous(limits=c(0, 6), breaks = c(0, 2, 4, 6))+
  ggplot2::scale_y_continuous(limits=c(.51, 1), breaks = c(.60, .80, 1), labels=c("60%", "80%", "100%"))


print(plotadrc_ec)

#PLOTTING WITH INFERIOR TEMPORAL

#Plot survival curves by biomarker group
fit <- survfit(Surv(Cox_followuptime, Converter_after_tauPET) ~ Biomarker_group_IT, data = data)
plotadrc_it <- autoplot(fit, conf.int = FALSE, censor = FALSE, surv.size = 2) +
  ggplot2::theme(panel.background = element_blank(), panel.grid.major = element_line(colour = "gray94", size = 0.25), axis.line = element_line(colour = "black"),
                                       text = element_text(family = "sans", size = 24), plot.title = element_text(hjust = 0.5, size = 28),
                                       axis.text = element_text(size = 24), legend.position = "none", legend.text = element_text(size = 24), legend.key = element_rect(fill = NA, colour = NA)) +
  ggplot2::labs(x="Time (Years)", y="Survival without MCI", title = 'Knight ADRC', 
                colour = "Biomarker group (IT)") + ggplot2::scale_color_manual(values = c("#FC8D62", "#E78AC3", "#66C2A5", "#FFD92F"), labels = c("A+T+", "A+T-", "A-T+", "A-T-"))+
  ggplot2::scale_x_continuous(limits=c(0, 6), breaks = c(0, 2, 4, 6))+
  ggplot2::scale_y_continuous(limits=c(.50, 1), breaks = c(.60, .80, 1), labels=c("60%", "80%", "100%"))



print(plotadrc_it)

#PLOTTING WITH ANY

#Plot survival curves by biomarker group
fit <- survfit(Surv(Cox_followuptime, Converter_after_tauPET) ~ Biomarker_group_ANY, data = data)
plotadrc_any <- autoplot(fit, conf.int = FALSE, censor = FALSE, surv.size = 2)+
  ggplot2::theme(panel.background = element_blank(), panel.grid.major = element_line(colour = "gray94", size = 0.25), axis.line = element_line(colour = "black"),
                                       text = element_text(family = "sans", size = 24), plot.title = element_text(hjust = 0.5, size = 28),
                                       axis.text = element_text(size = 24), legend.position = "none", legend.text = element_text(size = 24), legend.key = element_rect(fill = NA, colour = NA)) +
  ggplot2::labs(x="Time (Years)", y="Survival without MCI", title = 'Knight ADRC', 
                colour = "Biomarker group (ANY)") + ggplot2::scale_color_manual(values = c("#FC8D62", "#E78AC3", "#66C2A5", "#FFD92F"), labels = c("A+T+", "A+T-", "A-T+", "A-T-"))+
  ggplot2::scale_x_continuous(limits=c(0, 6), breaks = c(0, 2, 4, 6))+
  ggplot2::scale_y_continuous(limits=c(.51, 1), breaks = c(.60, .80, 1), labels=c("60%", "80%", "100%"))
  

print(plotadrc_any)


#for combining figures
library(patchwork)

#Plot for supplementary

combined_supp_plot <- plotadrc_ec / plotadrc_it / plotadrc_any

print(combined_supp_plot)

ggsave("./5_25_22/20220525_suppfig_survival_padandADRC_N151_NoMMSE.png", plot = combined_supp_plot, device = "png", width = 20.41, height = 45.72, units = "cm", dpi = 500)


data <- data%>%
  mutate(followuptime_Groups = case_when(
    Converter_after_tauPET==0 & Cox_followuptime <= 2 ~ 2,
    Converter_after_tauPET==0 & Cox_followuptime > 2 & Cox_followuptime <= 4 ~ 4,
    Converter_after_tauPET==0 & Cox_followuptime > 4 & Cox_followuptime <= 6 ~ 6,
    Converter_after_tauPET==0 & Cox_followuptime > 6 ~ 8,
  ))

meta1<- data[which(data$Biomarker_group_metaROI=="1"),]
meta2<- data[which(data$Biomarker_group_metaROI=="2"),]
meta3<- data[which(data$Biomarker_group_metaROI=="3"),]
meta4<- data[which(data$Biomarker_group_metaROI=="4"),]
table(meta1$followuptime_Groups)
table(meta2$followuptime_Groups)
table(meta3$followuptime_Groups)
table(meta4$followuptime_Groups)


EC1<- data[which(data$Biomarker_group_EC=="1"),]
EC2<- data[which(data$Biomarker_group_EC=="2"),]
EC3<- data[which(data$Biomarker_group_EC=="3"),]
EC4<- data[which(data$Biomarker_group_EC=="4"),]
table(EC1$followuptime_Groups)
table(EC2$followuptime_Groups)
table(EC3$followuptime_Groups)
table(EC4$followuptime_Groups)

IT1<- data[which(data$Biomarker_group_IT=="1"),]
IT2<- data[which(data$Biomarker_group_IT=="2"),]
IT3<- data[which(data$Biomarker_group_IT=="3"),]
IT4<- data[which(data$Biomarker_group_IT=="4"),]
table(IT1$followuptime_Groups)
table(IT2$followuptime_Groups)
table(IT3$followuptime_Groups)
table(IT4$followuptime_Groups)

ANY1<- data[which(data$Biomarker_group_ANY=="1"),]
ANY2<- data[which(data$Biomarker_group_ANY=="2"),]
ANY3<- data[which(data$Biomarker_group_ANY=="3"),]
ANY4<- data[which(data$Biomarker_group_ANY=="4"),]
table(ANY1$followuptime_Groups)
table(ANY2$followuptime_Groups)
table(ANY3$followuptime_Groups)
table(ANY4$followuptime_Groups)

data <- data%>%
  mutate(followuptime_Groups = case_when(
    Converter_after_tauPET==1 & Cox_followuptime <= 2 ~ 2,
    Converter_after_tauPET==1 & Cox_followuptime > 2 & Cox_followuptime <= 4 ~ 4,
    Converter_after_tauPET==1 & Cox_followuptime > 4 & Cox_followuptime <= 6 ~ 6,
    Converter_after_tauPET==1 & Cox_followuptime > 6 ~ 8,
  ))

meta1<- data[which(data$Biomarker_group_metaROI=="1"),]
meta2<- data[which(data$Biomarker_group_metaROI=="2"),]
meta3<- data[which(data$Biomarker_group_metaROI=="3"),]
meta4<- data[which(data$Biomarker_group_metaROI=="4"),]
table(meta1$followuptime_Groups)
table(meta2$followuptime_Groups)
table(meta3$followuptime_Groups)
table(meta4$followuptime_Groups)


EC1<- data[which(data$Biomarker_group_EC=="1"),]
EC2<- data[which(data$Biomarker_group_EC=="2"),]
EC3<- data[which(data$Biomarker_group_EC=="3"),]
EC4<- data[which(data$Biomarker_group_EC=="4"),]
table(EC1$followuptime_Groups)
table(EC2$followuptime_Groups)
table(EC3$followuptime_Groups)
table(EC4$followuptime_Groups)

IT1<- data[which(data$Biomarker_group_IT=="1"),]
IT2<- data[which(data$Biomarker_group_IT=="2"),]
IT3<- data[which(data$Biomarker_group_IT=="3"),]
IT4<- data[which(data$Biomarker_group_IT=="4"),]
table(IT1$followuptime_Groups)
table(IT2$followuptime_Groups)
table(IT3$followuptime_Groups)
table(IT4$followuptime_Groups)

ANY1<- data[which(data$Biomarker_group_ANY=="1"),]
ANY2<- data[which(data$Biomarker_group_ANY=="2"),]
ANY3<- data[which(data$Biomarker_group_ANY=="3"),]
ANY4<- data[which(data$Biomarker_group_ANY=="4"),]
table(ANY1$followuptime_Groups)
table(ANY2$followuptime_Groups)
table(ANY3$followuptime_Groups)
table(ANY4$followuptime_Groups)

