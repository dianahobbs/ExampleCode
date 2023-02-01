#### Cefepime data project 
### Setting up the data and running analyses
## Author: Diana Hobbs
## Date: May 2022

#### Set Up Packages and Load Data ####
pacman::p_load(tidyverse,jtools,psych,cowplot,multcomp,gmodels,lubridate)

data <- read.csv("./OriginalData.csv", header = TRUE)

data <- data%>%
  rename(Age=Age..................................years., 
         Sex=Sex..............Male...1.,
         CefeExtendedInfusion=Change.to.Cefepime.Extended.Infusion...............Yes...1.,
         CefeExtendedInfusion_Dose=Change.to.Cefepime.Extended.Infusion.Dose...............n...g.,
         ActBW_kg=Actual.Body.Weight...............................ABW...kg., 
         BaseCreatineClear=Baseline.Creatine.Clearance........................eCrCl.mL.min.,
         Tmax1=Tmax.Day.1................n..F.,
         Tmax2=Tmax.Day.2................n..F.,
         Tmax3=Tmax.Day.3................n..F.,
         Tmax4=Tmax.Day.4................n..F.,
         Tmax5=Tmax.Day.5................n..F.,
         Tmax6=Tmax.Day.6................n..F.,
         Tmax7=Tmax.Day.7................n..F.,
         TimetoDefervescence=Time.to.defervescence.........n...hours..only.if..24h.,
         GCSFAdminDuration=Duration.of.GCSF.administration............n...days.,
         DayssinceTherapy=If.yes.to.previous..days.since.last.therapy.........................n....days.,
         LengthofStay_days=Length.of.Stay.................n...days.,
         BaseSCrIncrease_0.3=Baseline.SCr.increased.0.3.mg.dL...............Yes.1.,
         BaseSCrIncrease_1.5=Baseline.SCr.increase.1.5x............Yes.1.,
         AKI=TOTAL.WITH.AKI,
         BaseProcalcitonin=Electrolyte.abnormalities..Baseline.Procalcitonin............1...yes., 
         #BaseProcalcitonin_level=Electrolyte.abnormalities..Baseline.Procalcitonin,
         BaseLDH=Electrolyte.abnormalities..Baseline.LDH...............1...yes.,
         #BaseLDH_level=Electrolyte.abnormalities..Baseline.LDH,
         BaseLA=Electrolyte.abnormalities..Baseline.LA...............1...yes.,
         #BaseLA_level=Electrolyte.abnormalities..Baseline.LA,
         BaseK=Electrolyte.abnormalities..Baseline.Potassium.............1...yes.,
         AdminGCSF=Patient.administered...givent.GCSF.......................Yes...1.,
         CultureCollect=Microbiology.Culture.Collected........................Yes...1.,
         Blood=Growth.from.Blood.Culture, 
         Sputum=Growth.from.Sputum.Culture,
         Urine=Growth.from.Urine.Culture,
         Misc=Miscellaneous.Culture.......................ie.intrabdominal..abcess.,
         BUG=Microbiology.Culture.growth..............BUG.,
         Mortality_30day=X30.day.mortality......Yes...1.,
         AntimicroEsclation=Antimicrobial.escalation....................Yes..1.,
         SIRS_CriteriaMet=SIRS.Criteria.Met..................Yes...1.,
         TherapyDuration_days=Duration.of.Therapy...................n...days.,
         ConfirmedInfection=Confirmed.or.suspected.site.of.infection................Yes...1.,
         RecordedFever=Fever.recorded...documented.as.home.reading....................n..F.,
         ReceivedCefeLD=Received.Cefepime.LD...Y.1.,
         DurationNeutropenia=Neutropenia.duration...days.)%>%
  dplyr::select(Dates,Age,Sex,CefeExtendedInfusion,CefeExtendedInfusion_Dose,ActBW_kg,BaseCreatineClear,
         Tmax1,Tmax2,Tmax3,Tmax4,Tmax5,Tmax6,Tmax7,TimetoDefervescence,GCSFAdminDuration,DayssinceTherapy,
         LengthofStay_days,BaseSCrIncrease_0.3,BaseSCrIncrease_1.5,AKI,BaseProcalcitonin,BaseLDH,BaseLA,BaseK,
         AdminGCSF,CultureCollect,Blood,Sputum,Urine,Misc,BUG,Mortality_30day,AntimicroEsclation,SIRS_CriteriaMet,
         TherapyDuration_days,ConfirmedInfection,RecordedFever,ReceivedCefeLD,DurationNeutropenia)
        

#MicrobioCount=!!!!!
         
# Load packages for survival curve 
pacman::p_load(survival,survminer,GGally,ggplot2)    

data$Dates <- as.Date(data$Dates, "%m/%d/%Y")
data$Date_Yr <- format(as.Date(data$Dates, format="%d/%m/%Y"),"%Y")

# using time period from april of 2016
data <- data%>%
  mutate(Accurate_TimePeriod = case_when(Dates <= "2016-04-01" ~ 0,
                                         Dates > "2016-04-01" ~ 1))
data <- data%>%
  mutate(TimePeriod = case_when(Date_Yr <= 2016 ~ 0,
                                Date_Yr > 2016 ~ 1))
 
Defervescence <- data%>%dplyr::select(CefeExtendedInfusion,BaseLDH,BaseProcalcitonin,SIRS_CriteriaMet,AdminGCSF,ReceivedCefeLD,TimetoDefervescence,TimePeriod,Accurate_TimePeriod)%>%
           mutate(BrokeFever_24=case_when(TimetoDefervescence <= 24 ~ 1,
                                          TimetoDefervescence > 24 ~ 0),
                  BrokeFever_48=case_when(TimetoDefervescence > 24 & TimetoDefervescence <=48 ~ 1,
                                          TimetoDefervescence > 48 ~ 0),
                  BrokeFever_72=case_when(TimetoDefervescence > 48 & TimetoDefervescence <=72 ~ 1,
                                          TimetoDefervescence > 72 ~ 0),
                  BrokeFever_96=case_when(TimetoDefervescence > 72 & TimetoDefervescence <=96 ~ 1,
                                          TimetoDefervescence > 96 ~ 0),
                  BrokeFever_120=case_when(TimetoDefervescence > 96 & TimetoDefervescence <=120 ~ 1,
                                          TimetoDefervescence > 120 ~ 0),
                  BrokeFever_144=case_when(TimetoDefervescence > 120 & TimetoDefervescence <=144 ~ 1,
                                          TimetoDefervescence > 144 ~ 0),
                  BrokeFever_168=case_when(TimetoDefervescence > 144 & TimetoDefervescence <=168 ~ 1,
                                          TimetoDefervescence > 168 ~ 0))

Defervescence[is.na(Defervescence)] <- 1

# Descriptive Stats and Chi Square for Groups
Defervescence%>%group_by(CefeExtendedInfusion)%>% # Get the median time to defervescence per Tx group
  summarize(Median_TtoD=median(TimetoDefervescence))

Tx_0 <- Defervescence%>%filter(CefeExtendedInfusion==0)
Tx_1 <- Defervescence%>%filter(CefeExtendedInfusion==1)
quantile(Tx_0$TimetoDefervescence)
quantile(Tx_1$TimetoDefervescence)
t.test(Defervescence$TimetoDefervescence ~ Defervescence$CefeExtendedInfusion, var.equal = TRUE)

SIRS_only <- Defervescence%>%filter(SIRS_CriteriaMet==1)
SIRS_only
t.test(SIRS_only$TimetoDefervescence ~ SIRS_only$CefeExtendedInfusion, var.equal = TRUE)

Defervescence%>%group_by(CefeExtendedInfusion,BrokeFever_168)%>%
  summarize(count=n())
CrossTable(Defervescence$CefeExtendedInfusion, Defervescence$BrokeFever_168)
chisq.test(Defervescence$CefeExtendedInfusion, Defervescence$BrokeFever_168)

Defervescence%>%group_by(CefeExtendedInfusion,BrokeFever_144)%>%
  summarize(count=n())
CrossTable(Defervescence$CefeExtendedInfusion, Defervescence$BrokeFever_144)
chisq.test(Defervescence$CefeExtendedInfusion, Defervescence$BrokeFever_144)

Defervescence%>%group_by(CefeExtendedInfusion,BrokeFever_120)%>%
  summarize(count=n())
CrossTable(Defervescence$CefeExtendedInfusion, Defervescence$BrokeFever_120)
chisq.test(Defervescence$CefeExtendedInfusion, Defervescence$BrokeFever_120)

Defervescence%>%group_by(CefeExtendedInfusion,BrokeFever_96)%>%
  summarize(count=n())
CrossTable(Defervescence$CefeExtendedInfusion, Defervescence$BrokeFever_96)
chisq.test(Defervescence$CefeExtendedInfusion, Defervescence$BrokeFever_96)

Defervescence%>%group_by(CefeExtendedInfusion,BrokeFever_72)%>%
  summarize(count=n())
CrossTable(Defervescence$CefeExtendedInfusion, Defervescence$BrokeFever_72)
chisq.test(Defervescence$CefeExtendedInfusion, Defervescence$BrokeFever_72)

Defervescence%>%group_by(CefeExtendedInfusion,BrokeFever_48)%>%
  summarize(count=n())
CrossTable(Defervescence$CefeExtendedInfusion, Defervescence$BrokeFever_48)
chisq.test(Defervescence$CefeExtendedInfusion, Defervescence$BrokeFever_48)

Defervescence%>%group_by(CefeExtendedInfusion,BrokeFever_24)%>%
  summarize(count=n())
CrossTable(Defervescence$CefeExtendedInfusion, Defervescence$BrokeFever_24)
chisq.test(Defervescence$CefeExtendedInfusion, Defervescence$BrokeFever_24)


transform(Defervescence, percent = ave(BrokeFever_144, CefeExtendedInfusion, FUN = prop.table))


table(Defervescence$BrokeFever_168 ~ Defervescence$CefeExtendedInfusion)

median(Defervescence$CefeExtendedInfusion, Defervescence$TimetoDefervescence)

CrossTable(Defervescence$CefeExtendedInfusion, Defervescence$BrokeFever_144)
chisq.test(Defervescence$CefeExtendedInfusion, Defervescence$BrokeFever_144)

CrossTable(Defervescence$CefeExtendedInfusion, Defervescence$BrokeFever_72)
chisq.test(Defervescence$CefeExtendedInfusion, Defervescence$BrokeFever_72)


CrossTable(data$ReceivedTMP, data$X30Day_Readmission)
chisq.test(data$ReceivedTMP, data$X30Day_Readmission)

CrossTable(data$ReceivedTMP, data$CarbapenemBeforeIndex)
chisq.test(data$ReceivedTMP, data$CarbapenemBeforeIndex)

CrossTable(data$ReceivedTMP, data$AntimicrobialBeforeIndex)
chisq.test(data$ReceivedTMP, data$AntimicrobialBeforeIndex)

mean(data$AvgAdjTMPDoseBodyWeight)

mean(data$DoxyDays, rm.na = TRUE)
mean(data$CeftazidimeDays)
mean(data$LevofloxacinDays)
mean(data$TMPDays)

describe.by(data$Received_Doxycycline)


# Survival Curves for each time point of defervescence 
require("survival")
fit <- survfit(Surv(Defervescence$TimetoDefervescence)~Defervescence$CefeExtendedInfusion, data=Defervescence)
ggsurvplot(
  fit,
  data=Defervescence, 
  size=1, 
  palette=
    c("#355BB7", "#EE7D36"), 
  xlab="Time to defervescence (hours)",
  ylab="Proportion of patients remaining febrile",
  title="Defervescence by 168 hours",
  surv.median.line = "hv",
  break.x.by = 24,
  xlim = c(0,168),
  pval=TRUE,
  risk.table=TRUE, 
  risk.tablel.col="strata",
  legend.labs=
    c("Standard Infusion", "Extended Infusion"), 
  risk.table.height=0.25,
  risk.table.y.text = FALSE, 
  ggtheme = theme_bw()
)
summary(fit)
survdiff(Surv(Defervescence$TimetoDefervescence)~Defervescence$CefeExtendedInfusion, rho=0)
survdiff(Surv(Defervescence$TimetoDefervescence)~Defervescence$CefeExtendedInfusion, rho=1)
coxph(Surv(Defervescence$TimetoDefervescence)~Defervescence$CefeExtendedInfusion, data=Defervescence)


fit24 <- survfit(Surv(Defervescence$TimetoDefervescence, Defervescence$BrokeFever_24)~Defervescence$CefeExtendedInfusion, data=Defervescence)
plot2<-ggsurvplot(
  fit24,
  data=Defervescence, 
  size=1, 
  palette=
    c("#355BB7", "#EE7D36"), 
  xlab="Time to defervescence (hours)",
  ylab="Proportion of patients remaining febrile",
  title="Defervescence by 24 hours",
  surv.median.line = "hv",
  break.x.by = 24,
  xlim = c(0,24),
  pval=TRUE,
  risk.table=TRUE, 
  risk.tablel.col="strata",
  legend.labs=
    c("Standard Infusion", "Extended Infusion"), 
  risk.table.height=0.25,
  risk.table.y.text = FALSE,
  ggtheme = theme_bw()
)
summary(fit24)
survdiff(Surv(Defervescence$TimetoDefervescence, Defervescence$BrokeFever_24)~Defervescence$CefeExtendedInfusion, rho=0)
coxph(Surv(Defervescence$TimetoDefervescence, Defervescence$BrokeFever_24)~Defervescence$CefeExtendedInfusion, data=Defervescence)


fit48 <- survfit(Surv(Defervescence$TimetoDefervescence, Defervescence$BrokeFever_48)~Defervescence$CefeExtendedInfusion, data=Defervescence)
plot3<-ggsurvplot(
  fit48,
  data=Defervescence, 
  size=1, 
  palette=
    c("#355BB7", "#EE7D36"), 
  xlab="Time to defervescence (hours)",
  ylab="Proportion of patients remaining febrile",
  title="Defervescence by 48 hours",
  surv.median.line = "hv",
  break.x.by = 24,
  xlim = c(0,48),
  pval=TRUE,
  risk.table=TRUE, 
  risk.tablel.col="strata",
  legend.labs=
    c("Standard Infusion", "Extended Infusion"), 
  risk.table.height=0.25,
  risk.table.y.text = FALSE,
  ggtheme = theme_bw()
)
summary(fit48)
survdiff(Surv(Defervescence$TimetoDefervescence, Defervescence$BrokeFever_48)~Defervescence$CefeExtendedInfusion, rho=0)
coxph(Surv(Defervescence$TimetoDefervescence, Defervescence$BrokeFever_48)~Defervescence$CefeExtendedInfusion, data=Defervescence)


fit72 <- survfit(Surv(Defervescence$TimetoDefervescence, Defervescence$BrokeFever_72)~Defervescence$CefeExtendedInfusion, data=Defervescence)
plot4<-ggsurvplot(
  fit72,
  data=Defervescence, 
  size=1, 
  palette=
    c("#355BB7", "#EE7D36"), 
  xlab="Time to defervescence (hours)",
  ylab="Proportion of patients remaining febrile",
  title="Defervescence by 72 hours",
  surv.median.line = "hv",
  break.x.by = 24,
  xlim = c(0,72),
  pval=TRUE,
  risk.table=TRUE, 
  risk.tablel.col="strata",
  legend.labs=
    c("Standard Infusion", "Extended Infusion"), 
  risk.table.height=0.25,
  risk.table.y.text = FALSE,
  ggtheme = theme_bw()
)
summary(fit72)
survdiff(Surv(Defervescence$TimetoDefervescence, Defervescence$BrokeFever_72)~Defervescence$CefeExtendedInfusion, rho=0)
coxph(Surv(Defervescence$TimetoDefervescence, Defervescence$BrokeFever_72)~Defervescence$CefeExtendedInfusion, data=Defervescence)


fit96 <- survfit(Surv(Defervescence$TimetoDefervescence, Defervescence$BrokeFever_96)~Defervescence$CefeExtendedInfusion, data=Defervescence)
ggsurvplot(
  fit96,
  data=Defervescence, 
  size=1, 
  palette=
    c("#355BB7", "#EE7D36"), 
  xlab="Time to defervescence (hours)",
  ylab="Proportion of patients remaining febrile",
  title="Defervescence by 96 hours",
  surv.median.line = "hv",
  break.x.by = 24,
  xlim = c(0,96),
  pval=TRUE,
  risk.table=TRUE, 
  risk.tablel.col="strata",
  legend.labs=
    c("Standard Infusion", "Extended Infusion"), 
  risk.table.height=0.25,
  risk.table.y.text = FALSE,
  ggtheme = theme_bw()
)
summary(fit96)
survdiff(Surv(Defervescence$TimetoDefervescence, Defervescence$BrokeFever_96)~Defervescence$CefeExtendedInfusion, rho=0)
coxph(Surv(Defervescence$TimetoDefervescence, Defervescence$BrokeFever_96)~Defervescence$CefeExtendedInfusion, data=Defervescence)


fit120 <- survfit(Surv(Defervescence$TimetoDefervescence, Defervescence$BrokeFever_120)~Defervescence$CefeExtendedInfusion, data=Defervescence)
ggsurvplot(
  fit120,
  data=Defervescence, 
  size=1, 
  palette=
    c("#355BB7", "#EE7D36"), 
  xlab="Time to defervescence (hours)",
  ylab="Proportion of patients remaining febrile",
  title="Defervescence by 120 hours",
  surv.median.line = "hv",
  break.x.by = 24,
  xlim = c(0,120),
  pval=TRUE,
  risk.table=TRUE, 
  risk.tablel.col="strata",
  legend.labs=
    c("Standard Infusion", "Extended Infusion"), 
  risk.table.height=0.25,
  risk.table.y.text = FALSE,
  ggtheme = theme_bw()
)
summary(fit120)
survdiff(Surv(Defervescence$TimetoDefervescence, Defervescence$BrokeFever_120)~Defervescence$CefeExtendedInfusion, rho=0)
coxph(Surv(Defervescence$TimetoDefervescence, Defervescence$BrokeFever_120)~Defervescence$CefeExtendedInfusion, data=Defervescence)


fit144 <- survfit(Surv(Defervescence$TimetoDefervescence, Defervescence$BrokeFever_144)~Defervescence$CefeExtendedInfusion, data=Defervescence)
ggsurvplot(
  fit144,
  data=Defervescence, 
  size=1, 
  palette=
    c("#355BB7", "#EE7D36"), 
  xlab="Time to defervescence (hours)",
  ylab="Proportion of patients remaining febrile",
  title="Defervescence by 144 hours",
  surv.median.line = "hv",
  break.x.by = 24,
  xlim = c(0,144),
  pval=TRUE,
  risk.table=TRUE, 
  risk.tablel.col="strata",
  legend.labs=
    c("Standard Infusion", "Extended Infusion"), 
  risk.table.height=0.25,
  risk.table.y.text = FALSE,
  ggtheme = theme_bw()
)
summary(fit144)
survdiff(Surv(Defervescence$TimetoDefervescence, Defervescence$BrokeFever_144)~Defervescence$CefeExtendedInfusion, rho=0)
coxph(Surv(Defervescence$TimetoDefervescence, Defervescence$BrokeFever_144)~Defervescence$CefeExtendedInfusion, data=Defervescence)


fit168 <- survfit(Surv(Defervescence$TimetoDefervescence, Defervescence$BrokeFever_168)~Defervescence$CefeExtendedInfusion, data=Defervescence)
plot1<-ggsurvplot(
  fit168,
  data=Defervescence, 
  size=1, 
  palette=
    c("#355BB7", "#EE7D36"), 
  xlab="Time to defervescence (hours)",
  ylab="Proportion of patients remaining febrile",
  title="Defervescence by 168 hours",
  surv.median.line = "hv",
  break.x.by = 24,
  xlim = c(0,168),
  pval=TRUE,
  risk.table=TRUE, 
  risk.tablel.col="strata",
  legend.labs=
    c("Standard Infusion", "Extended Infusion"), 
  risk.table.height=0.25,
  risk.table.y.text = FALSE,
  ggtheme = theme_bw()
)
summary(fit168)
survdiff(Surv(Defervescence$TimetoDefervescence, Defervescence$BrokeFever_168)~Defervescence$CefeExtendedInfusion, rho=0)
coxph(Surv(Defervescence$TimetoDefervescence, Defervescence$BrokeFever_168)~Defervescence$CefeExtendedInfusion, data=Defervescence)


###### Reviewer Comments Regarding Dates 
data$Dates <- as.Date(data$Dates, "%m/%d/%Y")
data$Date_Yr <- format(as.Date(data$Dates, format="%d/%m/%Y"),"%Y")

# using time period just from 2016
data <- data%>%
  mutate(TimePeriod = case_when(Date_Yr <= 2016 ~ 0,
                                Date_Yr > 2016 ~ 1))
CrossTable(data$TimePeriod, data$CefeExtendedInfusion)
chisq.test(data$TimePeriod, data$CefeExtendedInfusion)

# using time period from april of 2016
data <- data%>%
  mutate(Accurate_TimePeriod = case_when(Dates <= "2016-04-01" ~ 0,
                                   Dates > "2016-04-01" ~ 1))
CrossTable(data$Accurate_TimePeriod, data$CefeExtendedInfusion)
chisq.test(data$Accurate_TimePeriod, data$CefeExtendedInfusion)


SI<-data%>%filter(CefeExtendedInfusion==0)
EI<-data%>%filter(CefeExtendedInfusion==1)
range(SI$Dates)
range(EI$Dates)

# 4 figures arranged in 2 rows and 2 columns
cowplot::plot_grid(plot1$plot, plot2$plot, plot3$plot, plot4$plot,
                   nrow=2, ncol=2, labels = c('A', 'B', 'C', 'D'))



###### Hierarchal Linear Regression ANalysis 

a<-lm(TimetoDefervescence ~ CefeExtendedInfusion, data=data)
summary(a)

b<-lm(TimetoDefervescence ~ CefeExtendedInfusion + Accurate_TimePeriod, data=data)
summary(b)

c<-lm(TimetoDefervescence ~ CefeExtendedInfusion + BaseLDH + BaseProcalcitonin + SIRS_CriteriaMet + AdminGCSF + ReceivedCefeLD, data=data) 
summary(c)

d<-lm(TimetoDefervescence ~ CefeExtendedInfusion + BaseLDH + BaseProcalcitonin + SIRS_CriteriaMet + AdminGCSF + ReceivedCefeLD + Accurate_TimePeriod, data=data) 
summary(d)

reviewer_remove<-lm(TimetoDefervescence ~ CefeExtendedInfusion + SIRS_CriteriaMet , data=data) #+ ReceivedCefeLD
summary(reviewer_remove)

reviewer_remove1<-lm(TimetoDefervescence ~ CefeExtendedInfusion + SIRS_CriteriaMet  + Accurate_TimePeriod, data=data) #+ ReceivedCefeLD 
summary(reviewer_remove1)

e<-lm(BrokeFever_72 ~ CefeExtendedInfusion, data=Defervescence)
summary(e)

f<-lm(BrokeFever_72 ~ CefeExtendedInfusion + Accurate_TimePeriod, data=Defervescence)
summary(f)

g<-lm(BrokeFever_72 ~ CefeExtendedInfusion + BaseLDH + BaseProcalcitonin + SIRS_CriteriaMet + AdminGCSF + ReceivedCefeLD, data=Defervescence) 
summary(g)

h<-lm(BrokeFever_72 ~ CefeExtendedInfusion + BaseLDH + BaseProcalcitonin + SIRS_CriteriaMet + AdminGCSF + ReceivedCefeLD + Accurate_TimePeriod, data=Defervescence) 
summary(h)


i<-lm(BrokeFever_48 ~ CefeExtendedInfusion, data=Defervescence)
summary(i)

j<-lm(BrokeFever_48 ~ CefeExtendedInfusion + Accurate_TimePeriod, data=Defervescence)
summary(j)

k<-lm(BrokeFever_48 ~ CefeExtendedInfusion + BaseLDH + BaseProcalcitonin + SIRS_CriteriaMet + AdminGCSF + ReceivedCefeLD, data=Defervescence) 
summary(k)

l<-lm(BrokeFever_48 ~ CefeExtendedInfusion + BaseLDH + BaseProcalcitonin + SIRS_CriteriaMet + AdminGCSF + ReceivedCefeLD + Accurate_TimePeriod, data=Defervescence) 
summary(l)


m<-lm(BrokeFever_24 ~ CefeExtendedInfusion, data=Defervescence)
summary(m)

n<-lm(BrokeFever_24 ~ CefeExtendedInfusion + Accurate_TimePeriod, data=Defervescence)
summary(n)

o<-lm(BrokeFever_24 ~ CefeExtendedInfusion + BaseLDH + BaseProcalcitonin + SIRS_CriteriaMet + AdminGCSF + ReceivedCefeLD, data=Defervescence) 
summary(o)

p<-lm(BrokeFever_24 ~ CefeExtendedInfusion + BaseLDH + BaseProcalcitonin + SIRS_CriteriaMet + AdminGCSF + ReceivedCefeLD + Accurate_TimePeriod, data=Defervescence) 
summary(p)



--------------------------
  
##### Adding TimePeriod (2014-2017 v 2016-2020)
fit24 <- survfit(Surv(Defervescence$TimetoDefervescence, Defervescence$BrokeFever_24)~Defervescence$CefeExtendedInfusion + Defervescence$Accurate_TimePeriod, data=Defervescence)
ggsurvplot(
  fit24,
  data=Defervescence, 
  size=1, 
  xlab="Time to defervescence (hours)",
  ylab="Proportion of patients remaining febrile",
  title="Defervescence by 24 hours",
  #surv.median.line = "hv",
  break.x.by = 24,
  xlim = c(0,24),
  pval=TRUE,
  risk.table=TRUE, 
  risk.tablel.col="strata",
  legend.labs=
    c("Standard Infusion/First Time Period", "Standard Infusion/Second Time Period", "Extended Infusion/First Time Period", "Extended Infusion/Second Time Period"), 
  risk.table.height=0.25,
  risk.table.y.text = TRUE,
  ggtheme = theme_bw()
)
summary(fit24)
survdiff(Surv(Defervescence$TimetoDefervescence, Defervescence$BrokeFever_24)~Defervescence$CefeExtendedInfusion + Defervescence$Accurate_TimePeriod, rho=0)


fit48 <- survfit(Surv(Defervescence$TimetoDefervescence, Defervescence$BrokeFever_48)~Defervescence$CefeExtendedInfusion + Defervescence$Accurate_TimePeriod, data=Defervescence)
ggsurvplot(
  fit48,
  data=Defervescence, 
  size=1, 
  xlab="Time to defervescence (hours)",
  ylab="Proportion of patients remaining febrile",
  title="Defervescence by 48 hours",
  surv.median.line = "hv",
  break.x.by = 24,
  xlim = c(0,48),
  pval=TRUE,
  risk.table=TRUE, 
  risk.tablel.col="strata",
  legend.labs=
    c("Standard Infusion/First Time Period", "Standard Infusion/Second Time Period", "Extended Infusion/First Time Period", "Extended Infusion/Second Time Period"), 
  risk.table.height=0.25,
  risk.table.y.text = TRUE,
  ggtheme = theme_bw()
)
summary(fit48)
survdiff(Surv(Defervescence$TimetoDefervescence, Defervescence$BrokeFever_48)~Defervescence$CefeExtendedInfusion + Defervescence$Accurate_TimePeriod, rho=0)

fit72 <- survfit(Surv(TimetoDefervescence, BrokeFever_72)~CefeExtendedInfusion + BaseLDH + BaseProcalcitonin + SIRS_CriteriaMet + AdminGCSF + ReceivedCefeLD + Accurate_TimePeriod, data=Defervescence)
ggsurvplot(
  fit72,
  data=Defervescence, 
  size=1, 
  xlab="Time to defervescence (hours)",
  ylab="Proportion of patients remaining febrile",
  title="Defervescence by 72 hours",
  surv.median.line = "hv",
  break.x.by = 24,
  xlim = c(0,72),
  pval=TRUE,
  risk.table=TRUE, 
  risk.tablel.col="strata",
  #legend.labs=
  #c("Standard Infusion/First Time Period", "Standard Infusion/Second Time Period", "Extended Infusion/First Time Period", "Extended Infusion/Second Time Period"), 
  risk.table.height=0.25,
  risk.table.y.text = TRUE,
  ggtheme = theme_bw()
)
summary(fit72)
survdiff(Surv(Defervescence$TimetoDefervescence, Defervescence$BrokeFever_72)~Defervescence$CefeExtendedInfusion + Defervescence$BaseLDH + Defervescence$BaseProcalcitonin + Defervescence$SIRS_CriteriaMet + Defervescence$AdminGCSF + Defervescence$ReceivedCefeLD + Defervescence$Accurate_TimePeriod,  rho=0)


fit96 <- survfit(Surv(Defervescence$TimetoDefervescence, Defervescence$BrokeFever_96)~Defervescence$CefeExtendedInfusion + Defervescence$Accurate_TimePeriod, data=Defervescence)
ggsurvplot(
  fit96,
  data=Defervescence, 
  size=1, 
  xlab="Time to defervescence (hours)",
  ylab="Proportion of patients remaining febrile",
  title="Defervescence by 96 hours",
  surv.median.line = "hv",
  break.x.by = 24,
  xlim = c(0,96),
  pval=TRUE,
  risk.table=TRUE, 
  risk.tablel.col="strata",
  legend.labs=
    c("Standard Infusion/First Time Period", "Standard Infusion/Second Time Period", "Extended Infusion/First Time Period", "Extended Infusion/Second Time Period"), 
  risk.table.height=0.25,
  risk.table.y.text = TRUE,
  ggtheme = theme_bw()
)
summary(fit96)
survdiff(Surv(Defervescence$TimetoDefervescence, Defervescence$BrokeFever_96)~Defervescence$CefeExtendedInfusion + Defervescence$Accurate_TimePeriod, rho=0)


fit120 <- survfit(Surv(Defervescence$TimetoDefervescence, Defervescence$BrokeFever_120)~Defervescence$CefeExtendedInfusion + Defervescence$Accurate_TimePeriod, data=Defervescence)
ggsurvplot(
  fit120,
  data=Defervescence, 
  size=1, 
  xlab="Time to defervescence (hours)",
  ylab="Proportion of patients remaining febrile",
  title="Defervescence by 120 hours",
  surv.median.line = "hv",
  break.x.by = 24,
  xlim = c(0,120),
  pval=TRUE,
  risk.table=TRUE, 
  risk.tablel.col="strata",
  legend.labs=
    c("Standard Infusion/First Time Period", "Standard Infusion/Second Time Period", "Extended Infusion/First Time Period", "Extended Infusion/Second Time Period"), 
  risk.table.height=0.25,
  risk.table.y.text = TRUE,
  ggtheme = theme_bw()
)
summary(fit120)
survdiff(Surv(Defervescence$TimetoDefervescence, Defervescence$BrokeFever_120)~Defervescence$CefeExtendedInfusion + Defervescence$Accurate_TimePeriod, rho=0)


fit144 <- survfit(Surv(Defervescence$TimetoDefervescence, Defervescence$BrokeFever_144)~Defervescence$CefeExtendedInfusion + Defervescence$Accurate_TimePeriod, data=Defervescence)
ggsurvplot(
  fit144,
  data=Defervescence, 
  size=1, 
  xlab="Time to defervescence (hours)",
  ylab="Proportion of patients remaining febrile",
  title="Defervescence by 144 hours",
  surv.median.line = "hv",
  break.x.by = 24,
  xlim = c(0,144),
  pval=TRUE,
  risk.table=TRUE, 
  risk.tablel.col="strata",
  legend.labs=
    c("Standard Infusion/First Time Period", "Standard Infusion/Second Time Period", "Extended Infusion/First Time Period", "Extended Infusion/Second Time Period"), 
  risk.table.height=0.25,
  risk.table.y.text = TRUE,
  ggtheme = theme_bw()
)
summary(fit144)
survdiff(Surv(Defervescence$TimetoDefervescence, Defervescence$BrokeFever_144)~Defervescence$CefeExtendedInfusion + Defervescence$Accurate_TimePeriod, rho=0)


fit168 <- survfit(Surv(Defervescence$TimetoDefervescence, Defervescence$BrokeFever_168)~Defervescence$CefeExtendedInfusion , data=Defervescence)
ggsurvplot(
  fit168,
  data=Defervescence, 
  size=1, 
  xlab="Time to defervescence (hours)",
  ylab="Proportion of patients remaining febrile",
  title="Defervescence by 168 hours",
  surv.median.line = "hv",
  break.x.by = 24,
  xlim = c(0,168),
  pval=TRUE,
  risk.table=TRUE, 
  risk.tablel.col="strata",
  #legend.labs=
  # c("Standard Infusion/First Time Period", "Standard Infusion/Second Time Period", "Extended Infusion/First Time Period", "Extended Infusion/Second Time Period"), 
  risk.table.height=0.25,
  risk.table.y.text = TRUE,
  ggtheme = theme_bw()
)
summary(fit168)
survdiff(Surv(Defervescence$TimetoDefervescence, Defervescence$BrokeFever_168)~Defervescence$CefeExtendedInfusion + Defervescence$TimePeriod, rho=0)
survdiff(Surv(Defervescence$TimetoDefervescence, Defervescence$BrokeFever_168)~Defervescence$CefeExtendedInfusion + Defervescence$TimePeriod, rho=1)


