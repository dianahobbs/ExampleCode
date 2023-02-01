#### Steno data project 
### Setting up the data and running analyses
## Author: Diana Hobbs
## Date: April 2022

#### Set Up Packages and Load Data ####
pacman::p_load(tidyverse,jtools,psych,cowplot,multcomp,gmodels)

data <- read.csv("./StenoData.csv", header = TRUE)

data <- data%>%
  rename(Age=Age..yr..at.time.of.admission, Sex=Gender, ICU=ICU.Admin..Y.1., MechVent=Mech.Vent..Y.1., ReceivedTMP=Received.TMP..Y.1., APACHE=APACHE.II.Score,
         CCI=CCI.Score, AddSteno48=Added.steno.abx..ceftaz..FQ..TMP.SMX..TC..after.48h.on.origil.abx.that.covered.but.within.10.days.of.originsl.abx..Y.N, 
         Mortality30Day=X30.day.mortality..y.1., Mortality31to90=X31.90.day.mortality..y.1., CarbapenemBeforeIndex=Carbapenem.b4.index, AntimicrobialBeforeIndex=Broad.GN.coverage.started.b4.index, 
         BaselineCreatineClearance_Female_IBW=Baseline.Creatine.Clearance..IBW.........................eCrCl.mL.min...Female,
         BaselineCreatineClearance_Male_IBW=Baseline.Creatine.Clearance..IBW.........................eCrCl.mL.min...Male,
         BaselineCreatineClearance_Female_ADJ=Baseline.Creatine.Clearance..AdjBW.........................eCrCl.mL.min..Female,
         BaselineCreatineClearance_Male_ADJ=Baseline.Creatine.Clearance..AdjBW.........................eCrCl.mL.min..Male,
         DoxyHrs=Doxy.duration.hrs, CeftazidimeHrs=Ceftazidime.duration.hrs, LevofloxacinHrs=Levo.duration.hrs, MinocyclineHrs=Minocycline.duration.hrs, Readmission30Day=Re.admitted.within.30.days.....Yes.1.,
         CiprHrs=Cipro.duration.hrs, TMPHrs=TMP.duration..hrs., TreatmentFail=Treatment.failure..Yes..1, AddSteno72=Steno.on.repeat.cx...72h.from.index...Y.N...non.taken., 
         InHouseMortality=In.house.mortality..Y.1..n.0, ActualBodyWeight_kg=Actual.Body.Weight..kg., IBW_kg=IBW..kg., AdjBW_kg=AdjBW..kg., TMPDoseAdjMgKgD=Adj.body.wt.mg.kg.d)%>%
  replace_na(list(DoxyHrs=0, CeftazidimeHrs=0, LevofloxacinHrs=0, MinocyclineHrs=0, CiprHrs=0, TMPHrs=0, TreatmentFail=0, AddSteno72=0, InHouseMortality=0,
                  BaselineCreatineClearance_Female_IBW=0, BaselineCreatineClearance_Male_IBW=0, BaselineCreatineClearance_Female_ADJ=0, BaselineCreatineClearance_Male_ADJ=0))%>%
         # Create the maximum amount of hours on alternative treatment
  mutate(AltTxHrs = pmax(DoxyHrs,CeftazidimeHrs,LevofloxacinHrs,MinocyclineHrs,CiprHrs),
         # max hours on TMP vs alt so we can remove anyone with less than 48 hrs worth of treatment
         TMPvALT_Duration_hrs = pmax(TMPHrs,AltTxHrs),
         # orig composite = added drug (tx fail) & positive repeat respiratory culture 72 hrs after index &/or In house mortality
         OriginalOutcome = case_when(TreatmentFail==0 & AddSteno72==0 & InHouseMortality==0 ~ 0,
                                     TreatmentFail==1 | AddSteno72==1 | InHouseMortality==1 ~ 1),
         # new composite = positive repeat respiratory culture 72 hrs after index &/or In house mortality
         NewPrimaryOutcome = case_when(AddSteno72==0 & InHouseMortality==0 ~ 0,
                                       AddSteno72==1 | InHouseMortality==1 ~ 1),
         # Combine two columns to get IBW and ADJ CC
         CreatineClearance_IBW = BaselineCreatineClearance_Female_IBW + BaselineCreatineClearance_Male_IBW,
         CreatineClearance_ADJ = BaselineCreatineClearance_Female_ADJ + BaselineCreatineClearance_Male_ADJ,
         # Adjust the weight so we can figure out correct TMP dose to use 
         Weight_to_adjust_TMPdose = case_when(ActualBodyWeight_kg > (1.2*IBW_kg) ~ AdjBW_kg,
                                              ActualBodyWeight_kg < IBW_kg ~ ActualBodyWeight_kg,
                                              (ActualBodyWeight_kg > IBW_kg) & (ActualBodyWeight_kg < (1.2*IBW_kg)) ~ ActualBodyWeight_kg),
         # Calculations to get the daily dose mg/d
         TMPDailyDose = TMPDoseAdjMgKgD * AdjBW_kg,
         # Figure out what TMP dose is based on above weight (will also need to edit for creatine clearance)
         TMPdose_by_weight = ((TMPDailyDose)/Weight_to_adjust_TMPdose),
         # Figure out what TMP dose is based on creatine clearance now 
         TMPdose_by_creatine = case_when(
           CreatineClearance_ADJ > 30 ~ TMPdose_by_weight,
           CreatineClearance_ADJ < 15 ~ (TMPdose_by_weight + (TMPdose_by_weight*0.66)),
           CreatineClearance_ADJ > 15 & CreatineClearance_ADJ < 30 ~ (TMPdose_by_weight + (TMPdose_by_weight*0.33))),
         HighvLowDose = case_when(ReceivedTMP==1 & TMPdose_by_creatine>=10 ~ 1,
                                  ReceivedTMP==1 & TMPdose_by_creatine<10 ~ 0))%>%
  dplyr::select(Age, Sex, ICU, MechVent, TreatmentFail, AddSteno48, AddSteno72, InHouseMortality, Mortality30Day, Mortality31to90, Readmission30Day, 
                OriginalOutcome, NewPrimaryOutcome, APACHE, CCI, CarbapenemBeforeIndex, AntimicrobialBeforeIndex,  AltTxHrs, ReceivedTMP, TMPHrs, 
                TMPvALT_Duration_hrs, ActualBodyWeight_kg, AdjBW_kg, IBW_kg, CreatineClearance_IBW, CreatineClearance_ADJ, Weight_to_adjust_TMPdose, 
                TMPdose_by_weight, TMPdose_by_creatine, HighvLowDose)
  
data1 <- data%>%filter(TMPvALT_Duration_hrs >=48)
data2 <- data%>%filter(TMPvALT_Duration_hrs >48)
       
         
#Create Z score for TmpDose based on Adj & clear & remove 3SD above mean
data <- data %>%
  mutate(., zscore = (TMPdose_by_creatine - mean(TMPdose_by_creatine, na.rm=TRUE))/sd(TMPdose_by_creatine, na.rm=TRUE))%>%
  filter(zscore <= "3" | is.na(zscore))

data1 <- data1 %>%
  mutate(., zscore = (TMPdose_by_creatine - mean(TMPdose_by_creatine, na.rm=TRUE))/sd(TMPdose_by_creatine, na.rm=TRUE))%>%
  filter(zscore <= "3" | is.na(zscore))

data2 <- data2 %>%
  mutate(., zscore = (TMPdose_by_creatine - mean(TMPdose_by_creatine, na.rm=TRUE))/sd(TMPdose_by_creatine, na.rm=TRUE))%>%
  filter(zscore <= "3" | is.na(zscore))         

##### Analyses #####

# Using the data that had >= 48 hr worth of treatment AND dropped 3 SD from mean adjusted dose 
rm(data,data2)
table(data1$ReceivedTMP)

CrossTable(data1$ReceivedTMP, data1$MechVent)
chisq.test(data1$ReceivedTMP, data1$MechVent)         

CrossTable(data1$ReceivedTMP, data1$TreatmentFail)
chisq.test(data1$ReceivedTMP, data1$TreatmentFail)

CrossTable(data1$ReceivedTMP, data1$AddSteno48)
chisq.test(data1$ReceivedTMP, data1$AddSteno48)

CrossTable(data1$ReceivedTMP, data1$AddSteno72)
chisq.test(data1$ReceivedTMP, data1$AddSteno72)

CrossTable(data1$ReceivedTMP, data1$InHouseMortality)
chisq.test(data1$ReceivedTMP, data1$InHouseMortality)


table(data1$HighvLowDose)

CrossTable(data1$HighvLowDose, data1$TreatmentFail)
chisq.test(data1$HighvLowDose, data1$TreatmentFail)

CrossTable(data1$HighvLowDose, data1$AddSteno48)
chisq.test(data1$HighvLowDose, data1$AddSteno48)

CrossTable(data1$HighvLowDose, data1$AddSteno72)
chisq.test(data1$HighvLowDose, data1$AddSteno72)

CrossTable(data1$HighvLowDose, data1$InHouseMortality)
chisq.test(data1$HighvLowDose, data1$InHouseMortality)

CrossTable(data1$HighvLowDose, data1$OriginalOutcome)
chisq.test(data1$HighvLowDose, data1$OriginalOutcome)

CrossTable(data1$HighvLowDose, data1$NewPrimaryOutcome)
chisq.test(data1$HighvLowDose, data1$NewPrimaryOutcome)

High <- data1%>%filter(HighvLowDose==1)
Low <- data1%>%filter(HighvLowDose==0)

mean(High$TMPdose_by_creatine)
sd(High$TMPdose_by_creatine)

mean(Low$TMPdose_by_creatine)
sd(Low$TMPdose_by_creatine)

data3 <- data1%>%mutate(ReceivedAltTx=case_when(AltTxHrs == 0 ~ 0, AltTxHrs > 0 ~ 1),
                        ReceivedOnlyTMP=case_when(ReceivedAltTx == 0 & ReceivedTMP == 1 ~ 1,
                                                  ReceivedTMP == 0 ~ 0,
                                                  ReceivedAltTx == 1 & ReceivedTMP == 1 ~ 2),
                        OnlyTMPvsTMPandAlt = case_when(ReceivedOnlyTMP==0 ~ 0, # did not receive TMP
                                                       ReceivedOnlyTMP==2 ~ 1)) # Received both TMP and Alt 
CrossTable(data3$ReceivedAltTx ~ data3$ReceivedOnlyTMP)
chisq.test(data3$ReceivedAltTx ~ data3$ReceivedOnlyTMP)
