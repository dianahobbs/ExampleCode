#### Propensity Score
### Setting up the data and running analyses
## Author: Diana Hobbs
## Date: November 2022

#### Set Up Packages and Load Data ####
pacman::p_load(tidyverse,dplyr,tableone,MatchIt,lmtest, sandwich, cobalt,rbounds,binarysens,apa,ggrepel)
remotes::install_github("ddsjoberg/gtsummary", build=F)
library("gtsummary")
 
data1 <- read.csv("./data.csv", header = TRUE)
data11 <- read.csv("./RecUTI_Success_Fig2.csv")
data1 <- merge(data1, data11, by="record_id", all=T)

data2 <- read.csv("./clinres_11.30.csv", header = TRUE)%>%
  select(record_id, Time.to.Clinical.Resolution..hrs.)
data <- merge(data1, data2, by="record_id"); rm(data1,data11,data2)

data3 <- read.csv("./Table1_12.6.csv", header = TRUE)%>%
  rename(Polyuria = Polyuria.urgency, CVA = CVA.flank.pain, Suprapubic = Suprapubic.abd.pain,
         K.pneumoniae=K..pneumoniae, P.mirabilis=P..mirabilis)%>%
  select(record_id, Dysuria, Polyuria, Hematuria, CVA, Suprapubic, Hypotension, K.pneumoniae, P.mirabilis)
data <- merge(data, data3, by="record_id"); rm(data3)

data <- data%>%mutate(Race.Caucasian = case_when(Race == 1 ~ "Caucasian",
                                                 Race != 1 ~ "Not Caucasian"),
                      Race.AA = case_when(Race == 2 ~ "African American",
                                          Race != 2 ~ "Not African American"),
                      Race.Other = case_when(Race == 3 ~ "Other",
                                             Race != 3 ~ "Not Other"))%>%
               rename(Sex = Male,
                      Arm = Treatment.Arm,
                      ICU = ICU.admission,
                      CCI = CCI.score, 
                      UncompUTI = Uncomplicated.UTI, 
                      CompUTI = Complicated.UTI,
                      Catheter = Catheter.related,
                      Pathogens = Number.of.urinary.pathogens,
                      Symptoms = Symptoms.Present,
                      ClinSuccess = Primary...Clin.Success,
                      SympRes_48 = Symptom.or.WBC.resolution..48.hrs,
                      TempRes_48 = Temp.resolution...48.hrs,
                      AbsentUTI_6mo = Absence.of.6mo.ESBL.UTI,
                      Mortality = In.hospital.mortality, 
                      Mortality_30 = X30d.Mortality,
                      HospLOS = Hospital.LOS..days.,
                      ICULOS = ICU.LOS..days.,
                      Time_to_ClinRes = Time.to.Clinical.Resolution..hrs.,
                      UrologicAbnormal = Urologic.Abnormalities)%>%
                mutate(Arm = case_when(Arm==1~0, Arm==2~1))%>%                  
                dplyr::select(Arm, Age, Sex, Race.Caucasian, Race.AA, Race.Other, ICU, CCI, APACHE, UncompUTI, CompUTI,Pyelonephritis, Catheter, Pathogens, Bacteremia, Symptoms, # End of covariates
                        ClinSuccess, SympRes_48, TempRes_48,AbsentUTI_6mo, Mortality, Mortality_30, HospLOS, ICULOS, Time_to_ClinRes, # End of outcomes
                        UrologicAbnormal, Immunocompromised, Dysuria, Polyuria, Hematuria, CVA, Suprapubic, Hypotension, E.coli, K.pneumoniae, P.mirabilis, RecurrentUTI_Success)

##### STATISTICS ON ORIGINAL DATASET #####
data1 <- data%>%mutate(Sex = case_when(Sex==1~"Male",
                                       Sex==0~"Female"),
                       Arm = case_when(Arm==0~"Carbapenem",
                                       Arm==1~"PipTazo"),
                       ICU = case_when(ICU==0~"no",
                                       ICU==1~"yes"),
                       UncompUTI = case_when(UncompUTI==0~"no",
                                             UncompUTI==1~"yes"),
                       CompUTI = case_when(CompUTI==0~"no",
                                           CompUTI==1~"yes"),
                       Pyelonephritis = case_when(Pyelonephritis==0~"no",
                                                  Pyelonephritis==1~"yes"),
                       Catheter = case_when(Catheter==0~"no",
                                            Catheter==1~"yes"),
                       Bacteremia = case_when(Bacteremia==0~"no",
                                              Bacteremia==1~"yes"),
                       Symptoms = case_when(Symptoms==0~"no",
                                            Symptoms==1~"yes"),
                       ClinSuccess = case_when(ClinSuccess==0~"no",
                                               ClinSuccess==1~"yes"),
                       SympRes_48 = case_when(SympRes_48==0~"no",
                                              SympRes_48==1~"yes"),
                       TempRes_48 = case_when(TempRes_48==0~"no",
                                              TempRes_48==1~"yes"),
                       AbsentUTI_6mo = case_when(AbsentUTI_6mo==0~"no",
                                                 AbsentUTI_6mo==1~"yes"),
                       Mortality = case_when(Mortality==0~"no",
                                             Mortality==1~"yes"),
                       Mortality_30 = case_when(Mortality_30==0~"no",
                                                Mortality_30==1~"yes"),
                       Pathogens = as.numeric(Pathogens),
                       UrologicAbnormal = case_when(UrologicAbnormal==0~"no",
                                                    UrologicAbnormal==1~"yes"),
                       Immunocompromised = case_when(Immunocompromised==0~"no",
                                                     Immunocompromised==1~"yes"),
                       Dysuria = case_when(Dysuria==0~"no",
                                           Dysuria==1~"yes"),
                       Polyuria = case_when(Polyuria==0~"no",
                                            Polyuria==1~"yes"),
                       Hematuria = case_when(Hematuria==0~"no",
                                             Hematuria==1~"yes"),
                       CVA = case_when(CVA==0~"no",
                                       CVA==1~"yes"),
                       Suprapubic = case_when(Suprapubic==0~"no",
                                              Suprapubic==1~"yes"),
                       Hypotension = case_when(Hypotension==0~"no",
                                               Hypotension==1~"yes"),
                       E.coli = case_when(E.coli==0~"no",
                                          E.coli==1~"yes"),
                       K.pneumoniae = case_when(K.pneumoniae==0~"no",
                                                K.pneumoniae==1~"yes"),
                       P.mirabilis = case_when(P.mirabilis==0~"no",
                                               P.mirabilis==1~"yes"),
                       RecurrentUTI_Success = case_when(RecurrentUTI_Success==0~"no",
                                                        RecurrentUTI_Success==1~"yes"))
data1$Pathogens <- as.numeric(as.character(data1$Pathogens))

# Descriptive Stats / demographic table
pacman::p_load(haven,gt,janitor,psych,apaTables,broom,modelsummary)

demo <- data1%>%
  dplyr::select(Arm, Sex, Race.Caucasian, Race.AA, Race.Other, Age, UrologicAbnormal, Immunocompromised, Bacteremia, Symptoms, Dysuria, Polyuria, 
                Hematuria, CVA, Suprapubic, Hypotension, ICU, CCI, APACHE, UncompUTI, CompUTI, Pyelonephritis, Catheter, E.coli, K.pneumoniae, P.mirabilis,
                ClinSuccess, SympRes_48, TempRes_48,AbsentUTI_6mo, Time_to_ClinRes, HospLOS, ICULOS, Mortality, Mortality_30, Pathogens)%>%
  tbl_summary(by=Arm,
              type = list(Pathogens = "continuous"),
              statistic=list(all_continuous()~ "{median}({sd})",
                             all_categorical()~"{n}({p}%)"),
              digits = all_continuous()~2,
              label = c(Sex = "Sex", 
                        Race.Caucasian = "Caucasian",
                        Race.AA = "African American",
                        Race.Other = "Other",
                        Age = "Age (yrs)", 
                        UrologicAbnormal = "Urologic Abnormalities",
                        Immunocompromised = "Immunocompromised",
                        Bacteremia = "Bacteremia", 
                        Symptoms = "Symptoms Present", 
                        Dysuria = "Dysuria", 
                        Polyuria = "Polyruia/urgency", 
                        Hematuria = "Hematuria", 
                        CVA = "CVA/flank pain", 
                        Suprapubic = "Suprapubic/abdominal pain", 
                        Hypotension = "Hypotension", 
                        ICU = "ICU Admission",
                        CCI = "CCI", 
                        APACHE = "APACHE", 
                        UncompUTI = "Uncomplicated UTI", 
                        CompUTI = "Complicated UTI",
                        Pyelonephritis = "Pyelonephritis", 
                        Catheter = "Catheter-Related", 
                        E.coli = "E.coli", 
                        K.pneumoniae = "K.pneumoniae", 
                        P.mirabilis = "P.mirabilis",
                        ClinSuccess = "Primary Clinical Success", 
                        SympRes_48 = "Symptom or WBC Resolution within 48 Hrs", 
                        TempRes_48 = "Temperature Resolution within 48 Hrs",
                        AbsentUTI_6mo = "Absence of ESBL UTI after 6 Months", 
                        Time_to_ClinRes = "Time to Clinical Resolution",
                        HospLOS = "Hospital Length of Stay (Days)", 
                        ICULOS = "ICU Length of Stay (Days)",
                        Mortality = "In-Hospital Mortality", 
                        Mortality_30 = "Mortality within 30 Days",
                        Pathogens ~ "Number of Urinary Pathogens"))%>%
  bold_labels()%>%
  italicize_levels()%>%
  #add_difference(conf.level = 0.95,pvalue_fun = ~style_pvalue(.x, digits = 2))
  add_p(pvalue_fun = ~style_pvalue(.x, digits = 2))
demo

# Table 1 Variables
a <- table(data1$Arm, data1$Sex)
a
chisq.test(data1$Arm, data1$Sex, correct=F)

b <- table(data1$Arm, data1$Race.Caucasian)
b
chisq.test(data1$Arm, data1$Race.Caucasian, correct=F)

c <- table(data1$Arm, data1$Race.AA)
c
chisq.test(data1$Arm, data1$Race.AA, correct=F)

d <- table(data1$Arm, data1$Race.Other)
d
chisq.test(data1$Arm, data1$Race.Other, correct=F)

wilcox.test(data1$Age ~ data1$Arm)
data1%>%group_by(Arm)%>%summarise(Median=median(Age), IQR_25 = quantile(Age, probs = 0.25), IQR_75 = quantile(Age, probs = 0.75))

e <- table(data1$Arm, data1$UrologicAbnormal)
e
chisq.test(data1$Arm, data1$UrologicAbnormal, correct=F)

f <- table(data1$Arm, data1$Immunocompromised)
f
chisq.test(data1$Arm, data1$Immunocompromised, correct=F)

g <- table(data1$Arm, data1$Bacteremia)
g
chisq.test(data1$Arm, data1$Bacteremia, correct=F)

h <- table(data1$Arm, data1$Symptoms)
h
chisq.test(data1$Arm, data1$Symptoms, correct=F)

i <- table(data1$Arm, data1$Dysuria)
i
chisq.test(data1$Arm, data1$Dysuria, correct=F)

j <- table(data1$Arm, data1$Polyuria)
j
chisq.test(data1$Arm, data1$Polyuria, correct=F)

k <- table(data1$Arm, data1$Hematuria)
k
chisq.test(data1$Arm, data1$Hematuria, correct=F)

l <- table(data1$Arm, data1$CVA)
l
chisq.test(data1$Arm, data1$CVA, correct=F)

m <- table(data1$Arm, data1$Suprapubic)
m
chisq.test(data1$Arm, data1$Suprapubic, correct=F)

n <- table(data1$Arm, data1$Hypotension)
n
chisq.test(data1$Arm, data1$Hypotension, correct=F)

o <- table(data1$Arm, data1$ICU)
o
chisq.test(data1$Arm, data1$ICU, correct=F)

wilcox.test(data1$CCI ~ data1$Arm)
data1%>%group_by(Arm)%>%summarise(Median=median(CCI), IQR_25 = quantile(CCI, probs = 0.25), IQR_75 = quantile(CCI, probs = 0.75))

wilcox.test(data1$APACHE ~ data1$Arm, na.rm=T)
data1%>%group_by(Arm)%>%summarise(Median=median(APACHE, na.rm=T), IQR_25 = quantile(APACHE, na.rm=T, probs = 0.25), IQR_75 = quantile(APACHE, na.rm=T, probs = 0.75))

p <- table(data1$Arm, data1$UncompUTI)
p
chisq.test(data1$Arm, data1$UncompUTI, correct=F)

q <- table(data1$Arm, data1$CompUTI)
q
chisq.test(data1$Arm, data1$CompUTI, correct=F)

r <- table(data1$Arm, data1$Pyelonephritis)
r
chisq.test(data1$Arm, data1$Pyelonephritis, correct=F)

s <- table(data1$Arm, data1$Catheter)
s
chisq.test(data1$Arm, data1$Catheter, correct=F)

t <- table(data1$Arm, data1$E.coli)
t
chisq.test(data1$Arm, data1$E.coli, correct=F)

u <- table(data1$Arm, data1$K.pneumoniae)
u
chisq.test(data1$Arm, data1$K.pneumoniae, correct=F)

v <- table(data1$Arm, data1$P.mirabilis)
v
chisq.test(data1$Arm, data1$P.mirabilis, correct=F)

# Table 3 Variables
pacman::p_load(catfun,epiR)
w <- table(data1$Arm, data1$ClinSuccess)
w
chisq.test(data1$Arm, data1$ClinSuccess, correct=F)
epi.2by2(w, conf.level = 0.95)

x <- table(data1$Arm, data1$SympRes_48)
x
chisq.test(data1$Arm, data1$SympRes_48, correct=F)
epi.2by2(x, conf.level = 0.95)

y <- table(data1$Arm, data1$TempRes_48)
y
chisq.test(data1$Arm, data1$TempRes_48, correct=F)
epi.2by2(y, conf.level = 0.95)

z <- table(data1$Arm, data1$AbsentUTI_6mo)
z
chisq.test(data1$Arm, data1$AbsentUTI_6mo, correct=F)
epi.2by2(z, conf.level = 0.95)

# Table 4 Variables
wilcox.test(data1$Time_to_ClinRes ~ data1$Arm, na.rm=T)
data1%>%group_by(Arm)%>%summarise(Median=median(Time_to_ClinRes, na.rm=T), IQR_25 = quantile(Time_to_ClinRes, na.rm=T, probs = 0.25), IQR_75 = quantile(Time_to_ClinRes, na.rm=T, probs = 0.75))

wilcox.test(data1$HospLOS ~ data1$Arm)
data1%>%group_by(Arm)%>%summarise(Median=median(HospLOS), IQR_25 = quantile(HospLOS, probs = 0.25), IQR_75 = quantile(HospLOS, probs = 0.75))

wilcox.test(data1$ICULOS ~ data1$Arm, na.rm=T)
data1%>%group_by(Arm)%>%summarise(Median=median(ICULOS, na.rm=T), IQR_25 = quantile(ICULOS, na.rm=T, probs = 0.25), IQR_75 = quantile(ICULOS, na.rm=T, probs = 0.75))

table(data1$Arm, data1$Mortality)
aa <- fisher.test(data1$Arm, data1$Mortality)
aa$p.value

table(data1$Arm, data1$Mortality_30)
bb <- fisher.test(data1$Arm, data1$Mortality_30)
bb$p.value


###### Compute covariance imbalances #####
# Std Difference 
treated <- data$Arm==1
cov <- data[,c("Age","Sex","ICU","CCI","UncompUTI",'CompUTI','Pyelonephritis','Catheter','Pathogens','Bacteremia','Symptoms')]
std.diff <-apply(cov,2,function(x) 100*(mean(x[treated])-mean(x[!treated]))/(sqrt(0.5*(var(x[treated])+var(x[!treated])))))
abs(std.diff) # ICU admission >25% and is suspect 

# Chi-square tests
pacman::p_load(RItools)
xBalance(Arm ~ Age + Sex + ICU + CCI + UncompUTI + CompUTI + Pyelonephritis + Catheter + Pathogens + Bacteremia + Symptoms,
         data=data, report = c("chisquare.test")) ## overall X^2(10)=10.4, p = 0.41
# Estimate the propensity score 
ps <- glm(Arm ~ Age + Sex + ICU + CCI + UncompUTI + CompUTI + Pyelonephritis + Catheter + Pathogens + Bacteremia + Symptoms, 
           data=data, family = "binomial")
summary(ps)
# Attach PS to datafile 
data$psvalue <- predict(ps, type="response")
#back to back histogram 
pacman::p_load(Hmisc)
histbackback(split(data$psvalue,
                   data$Arm),main="Propensity Score before Matching",
             xlab=c("carbapenem","piptazo"))

########## Match the treatment arms #####
# Match each carbapenem treatment with piptazo based on covariates
# Optimal matchin g
pacman::p_load(optmatch)
opt.match <- matchit(Arm ~ Age + Sex + ICU + CCI + UncompUTI + CompUTI + Pyelonephritis + Catheter + Pathogens + Bacteremia + Symptoms,
                     data = data,
                     method = "optimal",
                     distance = "logit", # specifies that generalized linear model calculates the PS based on all the covariates
                     ratio = 1, # specifies one-to-one matching is used (higher ratio = fewer discarded observations)
                     replace = F) # specifies matching is done without replacement 
summary(opt.match)
matched_data <- match.data(opt.match)
head(matched_data)

# Plot the imbalance between carbapenem and piptazo
plot(opt.match, type = "jitter", interactive = F)
plot(summary(opt.match), abs = F)

pacman::p_load(cobalt)
bal.tab(opt.match,m.threshold=0.1,un=T)
bal.tab(opt.match,v.threshold=2)
bal.plot(opt.match, 
         var.name="Age",
         which = "both",
         grid = F)
love.plot(bal.tab(opt.match, m.threshold=0.1),
          stat = "mean.diffs",
          grid = F, 
          stars = "raw",
          abs = F)

pacman::p_load(tableone)
table <- CreateTableOne(vars = c("Age",'Sex','ICU',"CCI",'UncompUTI','CompUTI','Pyelonephritis','Catheter','Pathogens','Bacteremia','Symptoms',
                                 'ClinSuccess','SympRes_48','TempRes_48','AbsentUTI_6mo','Mortality','Mortality_30'),
                        data=matched_data,
                        factorVars = c('Sex','ICU','UncompUTI','CompUTI','Pyelonephritis','Catheter','Bacteremia','Symptoms',
                                       'ClinSuccess',"SympRes_48",'TempRes_48','AbsentUTI_6mo','Mortality','Mortality_30'),
                        strata="Arm",
                        smd=T)

table <- print(table,smd=T,showAllLevels=T,noSpaces=T,printToggle=F)
write.csv(table,"./Tables/table_after_matching.csv")


# Compute indices of covariate imbalance after matching 
# Std Difference 
treated1 <- (matched_data$Arm==1)
cov1 <- data[,c("Age","Sex","ICU","CCI","UncompUTI",'CompUTI','Pyelonephritis','Catheter','Pathogens','Bacteremia','Symptoms')]
std.diff1 <-apply(cov1,2,function(x) 100*(mean(x[treated1])-mean(x[!treated1]))/(sqrt(0.5*(var(x[treated1])+var(x[!treated1])))))
abs(std.diff1) # none >25, none suspect anymore 
# Chi-square tests
xBalance(Arm ~ Age + Sex + ICU + CCI + UncompUTI + CompUTI + Pyelonephritis + Catheter + Pathogens + Bacteremia + Symptoms,
         data=matched_data, report = c("chisquare.test")) ## overall X^2(10)=2.02, p = 0.996
histbackback(split(matched_data$psvalue,
                   matched_data$Arm),main="Propensity Score after Matching",
             xlab=c("carbapenem","piptazo"))


###### STATS ON MATCHED DATASET #####
# Extract matched data
matched_data <- match.data(opt.match)
head(matched_data)

matched_data1 <- matched_data%>%mutate(Sex = case_when(Sex==1~"Male",
                                                       Sex==0~"Female"),
                                       Arm = case_when(Arm==0~"Carbapenem",
                                                       Arm==1~"PipTazo"),
                                       ICU = case_when(ICU==0~"no",
                                                       ICU==1~"yes"),
                                       UncompUTI = case_when(UncompUTI==0~"no",
                                                             UncompUTI==1~"yes"),
                                       CompUTI = case_when(CompUTI==0~"no",
                                                           CompUTI==1~"yes"),
                                       Pyelonephritis = case_when(Pyelonephritis==0~"no",
                                                                  Pyelonephritis==1~"yes"),
                                       Catheter = case_when(Catheter==0~"no",
                                                            Catheter==1~"yes"),
                                       Bacteremia = case_when(Bacteremia==0~"no",
                                                              Bacteremia==1~"yes"),
                                       Symptoms = case_when(Symptoms==0~"no",
                                                            Symptoms==1~"yes"),
                                       ClinSuccess = case_when(ClinSuccess==0~"no",
                                                               ClinSuccess==1~"yes"),
                                       SympRes_48 = case_when(SympRes_48==0~"no",
                                                              SympRes_48==1~"yes"),
                                       TempRes_48 = case_when(TempRes_48==0~"no",
                                                              TempRes_48==1~"yes"),
                                       AbsentUTI_6mo = case_when(AbsentUTI_6mo==0~"no",
                                                                 AbsentUTI_6mo==1~"yes"),
                                       Mortality = case_when(Mortality==0~"no",
                                                             Mortality==1~"yes"),
                                       Mortality_30 = case_when(Mortality_30==0~"no",
                                                                Mortality_30==1~"yes"),
                                       Pathogens = as.numeric(Pathogens),
                                       UrologicAbnormal = case_when(UrologicAbnormal==0~"no",
                                                                    UrologicAbnormal==1~"yes"),
                                       Immunocompromised = case_when(Immunocompromised==0~"no",
                                                                     Immunocompromised==1~"yes"),
                                       Dysuria = case_when(Dysuria==0~"no",
                                                           Dysuria==1~"yes"),
                                       Polyuria = case_when(Polyuria==0~"no",
                                                            Polyuria==1~"yes"),
                                       Hematuria = case_when(Hematuria==0~"no",
                                                             Hematuria==1~"yes"),
                                       CVA = case_when(CVA==0~"no",
                                                       CVA==1~"yes"),
                                       Suprapubic = case_when(Suprapubic==0~"no",
                                                              Suprapubic==1~"yes"),
                                       Hypotension = case_when(Hypotension==0~"no",
                                                               Hypotension==1~"yes"),
                                       E.coli = case_when(E.coli==0~"no",
                                                          E.coli==1~"yes"),
                                       K.pneumoniae = case_when(K.pneumoniae==0~"no",
                                                                K.pneumoniae==1~"yes"),
                                       P.mirabilis = case_when(P.mirabilis==0~"no",
                                                               P.mirabilis==1~"yes"),
                                       RecurrentUTI_Success = case_when(RecurrentUTI_Success==0~"no",
                                                                        RecurrentUTI_Success==1~"yes"))
matched_data1$Pathogens <- as.numeric(as.character(matched_data1$Pathogens))


demo_matched <- matched_data1%>%
  dplyr::select(Arm, Sex, Race.Caucasian, Race.AA, Race.Other, Age, UrologicAbnormal, Immunocompromised, Bacteremia, Symptoms, Dysuria, Polyuria, 
                Hematuria, CVA, Suprapubic, Hypotension, ICU, CCI, APACHE, UncompUTI, CompUTI, Pyelonephritis, Catheter, E.coli, K.pneumoniae, P.mirabilis,
                ClinSuccess, SympRes_48, TempRes_48,AbsentUTI_6mo, Time_to_ClinRes, HospLOS, ICULOS, Mortality, Mortality_30, Pathogens)%>%
  tbl_summary(by=Arm,
              type = list(Pathogens = "continuous"),
              statistic=list(all_continuous()~ "{median}({sd})",
                             all_categorical()~"{n}({p}%)"),
              digits = all_continuous()~2,
              label = c(Sex = "Sex", 
                        Race.Caucasian = "Caucasian",
                        Race.AA = "African American",
                        Race.Other = "Other",
                        Age = "Age (yrs)", 
                        UrologicAbnormal = "Urologic Abnormalities",
                        Immunocompromised = "Immunocompromised",
                        Bacteremia = "Bacteremia", 
                        Symptoms = "Symptoms Present", 
                        Dysuria = "Dysuria", 
                        Polyuria = "Polyruia/urgency", 
                        Hematuria = "Hematuria", 
                        CVA = "CVA/flank pain", 
                        Suprapubic = "Suprapubic/abdominal pain", 
                        Hypotension = "Hypotension", 
                        ICU = "ICU Admission",
                        CCI = "CCI", 
                        APACHE = "APACHE", 
                        UncompUTI = "Uncomplicated UTI", 
                        CompUTI = "Complicated UTI",
                        Pyelonephritis = "Pyelonephritis", 
                        Catheter = "Catheter-Related", 
                        E.coli = "E.coli", 
                        K.pneumoniae = "K.pneumoniae", 
                        P.mirabilis = "P.mirabilis",
                        ClinSuccess = "Primary Clinical Success", 
                        SympRes_48 = "Symptom or WBC Resolution within 48 Hrs", 
                        TempRes_48 = "Temperature Resolution within 48 Hrs",
                        AbsentUTI_6mo = "Absence of ESBL UTI after 6 Months", 
                        Time_to_ClinRes = "Time to Clinical Resolution",
                        HospLOS = "Hospital Length of Stay (Days)", 
                        ICULOS = "ICU Length of Stay (Days)",
                        Mortality = "In-Hospital Mortality", 
                        Mortality_30 = "Mortality within 30 Days",
                        Pathogens ~ "Number of Urinary Pathogens"))%>%
  bold_labels()%>%
  italicize_levels()%>%
  #add_difference(conf.level = 0.95,pvalue_fun = ~style_pvalue(.x, digits = 2))
  add_p(pvalue_fun = ~style_pvalue(.x, digits = 2))
demo_matched

eee# Table 1 Variables
a <- table(matched_data1$Arm, matched_data1$Sex)
a
chisq.test(matched_data1$Arm, matched_data1$Sex, correct=F)

b <- table(matched_data1$Arm, matched_data1$Race.Caucasian)
b
chisq.test(matched_data1$Arm, matched_data1$Race.Caucasian, correct=F)

c <- table(matched_data1$Arm, matched_data1$Race.AA)
c
chisq.test(matched_data1$Arm, matched_data1$Race.AA, correct=F)

d <- table(matched_data1$Arm, matched_data1$Race.Other)
d
chisq.test(matched_data1$Arm, matched_data1$Race.Other, correct=F)

wilcox.test(matched_data1$Age ~ matched_data1$Arm)
matched_data1%>%group_by(Arm)%>%summarise(Median=median(Age), IQR_25 = quantile(Age, probs = 0.25), IQR_75 = quantile(Age, probs = 0.75))

e <- table(matched_data1$Arm, matched_data1$UrologicAbnormal)
e
chisq.test(matched_data1$Arm, matched_data1$UrologicAbnormal, correct=F)

f <- table(matched_data1$Arm, matched_data1$Immunocompromised)
f
chisq.test(matched_data1$Arm, matched_data1$Immunocompromised, correct=F)

g <- table(matched_data1$Arm, matched_data1$Bacteremia)
g
chisq.test(matched_data1$Arm, matched_data1$Bacteremia, correct=F)

h <- table(matched_data1$Arm, matched_data1$Symptoms)
h
chisq.test(matched_data1$Arm, matched_data1$Symptoms, correct=F)

i <- table(matched_data1$Arm, matched_data1$Dysuria)
i
chisq.test(matched_data1$Arm, matched_data1$Dysuria, correct=F)

j <- table(matched_data1$Arm, matched_data1$Polyuria)
j
chisq.test(matched_data1$Arm, matched_data1$Polyuria, correct=F)

k <- table(matched_data1$Arm, matched_data1$Hematuria)
k
chisq.test(matched_data1$Arm, matched_data1$Hematuria, correct=F)

l <- table(matched_data1$Arm, matched_data1$CVA)
l
chisq.test(matched_data1$Arm, matched_data1$CVA, correct=F)

m <- table(matched_data1$Arm, matched_data1$Suprapubic)
m
chisq.test(matched_data1$Arm, matched_data1$Suprapubic, correct=F)

n <- table(matched_data1$Arm, matched_data1$Hypotension)
n
chisq.test(matched_data1$Arm, matched_data1$Hypotension, correct=F)

o <- table(matched_data1$Arm, matched_data1$ICU)
o
chisq.test(matched_data1$Arm, matched_data1$ICU, correct=F)

wilcox.test(matched_data1$CCI ~ matched_data1$Arm)
matched_data1%>%group_by(Arm)%>%summarise(Median=median(CCI), IQR_25 = quantile(CCI, probs = 0.25), IQR_75 = quantile(CCI, probs = 0.75))

wilcox.test(matched_data1$APACHE ~ matched_data1$Arm, na.rm=T)
matched_data1%>%group_by(Arm)%>%summarise(Median=median(APACHE, na.rm=T), IQR_25 = quantile(APACHE, na.rm=T, probs = 0.25), IQR_75 = quantile(APACHE, na.rm=T, probs = 0.75))

p <- table(matched_data1$Arm, matched_data1$UncompUTI)
p
chisq.test(matched_data1$Arm, matched_data1$UncompUTI, correct=F)

q <- table(matched_data1$Arm, matched_data1$CompUTI)
q
chisq.test(matched_data1$Arm, matched_data1$CompUTI, correct=F)

r <- table(matched_data1$Arm, matched_data1$Pyelonephritis)
r
chisq.test(matched_data1$Arm, matched_data1$Pyelonephritis, correct=F)

s <- table(matched_data1$Arm, matched_data1$Catheter)
s
chisq.test(matched_data1$Arm, matched_data1$Catheter, correct=F)

t <- table(matched_data1$Arm, matched_data1$E.coli)
t
chisq.test(matched_data1$Arm, matched_data1$E.coli, correct=F)

u <- table(matched_data1$Arm, matched_data1$K.pneumoniae)
u
chisq.test(matched_data1$Arm, matched_data1$K.pneumoniae, correct=F)

v <- table(matched_data1$Arm, matched_data1$P.mirabilis)
v
chisq.test(matched_data1$Arm, matched_data1$P.mirabilis, correct=F)

# Table 3 Variables
pacman::p_load(catfun,epiR)
w <- table(matched_data1$Arm, matched_data1$ClinSuccess)
w
chisq.test(matched_data1$Arm, matched_data1$ClinSuccess, correct=F)
epi.2by2(w, conf.level = 0.95)

x <- table(matched_data1$Arm, matched_data1$SympRes_48)
x
chisq.test(matched_data1$Arm, matched_data1$SympRes_48, correct=F)
epi.2by2(x, conf.level = 0.95)

y <- table(matched_data1$Arm, matched_data1$TempRes_48)
y
chisq.test(matched_data1$Arm, matched_data1$TempRes_48, correct=F)
epi.2by2(y, conf.level = 0.95)

z <- table(matched_data1$Arm, matched_data1$AbsentUTI_6mo)
z
chisq.test(matched_data1$Arm, matched_data1$AbsentUTI_6mo, correct=F)
epi.2by2(z, conf.level = 0.95)

# Table 4 Variables
wilcox.test(matched_data1$Time_to_ClinRes ~ matched_data1$Arm, na.rm=T)
matched_data1%>%group_by(Arm)%>%summarise(Median=median(Time_to_ClinRes, na.rm=T), IQR_25 = quantile(Time_to_ClinRes, na.rm=T, probs = 0.25), IQR_75 = quantile(Time_to_ClinRes, na.rm=T, probs = 0.75))

wilcox.test(matched_data1$HospLOS ~ matched_data1$Arm)
matched_data1%>%group_by(Arm)%>%summarise(Median=median(HospLOS), IQR_25 = quantile(HospLOS, probs = 0.25), IQR_75 = quantile(HospLOS, probs = 0.75))

wilcox.test(matched_data1$ICULOS ~ matched_data1$Arm, na.rm=T)
matched_data1%>%group_by(Arm)%>%summarise(Median=median(ICULOS, na.rm=T), IQR_25 = quantile(ICULOS, na.rm=T, probs = 0.25), IQR_75 = quantile(ICULOS, na.rm=T, probs = 0.75))

table(matched_data1$Arm, matched_data1$Mortality)
aa <- fisher.test(matched_data1$Arm, matched_data1$Mortality)
aa$p.value

table(matched_data1$Arm, matched_data1$Mortality_30)
bb <- fisher.test(matched_data1$Arm, matched_data1$Mortality_30)
bb$p.value


# Run logistic regression model with clin success as outcome, and arm as only* predictor
res.clinsuccess <- glm(ClinSuccess ~ Arm + Age + Sex + ICU + CCI + UncompUTI + CompUTI + Pyelonephritis + Catheter + Pathogens + Bacteremia + Symptoms,
           data = matched_data,
           weights = matched_data$weights, 
           family = binomial)
summary(res.clinsuccess) # Coef for arm is -0.09, NON SIG, so no difference in effect of treatment on clinical success 
# Test the coefficient using cluster robust standard error
clinsuccess.coeftest <- coeftest(res.clinsuccess, vcov. = vcovCL, cluster = ~subclass)
# Calculate the CI based on cluster robust std error
clinsuccess.coefci <- coefci(res.clinsuccess, vcov. = vcovCL, cluster = ~subclass, level = 0.95)
# results indicate clin success among piptazo is 9% (0.09 points) lower than among carbapenem (95% CI [-0.66, 0.48])
write.csv(clinsuccess.coeftest, "./Tables/clinsuccess.coeftest.csv")
write.csv(clinsuccess.coefci, "./Tables/clinsuccess.coefci.csv")

table(data$Arm,data$TempRes_48)
res.SympRes_48 <- glm(SympRes_48 ~ Arm + Age + Sex + ICU + CCI + UncompUTI + CompUTI + Pyelonephritis + Catheter + Pathogens + Bacteremia + Symptoms,
                       data = matched_data,
                       weights = matched_data$weights, 
                       family = binomial)
summary(res.SympRes_48) # Coef for arm is -0.41, NON SIG, so no difference in effect of treatment on SympRes48 
# Test the coefficient using cluster robust standard error
SympRes_48.coeftest <- coeftest(res.SympRes_48, vcov. = vcovCL, cluster = ~subclass)
# Calculate the CI based on cluster robust std error
SympRes_48.coefci <- coefci(res.SympRes_48, vcov. = vcovCL, cluster = ~subclass, level = 0.95)
write.csv(SympRes_48.coeftest, "./Tables/SympRes_48.coeftest.csv")
write.csv(SympRes_48.coefci, "./Tables/SympRes_48.coefci.csv")


res.TempRes_48 <- glm(TempRes_48 ~ Arm + Age + Sex + ICU + CCI + UncompUTI + CompUTI + Pyelonephritis + Catheter + Pathogens + Bacteremia + Symptoms,
                      data = matched_data,
                      weights = matched_data$weights, 
                      family = binomial)
summary(res.TempRes_48)# Coef for arm is 1.43, IT IS SIG, so there is a difference in effect of treatment on TempRes48
# Test the coefficient using cluster robust standard error
TempRes_48.coeftest <- coeftest(res.TempRes_48, vcov. = vcovCL, cluster = ~subclass)
# Calculate the CI based on cluster robust std error
TempRes_48.coefci <- coefci(res.TempRes_48, vcov. = vcovCL, cluster = ~subclass, level = 0.95)
write.csv(TempRes_48.coeftest, "./Tables/TempRes_48.coeftest.csv")
write.csv(TempRes_48.coefci, "./Tables/TempRes_48.coefci.csv")


res.AbsentUTI_6mo <- glm(AbsentUTI_6mo ~ Arm + Age + Sex + ICU + CCI + UncompUTI + CompUTI + Pyelonephritis + Catheter + Pathogens + Bacteremia + Symptoms,
                      data = matched_data,
                      weights = matched_data$weights, 
                      family = binomial)
summary(res.AbsentUTI_6mo)# Coef for arm is -0.07, NON SIG, so no difference in effect of treatment on AbsentUTI_6mo 
# Test the coefficient using cluster robust standard error
AbsentUTI_6mo.coeftest <- coeftest(res.AbsentUTI_6mo, vcov. = vcovCL, cluster = ~subclass)
# Calculate the CI based on cluster robust std error
AbsentUTI_6mo.coefci <- coefci(res.AbsentUTI_6mo, vcov. = vcovCL, cluster = ~subclass, level = 0.95)
write.csv(AbsentUTI_6mo.coeftest, "./Tables/AbsentUTI_6mo.coeftest.csv")
write.csv(AbsentUTI_6mo.coefci, "./Tables/AbsentUTI_6mo.coefci.csv")


res.Mortality <- glm(Mortality ~ Arm + Age + Sex + ICU + CCI + UncompUTI + CompUTI + Pyelonephritis + Catheter + Pathogens + Bacteremia + Symptoms,
                         data = matched_data,
                         weights = matched_data$weights, 
                         family = binomial)
summary(res.Mortality)# Coef for arm is -26.78, NON SIG, so no difference in effect of treatment on Mortality 
# Test the coefficient using cluster robust standard error
Mortality.coeftest <- coeftest(res.Mortality, vcov. = vcovCL, cluster = ~subclass)
# Calculate the CI based on cluster robust std error
Mortality.coefci <- coefci(res.Mortality, vcov. = vcovCL, cluster = ~subclass, level = 0.95)
write.csv(Mortality.coeftest, "./Tables/Mortality.coeftest.csv")
write.csv(Mortality.coefci, "./Tables/Mortality.coefci.csv")


res.Mortality_30 <- glm(Mortality_30 ~ Arm + Age + Sex + ICU + CCI + UncompUTI + CompUTI + Pyelonephritis + Catheter + Pathogens + Bacteremia + Symptoms,
                     data = matched_data,
                     weights = matched_data$weights, 
                     family = binomial)
summary(res.Mortality_30)# Coef for arm is -316.13, NON SIG, so no difference in effect of treatment on Mortality_30 
# Test the coefficient using cluster robust standard error
Mortality_30.coeftest <- coeftest(res.Mortality_30, vcov. = vcovCL, cluster = ~subclass)
# Calculate the CI based on cluster robust std error
Mortality_30.coefci <- coefci(res.Mortality_30, vcov. = vcovCL, cluster = ~subclass, level = 0.95)
write.csv(Mortality_30.coeftest, "./Tables/Mortality_30.coeftest.csv")
write.csv(Mortality_30.coefci, "./Tables/Mortality_30.coefci.csv")


res.all_coef<- matrix(0, 6, 4)
colnames(res.all_coef)<-c("Estimate (Tx Arm)", "p", "2.5 CI", "97.5 CI")
rownames(res.all_coef)<-c("ClinSuccess","SympRes_48","TempRes_48", "AbsentUTI_6mo", "Mortality", "Mortality_30")
clinsuccess_CI <- read.csv("./Tables/clinsuccess.coefci.csv", header = TRUE)
sympres_CI <- read.csv("./Tables/SympRes_48.coefci.csv", header = TRUE)
tempres_CI <- read.csv("./Tables/TempRes_48.coefci.csv", header = TRUE)
absentUTI_CI <- read.csv("./Tables/AbsentUTI_6mo.coefci.csv", header = TRUE)
mortality_CI <- read.csv("./Tables/Mortality.coefci.csv", header = TRUE)
mortality30_CI <- read.csv("./Tables/Mortality_30.coefci.csv", header = TRUE)
res.all_coef[1,1]=res.clinsuccess$coefficients[2]
res.all_coef[1,2]=0.7691
res.all_coef[1,3]=clinsuccess_CI$X2.5..[2]
res.all_coef[1,4]=clinsuccess_CI$X97.5..[2]

res.all_coef[2,1]=res.SympRes_48$coefficients[2]
res.all_coef[2,2]=0.2031
res.all_coef[2,3]=sympres_CI$X2.5..[2]
res.all_coef[2,4]=sympres_CI$X97.5..[2]

res.all_coef[3,1]=res.TempRes_48$coefficients[2]
res.all_coef[3,2]=0.0399
res.all_coef[3,3]=tempres_CI$X2.5..[2]
res.all_coef[3,4]=tempres_CI$X97.5..[2]

res.all_coef[4,1]=res.AbsentUTI_6mo$coefficients[2]
res.all_coef[4,2]=0.84928
res.all_coef[4,3]=absentUTI_CI$X2.5..[2]
res.all_coef[4,4]=absentUTI_CI$X97.5..[2]

res.all_coef[5,1]=res.Mortality$coefficients[2]
res.all_coef[5,2]=0.383
res.all_coef[5,3]=mortality_CI$X2.5..[2]
res.all_coef[5,4]=mortality_CI$X97.5..[2]

res.all_coef[6,1]=res.Mortality_30$coefficients[2]
res.all_coef[6,2]=0.991
res.all_coef[6,3]=mortality30_CI$X2.5..[2]
res.all_coef[6,4]=mortality30_CI$X97.5..[2]
#Must open this file to put in likelihood ratio tests 
write.csv(res.all_coef,"./Tables/All_Coef.csv")



##### figure 2
data_fig2 <- data1%>%
  mutate(Bacteremia_Success = case_when(Bacteremia=="yes" & ClinSuccess=="no" ~ "no",
                                        Bacteremia=="yes" & ClinSuccess=="yes" ~ "yes"),
         ICU_Success = case_when(ICU=="yes" & ClinSuccess=="no" ~ "no",
                                 ICU=="yes" & ClinSuccess=="yes" ~ "yes"),
         UncompUTI_Success = case_when(UncompUTI=="yes" & ClinSuccess=="no" ~ "no",
                                       UncompUTI=="yes" & ClinSuccess=="yes" ~ "yes"),
         CompUTI_Success = case_when(CompUTI=="yes" & ClinSuccess=="no" ~ "no",
                                     CompUTI=="yes" & ClinSuccess=="yes" ~ "yes"),
         Pyelonephritis_Success = case_when(Pyelonephritis=="yes" & ClinSuccess=="no" ~ "no",
                                            Pyelonephritis=="yes" & ClinSuccess=="yes" ~ "yes"),
         Catheter_Success = case_when(Catheter=="yes" & ClinSuccess=="no" ~ "no",
                                      Catheter=="yes" & ClinSuccess=="yes" ~ "yes"))

matched_fig2 <- matched_data1%>%
  mutate(Bacteremia_Success = case_when(Bacteremia=="yes" & ClinSuccess=="no" ~ "no",
                                        Bacteremia=="yes" & ClinSuccess=="yes" ~ "yes"),
         ICU_Success = case_when(ICU=="yes" & ClinSuccess=="no" ~ "no",
                                 ICU=="yes" & ClinSuccess=="yes" ~ "yes"),
         UncompUTI_Success = case_when(UncompUTI=="yes" & ClinSuccess=="no" ~ "no",
                                       UncompUTI=="yes" & ClinSuccess=="yes" ~ "yes"),
         CompUTI_Success = case_when(CompUTI=="yes" & ClinSuccess=="no" ~ "no",
                                     CompUTI=="yes" & ClinSuccess=="yes" ~ "yes"),
         Pyelonephritis_Success = case_when(Pyelonephritis=="yes" & ClinSuccess=="no" ~ "no",
                                            Pyelonephritis=="yes" & ClinSuccess=="yes" ~ "yes"),
         Catheter_Success = case_when(Catheter=="yes" & ClinSuccess=="no" ~ "no",
                                      Catheter=="yes" & ClinSuccess=="yes" ~ "yes"))

demo_data_fig2 <- data_fig2%>%
  dplyr::select(Arm, Bacteremia_Success, ICU_Success, UncompUTI_Success, CompUTI_Success, Pyelonephritis_Success, Catheter_Success, RecurrentUTI_Success)%>%
  tbl_summary(by=Arm,
              statistic=list(all_categorical()~"{n}({p}%)"),
              label = c(Bacteremia_Success = "Bacteremia", 
                        ICU_Success = "ICU Admission",
                        UncompUTI_Success = "Uncomplicated UTI", 
                        CompUTI_Success = "Complicated UTI",
                        Pyelonephritis_Success = "Pyelonephritis", 
                        Catheter_Success = "Catheter-Related", 
                        RecurrentUTI_Success = "Recurrent UTI"))%>%
  bold_labels()%>%
  italicize_levels()%>%
  add_p(pvalue_fun = ~style_pvalue(.x, digits = 2))
demo_data_fig2


demo_matched_data_fig2 <- matched_fig2%>%
  dplyr::select(Arm, Bacteremia_Success, ICU_Success, UncompUTI_Success, CompUTI_Success, Pyelonephritis_Success, Catheter_Success, RecurrentUTI_Success)%>%
  tbl_summary(by=Arm,
              statistic=list(all_categorical()~"{n}({p}%)"),
              label = c(Bacteremia_Success = "Bacteremia", 
                        ICU_Success = "ICU Admission",
                        UncompUTI_Success = "Uncomplicated UTI", 
                        CompUTI_Success = "Complicated UTI",
                        Pyelonephritis_Success = "Pyelonephritis", 
                        Catheter_Success = "Catheter-Related", 
                        RecurrentUTI_Success = "Recurrent UTI"))%>%
  bold_labels()%>%
  italicize_levels()%>%
  add_p(pvalue_fun = ~style_pvalue(.x, digits = 2))
demo_matched_data_fig2


a <- table(data_fig2$Arm, data_fig2$Bacteremia_Success)
chisq.test(data_fig2$Arm, data_fig2$Bacteremia_Success, correct=F)
aa <- epi.2by2(a, conf.level = 0.95)
b <- table(matched_fig2$Arm, matched_fig2$Bacteremia_Success)
chisq.test(matched_fig2$Arm, matched_fig2$Bacteremia_Success, correct=F)
bb <- epi.2by2(b, conf.level = 0.95)

c <- table(data_fig2$Arm, data_fig2$ICU_Success)
chisq.test(data_fig2$Arm, data_fig2$ICU_Success, correct=F)
cc <- epi.2by2(c, conf.level = 0.95)
d <- table(matched_fig2$Arm, matched_fig2$ICU_Success)
chisq.test(matched_fig2$Arm, matched_fig2$ICU_Success, correct=F)
dd <- epi.2by2(d, conf.level = 0.95)

e <- table(data_fig2$Arm, data_fig2$UncompUTI_Success)
chisq.test(data_fig2$Arm, data_fig2$UncompUTI_Success, correct=F)
ee <- epi.2by2(e, conf.level = 0.95)
f <- table(matched_fig2$Arm, matched_fig2$UncompUTI_Success)
chisq.test(matched_fig2$Arm, matched_fig2$UncompUTI_Success, correct=F)
ff <- epi.2by2(f, conf.level = 0.95)

g <- table(data_fig2$Arm, data_fig2$CompUTI_Success)
chisq.test(data_fig2$Arm, data_fig2$CompUTI_Success, correct=F)
gg <- epi.2by2(g, conf.level = 0.95)
h <- table(matched_fig2$Arm, matched_fig2$CompUTI_Success)
chisq.test(matched_fig2$Arm, matched_fig2$CompUTI_Success, correct=F)
hh <- epi.2by2(h, conf.level = 0.95)

i <- table(data_fig2$Arm, data_fig2$Pyelonephritis_Success)
chisq.test(data_fig2$Arm, data_fig2$Pyelonephritis_Success, correct=F)
ii <- epi.2by2(i, conf.level = 0.95)
j <- table(matched_fig2$Arm, matched_fig2$Pyelonephritis_Success)
chisq.test(matched_fig2$Arm, matched_fig2$Pyelonephritis_Success, correct=F)
jj <- epi.2by2(j, conf.level = 0.95)

k <- table(data_fig2$Arm, data_fig2$Catheter_Success)
chisq.test(data_fig2$Arm, data_fig2$Catheter_Success, correct=F)
kk <- epi.2by2(k, conf.level = 0.95)
l <- table(matched_fig2$Arm, matched_fig2$Catheter_Success)
chisq.test(matched_fig2$Arm, matched_fig2$Catheter_Success, correct=F)
ll <- epi.2by2(l, conf.level = 0.95)

m <- table(data_fig2$Arm, data_fig2$RecurrentUTI_Success)
chisq.test(data_fig2$Arm, data_fig2$RecurrentUTI_Success, correct=F)
mm <- epi.2by2(m, conf.level = 0.95)
n <- table(matched_fig2$Arm, matched_fig2$RecurrentUTI_Success)
chisq.test(matched_fig2$Arm, matched_fig2$RecurrentUTI_Success, correct=F)
nn <- epi.2by2(n, conf.level = 0.95)

fig2 <- matrix(0, 7, 7)
colnames(fig2)<-c("value", "lower", "upper", "matched_value", "matched_lower", "matched_upper", "index")
rownames(fig2)<-c("Bacteremia","ICU","UncompUTI", "CompUTI", "Pyelonephritis", "Catheter", "RecurrentUTI")
fig2[1,1]=aa$massoc.summary[3,2]
fig2[1,2]=aa$massoc.summary[3,3]
fig2[1,3]=aa$massoc.summary[3,4]
fig2[1,4]=bb$massoc.summary[3,2]
fig2[1,5]=bb$massoc.summary[3,3]
fig2[1,6]=bb$massoc.summary[3,4]
fig2[1,7]=1

fig2[2,1]=cc$massoc.summary[3,2]
fig2[2,2]=cc$massoc.summary[3,3]
fig2[2,3]=cc$massoc.summary[3,4]
fig2[2,4]=dd$massoc.summary[3,2]
fig2[2,5]=dd$massoc.summary[3,3]
fig2[2,6]=dd$massoc.summary[3,4]
fig2[2,7]=2

fig2[3,1]=ee$massoc.summary[3,2]
fig2[3,2]=ee$massoc.summary[3,3]
fig2[3,3]=ee$massoc.summary[3,4]
fig2[3,4]=ff$massoc.summary[3,2]
fig2[3,5]=ff$massoc.summary[3,3]
fig2[3,6]=ff$massoc.summary[3,4]
fig2[3,7]=3

fig2[4,1]=gg$massoc.summary[3,2]
fig2[4,2]=gg$massoc.summary[3,3]
fig2[4,3]=gg$massoc.summary[3,4]
fig2[4,4]=hh$massoc.summary[3,2]
fig2[4,5]=hh$massoc.summary[3,3]
fig2[4,6]=hh$massoc.summary[3,4]
fig2[4,7]=4

fig2[5,1]=ii$massoc.summary[3,2]
fig2[5,2]=ii$massoc.summary[3,3]
fig2[5,3]=ii$massoc.summary[3,4]
fig2[5,4]=jj$massoc.summary[3,2]
fig2[5,5]=jj$massoc.summary[3,3]
fig2[5,6]=jj$massoc.summary[3,4]
fig2[5,7]=5

fig2[6,1]=kk$massoc.summary[3,2]
fig2[6,2]=kk$massoc.summary[3,3]
fig2[6,3]=kk$massoc.summary[3,4]
fig2[6,4]=ll$massoc.summary[3,2]
fig2[6,5]=ll$massoc.summary[3,3]
fig2[6,6]=ll$massoc.summary[3,4]
fig2[6,7]=6

fig2[7,1]=mm$massoc.summary[3,2]
fig2[7,2]=mm$massoc.summary[3,3]
fig2[7,3]=mm$massoc.summary[3,4]
fig2[7,4]=nn$massoc.summary[3,2]
fig2[7,5]=nn$massoc.summary[3,3]
fig2[7,6]=nn$massoc.summary[3,4]
fig2[7,7]=7

fig2 <- as.data.frame(fig2)

dev.off()
annotation <- data.frame(
  x = c(-25, 20),
  y = c(-0.0001, -0.0001),
  label = c("Favors carbapenem", "Favors TZP")
)
data_plot <- ggplot(fig2, aes(y = index, x = value)) +
  geom_point(shape = 15, size = 3) +  
  geom_linerange(aes(xmin=lower, xmax=upper)) +
  geom_vline(xintercept = 0, linetype = "dashed", cex = 0.25, alpha = 0.5) +
  geom_text(data=annotation, aes(x=x, y=y, label=label)) +
  scale_y_continuous(name = "", breaks=NULL, trans = "reverse") +
  scale_x_continuous(name = "Between group risk difference (95 CI), %", breaks=c(-40, -20, 0, 20, 40, 60, 80)) +
  theme_classic() +
  theme(axis.line.y = element_blank(),
        axis.text.x = element_text(size = (10), colour = "black")) 
data_plot 
ggsave(data_plot, filename = "full.tiff", dpi = 900, type = "cairo",width = 6.5, height = 6, units = "in",  antialias="none")

matched_data_plot <- ggplot(fig2, aes(y = index, x = matched_value)) +
  geom_point(shape = 15, size = 3) +  
  geom_linerange(aes(xmin=matched_lower, xmax=matched_upper)) +
  geom_vline(xintercept = 0, linetype = "dashed", cex = 0.25, alpha = 0.5) +
  geom_text(data=annotation, aes(x=x, y=y, label=label)) +
  scale_y_continuous(name = "", breaks=NULL, trans = "reverse") +
  scale_x_continuous(name = "Between group risk difference (95 CI), %", breaks=c(-40, -20, 0, 20, 40, 60, 80)) +
  theme_classic() +
  theme(axis.line.y = element_blank(),
        axis.text.x = element_text(size = (10), colour = "black"))
matched_data_plot 
ggsave(matched_data_plot, filename = "matched.tiff", dpi = 900, type = "cairo",width = 6.5, height = 6, units = "in",  antialias="none")





########## Sensitivity analysis for estimate #####
pacman::p_load(rbounds)
m.pairs <- cbind(data[row.names(opt.match$match.matrix), "ClinSuccess"], 
                 data[opt.match$match.matrix, "ClinSuccess"])%>%na.omit()
psens(x=m.pairs[,1],y=m.pairs[,2],Gamma=6,GammaInc=1)
m.pairs <- cbind(data[row.names(opt.match$match.matrix), "HospLOS"], 
                 data[opt.match$match.matrix, "HospLOS"])%>%na.omit()
psens(x=m.pairs[,1],y=m.pairs[,2],Gamma=6,GammaInc=1)
m.pairs <- cbind(data[row.names(opt.match$match.matrix), "ICULOS"], 
                 data[opt.match$match.matrix, "ICULOS"])%>%na.omit()
psens(x=m.pairs[,1],y=m.pairs[,2],Gamma=6,GammaInc=1)
m.pairs <- cbind(data[row.names(opt.match$match.matrix), "Time_to_ClinRes"], 
                 data[opt.match$match.matrix, "Time_to_ClinRes"])%>%na.omit()
psens(x=m.pairs[,1],y=m.pairs[,2],Gamma=6,GammaInc=1)
### WhYYYYY IS THIS PACKAGE NOT WORKING ASDLFH AETH
m.pairs <- cbind(data[row.names(opt.match$match.matrix), "SympRes_48"], 
                 data[opt.match$match.matrix, "SympRes_48"])%>%na.omit()
x<-sum((m.pairs[,1]==FALSE)&(m.pairs[,2]==TRUE))
y<-sum((m.pairs[,1]==TRUE)&(m.pairs[,2]==FALSE))
binarysens(x=x,y=y,Gamma = 15, GammaInc = 2)
m.pairs <- cbind(data[row.names(opt.match$match.matrix), "TempRes_48"], 
                 data[opt.match$match.matrix, "TempRes_48"])%>%na.omit()
x<-sum((m.pairs[,1]==FALSE)&(m.pairs[,2]==TRUE))
y<-sum((m.pairs[,1]==TRUE)&(m.pairs[,2]==FALSE))
binarysens(x=x,y=y,Gamma = 15, GammaInc = 2)
m.pairs <- cbind(data[row.names(opt.match$match.matrix), "AbsentUTI_6mo"], 
                 data[opt.match$match.matrix, "AbsentUTI_6mo"])%>%na.omit()
x<-sum((m.pairs[,1]==FALSE)&(m.pairs[,2]==TRUE))
y<-sum((m.pairs[,1]==TRUE)&(m.pairs[,2]==FALSE))
binarysens(x=x,y=y,Gamma = 15, GammaInc = 2)
m.pairs <- cbind(data[row.names(opt.match$match.matrix), "Mortality"], 
                 data[opt.match$match.matrix, "Mortality"])%>%na.omit()
x<-sum((m.pairs[,1]==FALSE)&(m.pairs[,2]==TRUE))
y<-sum((m.pairs[,1]==TRUE)&(m.pairs[,2]==FALSE))
binarysens(x=x,y=y,Gamma = 15, GammaInc = 2)
m.pairs <- cbind(data[row.names(opt.match$match.matrix), "Mortality_30"], 
                 data[opt.match$match.matrix, "Mortality_30"])%>%na.omit()
x<-sum((m.pairs[,1]==FALSE)&(m.pairs[,2]==TRUE))
y<-sum((m.pairs[,1]==TRUE)&(m.pairs[,2]==FALSE))
binarysens(x=x,y=y,Gamma = 15, GammaInc = 2)

