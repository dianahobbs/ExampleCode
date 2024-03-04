#### Title: AAIC 24 Abstract (replication) Brian's Tau Age/Sex using PVC Tau and ADNI Amyloid
### Author: Diana Hobbs
## Date: January 2024

rm(list = ls())
gc()

pacman::p_load(Hmisc, tidyverse, dplyr, ADNIMERGE, zoo, lubridate)
# install.packages("/Users/hobbsd/Downloads/ADNIMERGE_0.0.1.tar.gz", repos = NULL, type = "source")
# help(package = "ADNIMERGE")
# names(adnimerge)

# creating a centering with 'scale()' function to use later
center_scale <- function(x) {
  scale(x, scale = FALSE)
}

# Upload data from ADNI MERGE files 
data(adnimerge,ptdemog,cdr,mmse,dxsum,ucberkeleyav45,ucberkeleyav1451_pvc)

# Determine which variables I want to merge on 
# vars = c("RID","VISCODE")

## Demographics
# Select relevant demographic data 
demo <- ptdemog%>%mutate(Female=case_when(PTGENDER=="Female" ~ 1, PTGENDER=="Male" ~ 0),
                         sex=str_replace_all(PTGENDER, c("Male"="male","Female"="female")),
                         race.binary=case_when(PTRACCAT=="White" ~ 0, PTRACCAT!="White" ~ 1),
                         pseudo.birth.date=USERDATE-as.numeric(AGE*365.25))%>%
  rename(race=PTRACCAT,age=AGE)%>%group_by(RID)%>%
  fill(c(sex, Female, race, race.binary, pseudo.birth.date), .direction = 'down')%>%
  arrange(USERDATE)%>%filter(row_number()==1)%>%ungroup()%>%mutate(centered.age=center_scale(age))%>%
  dplyr::select(RID,sex,Female,race,race.binary,pseudo.birth.date,age)

## Clinical 
# Select relevant Dx data 
dx <- dxsum%>%rename(dx1=DIAGNOSIS,dx.date=USERDATE)%>%dplyr::select(RID,dx1,dx.date)

# Merge CDR and MMSE data and select relevant variables
clin <- merge(cdr, mmse, by="RID")%>%rename(cdr.date=USERDATE.x,mmse.date=USERDATE.y,cdr=CDGLOBAL,mmse=MMSCORE)%>%
  mutate(cdr.mmse.lag=abs(as.numeric(difftime(cdr.date,mmse.date,unit="weeks"))/52.25))%>%filter(between(cdr.mmse.lag,0,1.5))%>%
  dplyr::select(RID,cdr.date, mmse.date,cdr,mmse,cdr.mmse.lag)

# merge dx and clin
clin.dx <- merge(dx, clin, by="RID")%>%mutate(cdr.dx.lag=abs(as.numeric(difftime(cdr.date,dx.date,unit="weeks"))/52.25))%>%
  filter(between(cdr.dx.lag,0,1.5))

# Merge clin with demographics and dx 
demo.clin = merge(demo,clin.dx,by="RID")

# Upload amyloid information 
av45 <- ucberkeleyav45%>%mutate(centiloids=(SUMMARYSUVR_WHOLECEREBNORM*188.22)-189.16)%>%
  rename(av45.date=EXAMDATE,av45.sumsuvr=SUMMARYSUVR_WHOLECEREBNORM,av45.ampos1.11=SUMMARYSUVR_WHOLECEREBNORM_1.11CUTOFF)%>%
  filter(av45.ampos1.11==0)%>%dplyr::select(RID,av45.date,av45.sumsuvr,av45.ampos1.11,centiloids)

tau <- read.csv("/Users/hobbsd/Documents/GitHub/Brian_Requests/forBrian/SexAgeTau/ADNI_Rep/UCBERKELEY_TAUPVC_6MM_11Dec2023.csv")%>%
  rename(tau.date=SCANDATE,
         CHORPLEX=CHOROID_PLEXUS_SUVR,
         SSTBANK = CTX_BANKSSTS_SUVR,
         CAUDANTCNG = CTX_CAUDALANTERIORCINGULATE_SUVR,
         CAUDMIDFRN = CTX_CAUDALMIDDLEFRONTAL_SUVR,
         CUNEUS = CTX_CUNEUS_SUVR,
         ENTORHINAL  = CTX_ENTORHINAL_SUVR,
         FRNPOLE = CTX_FRONTALPOLE_SUVR,
         FUSIFORM = CTX_FUSIFORM_SUVR,
         INFERPRTL = CTX_INFERIORPARIETAL_SUVR,
         INFERTMP = CTX_INFERIORTEMPORAL_SUVR,
         INSULA = CTX_INSULA_SUVR,
         ISTHMUSCNG = CTX_ISTHMUSCINGULATE_SUVR,
         LATOCC = CTX_LATERALOCCIPITAL_SUVR,
         LATORBFRN = CTX_LATERALORBITOFRONTAL_SUVR,
         LINGUAL = CTX_LINGUAL_SUVR,
         MEDORBFRN = CTX_MEDIALORBITOFRONTAL_SUVR,
         MIDTMP = CTX_MIDDLETEMPORAL_SUVR,
         PARACNTRL = CTX_PARACENTRAL_SUVR,
         PARAHPCMPL = CTX_PARAHIPPOCAMPAL_SUVR,
         PARSOPCLRS = CTX_PARSOPERCULARIS_SUVR,
         PARSORBLS = CTX_PARSORBITALIS_SUVR,
         PARSTRNGLS = CTX_PARSTRIANGULARIS_SUVR,
         PERICLCRN = CTX_PERICALCARINE_SUVR,
         POSTCNTRL = CTX_POSTCENTRAL_SUVR,
         POSTCNG = CTX_POSTERIORCINGULATE_SUVR,
         PRECNTRL  = CTX_PRECENTRAL_SUVR,
         PRECUNEUS = CTX_PRECUNEUS_SUVR,
         ROSANTCNG = CTX_ROSTRALANTERIORCINGULATE_SUVR,
         ROSMIDFRN = CTX_ROSTRALMIDDLEFRONTAL_SUVR,
         SUPERFRN = CTX_SUPERIORFRONTAL_SUVR,
         SUPERPRTL = CTX_SUPERIORPARIETAL_SUVR,
         SUPERTMP = CTX_SUPERIORTEMPORAL_SUVR,
         SUPRAMRGNL = CTX_SUPRAMARGINAL_SUVR,
         TMPPOLE = CTX_TEMPORALPOLE_SUVR,
         TRANSTMP = CTX_TRANSVERSETEMPORAL_SUVR,
         AMYGDALA = AMYGDALA_SUVR,
         CAUD = CAUDATE_SUVR,
         HIPPOCAMPUS = HIPPOCAMPUS_SUVR,
         PALLIDUM = PALLIDUM_SUVR,
         PUTAMEN = PUTAMEN_SUVR,
         THALAMUS_PROPER = THALAMUS_PROPER_SUVR,
         VENTRALDC = VENTRALDC_SUVR)%>%
  dplyr::select(RID,tau.date,CHORPLEX,SSTBANK,CAUDANTCNG,CAUDMIDFRN,CUNEUS,ENTORHINAL,FRNPOLE,FUSIFORM,INFERPRTL,INFERTMP,INSULA,ISTHMUSCNG,LATOCC,LATORBFRN,LINGUAL,
                MEDORBFRN,MIDTMP,PARACNTRL,PARAHPCMPL,PARSOPCLRS,PARSORBLS,PARSTRNGLS,PERICLCRN,POSTCNTRL,POSTCNG,PRECNTRL,PRECUNEUS,ROSANTCNG,ROSMIDFRN,SUPERFRN,SUPERPRTL,SUPERTMP,
                SUPRAMRGNL,TMPPOLE,TRANSTMP,AMYGDALA,CAUD,HIPPOCAMPUS,PALLIDUM,PUTAMEN,THALAMUS_PROPER,VENTRALDC)

image <- merge(av45, tau, by="RID")%>%mutate(amy.tau.lag = abs(as.numeric(difftime(av45.date,tau.date,unit="weeks"))/52.25))%>%
  filter(between(amy.tau.lag,0,1.5))%>%group_by(RID)%>%arrange(amy.tau.lag)%>%filter(row_number()==1)

data <- merge(image,demo.clin,by="RID")%>%mutate(cdr.tau.lag = abs(as.numeric(difftime(cdr.date,tau.date,unit="weeks"))/52.25))%>%
  filter(between(cdr.tau.lag,0,1.5))%>%group_by(RID)%>%arrange(cdr.tau.lag)%>%filter(row_number()==1)

############ DO BRIAN'S ANALYSES FOR REPLICATION #####################
pacman::p_load(readxl,writexl,ggplot2,scales,cowplot,ggseg,MASS,sfsmisc,tidyverse) 

#READ IN DATA
tauage = data
write.csv(tauage, "/Users/hobbsd/Documents/GitHub/Brian_Requests/forBrian/SexAgeTau/AAIC_2024/PVC_v_nonPVC_NO-DIAN/data_ADNI_PVC.csv")
ggseg_list =read.csv("/Users/hobbsd/Documents/GitHub/Brian_Requests/forBrian/SexAgeTau/ggseg_namelist.csv")


ADNI_IDs <- read.csv("/Users/hobbsd/Documents/GitHub/Brian_Requests/forBrian/SexAgeTau/AAIC_2024/PVC_v_nonPVC_NO-DIAN/ADNI_IDs.csv")%>%dplyr::select(RID)
tauage <- inner_join(tauage, ADNI_IDs, by="RID")

mean(tauage$age)
sd(tauage$age)
range(tauage$age)

# reassign variable types
tauage$sex=as.factor(tauage$sex); tauage$Female=as.factor(tauage$Female); tauage$race.binary=as.factor(tauage$race.binary); 
tauage$sex <- relevel(tauage$sex, "male") # set reference as we expect women to have greater levels

# Create empty matrices to hold results 
#empty_matrix <- matrix(0, 42, 5)
empty_matrix <- data.frame(region = character(), B = numeric(), t = numeric(), p = numeric(), stringsAsFactors = FALSE)
model_names_reg <- c("model_age","model_sex","model_race","model_ab")
model_names_rob <- c("model_age_rob","model_sex_rob","model_race_rob","model_ab_rob")
model_names2_reg <- c("model_age2","model_sex2","model_race2","model_ab2","model_int")
model_names2_rob <- c("model_age2_rob","model_sex2_rob","model_race2_rob","model_ab2_rob", "model_int_rob")
models <- c(model_names_reg, model_names_rob, model_names2_reg, model_names2_rob)
# results_list <- lapply(models, function(model_name) empty_matrix)
result_matrices <- list()
for (model_name in models) {
  result_matrices[[model_name]] <- empty_matrix
}

# Function to create dataframes for lm and rlm
create_model_df <- function(region_name, model_fit, covariate_index, robust = FALSE) {
  if (!robust) {
    p_value_col <- summary(model_fit)$coefficients[covariate_index, 4]
  } else {
    p_value_col <- (f.robftest(model_fit, var = covariate_index))$p.value
  }
  
  data.frame(
    region = region_name,
    B = coef(model_fit)[covariate_index],
    t = summary(model_fit)$coefficients[covariate_index, 3],
    p = p_value_col
  )
}

# Loop through brain region columns and save values to matrices 
start_col <- which(colnames(tauage)=="CHORPLEX") 
end_col <- which(colnames(tauage)=="VENTRALDC") 
nreg <- end_col-start_col+1 #index loop
outcount <- 1 #start a count for writing out to matrices

# Loop through regions
for (i in start_col:end_col) {
  region_name <- colnames(tauage)[i]
  
  # Run lm and rlm for each region
  region_fit <- lm(get(region_name) ~ sex + centiloids + race.binary + age, data = tauage)
  robust_fit <- rlm(get(region_name) ~ sex + centiloids + race.binary + age, data = tauage)
  region_fit_int <- lm(get(region_name) ~ sex + centiloids + race.binary + age + age * sex, data = tauage)
  robust_fit_int <- rlm(get(region_name) ~ sex + centiloids + race.binary + age + age * sex, data = tauage)
  
  for (model_name in models) {
    covariate_index <- switch(model_name,
                              model_sex = 2,model_sex_rob = 2,model_sex2 = 2,model_sex2_rob = 2,
                              model_ab = 3,model_ab_rob = 3,model_ab2 = 3,model_ab2_rob = 3,
                              model_race = 4,model_race_rob = 4,model_race2 = 4,model_race2_rob = 4,
                              model_age = 5, model_age_rob = 5,model_age2 = 5,model_age2_rob = 5,
                              model_int = 6,model_int_rob = 6)
    
    if (model_name %in% model_names_reg) {
      result_matrices[[model_name]][outcount, ] <- create_model_df(region_name, region_fit, covariate_index)
    } else if (model_name %in% model_names2_reg) {
      result_matrices[[model_name]][outcount, ] <- create_model_df(region_name, region_fit_int, covariate_index)
    } else if (model_name %in% model_names_rob) {
      result_matrices[[model_name]][outcount, ] <- create_model_df(region_name, robust_fit, covariate_index, robust = TRUE)
    } else if (model_name %in% model_names2_rob) {
      result_matrices[[model_name]][outcount, ] <- create_model_df(region_name, robust_fit_int, covariate_index, robust = TRUE)
    }
  }
  
  outcount <- outcount + 1
}    

# Do BH corrections 
for (model_name in models) {
  result_matrices[[model_name]][, "p_corr"] <- as.numeric(p.adjust(result_matrices[[model_name]][, 4], method = "BH", n = length(result_matrices[[model_name]][, 4])))
}

# extract model data into new, tight matrices. Note that Choroid plexus isn't in this as it can't be plotted.
cort1=agrep("SSTBANK", result_matrices[["model_age"]][["region"]], ignore.case=TRUE, value=FALSE, max.distance=0.1 ) #gets position of first
cort2=agrep("TRANSTMP", result_matrices[["model_age"]][["region"]], ignore.case=TRUE, value=FALSE, max.distance=0.1 ) #gets position of last
sub1=agrep("AMYGDALA", result_matrices[["model_age"]][["region"]], ignore.case=TRUE, value=FALSE, max.distance=0.1 ) #gets position of first subcortical
sub2=agrep("VENTRALDC", result_matrices[["model_age"]][["region"]], ignore.case=TRUE, value=FALSE, max.distance=0.1 ) #gets position of last subcortical


extract_and_return <- function(model_matrix, cort1, cort2, sub1, sub2, model_type, robust = FALSE) {
  cort_name <- paste0(model_type, "_cort")
  sub_name <- paste0(model_type, "_sub")
  
  if (robust) {
    cort_name <- paste0(cort_name, "_rob")
    sub_name <- paste0(sub_name, "_rob")
  }
  
  cort_df <- model_matrix[cort1:cort2, ]
  sub_df <- model_matrix[sub1:sub2, ]
  
  # Return the data frames
  return(list(cort_name = cort_df, sub_name = sub_df))
}

# Initialize lists
cort_list <- list()
sub_list <- list()

# List of model names
next_models <- c("model_age","model_sex","model_race","model_ab","model_int", 
                 "model_age_rob","model_sex_rob",'model_race_rob',"model_ab_rob","model_int_rob")

for (model_name in next_models) {
  cort_sub_data <- extract_and_return(result_matrices[[model_name]], cort1, cort2, sub1, sub2, model_name)
  
  # Append the data frames to the lists
  cort_list[[model_name]] <- cort_sub_data$cort_name
  sub_list[[model_name]] <- cort_sub_data$sub_name
}

for (i in seq_along(cort_list)) {
  cort_list[[i]][, 1] <- ggseg_list$cortical
  cort_list[[i]][, 2:5] <- sapply(cort_list[[i]][, 2:5], as.numeric)
}

for (i in seq_along(sub_list)) {
  sub_list[[i]][, 1] <- ggseg_list[1:7, 2]
  sub_list[[i]][, 2:5] <- sapply(sub_list[[i]][, 2:5], as.numeric)
}

path <- "/Users/hobbsd/Documents/GitHub/Brian_Requests/forBrian/3D_Brain/deployment_models/ADNI_PVC_models/"
for (model_name in names(cort_list)) {
  model_data <- cort_list[[model_name]]
  csv_filename <- paste0("cort_", model_name, ".csv")
  write.csv(model_data, file = file.path(path, csv_filename), row.names = FALSE)
}

for (model_name in names(sub_list)) {
  model_data <- sub_list[[model_name]]
  csv_filename <- paste0("subcort_", model_name, ".csv")
  write.csv(model_data, file = file.path(path, csv_filename), row.names = FALSE)
}


# Initialize a list to store the significant matrices
cort_list_sig <- cort_list
sub_list_sig <- sub_list

# Define a function to subset significant values
subset_significant <- function(matrix_list) {
  # Loop through each matrix in the list
  for (i in seq_along(matrix_list)) {
    # Extract the current matrix
    current_matrix <- matrix_list[[i]]
    # Subset only the significant values
    significant_values <- current_matrix[current_matrix[, "p_corr"] < 0.05, ]
    # Save the significant values back into the list
    matrix_list[[i]] <- significant_values
  }
  return(matrix_list)
}

# Apply the function to cort_list and sub_list
cort_list_sig <- subset_significant(cort_list)
sub_list_sig <- subset_significant(sub_list)

##### Make age images #####
#This theme specifies parameters for all of the plots. Change here to change any of them
custom_theme = list(
  scale_fill_gradientn(colours=c("dodgerblue4","light blue","white","firebrick","goldenrod"), limits=c(-5,5), oob=squish),
  theme(text = element_text(colour = "white")),  theme(axis.text = element_blank(), axis.title = element_blank()), 
  theme(plot.background = element_rect(fill = "black",colour = NA)), theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), theme(panel.border = element_blank()), guides(fill="none", color="none"))

custom_theme_sub = list(
  scale_fill_gradientn(colours=c("dodgerblue4","light blue","white","firebrick","goldenrod"), limits=c(-5,5), oob=squish),
  theme(text = element_text(colour = "white")),  theme(axis.text = element_blank(), axis.title = element_blank()), 
  theme(plot.background = element_rect(fill = "black",colour = NA)), theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), theme(panel.border = element_blank()))

create_ggseg_plot <- function(model_list, model_name, is_cortical = TRUE) {
  title_ <- sub("model_", "", model_name)
  
  if (is_cortical) {
    gg_plot <- ggseg(.data = model_list[[model_name]], mapping = aes_string(fill = "t"),
                     position = "stacked", colour = "black", size = 0.5, hemisphere = "left") + 
      labs(title = title_) + custom_theme
  } else if (is_cortical == FALSE) {
    gg_plot <- ggseg(.data = model_list[[model_name]], mapping = aes_string(fill = "t"),
                     position = "stacked", colour = "black", size = 0.5, hemisphere = "left", atlas=aseg) + 
      labs(title = title_) + custom_theme_sub
  }
  
  return(gg_plot)
}

selected_mods = c(model_names_reg, "model_int")
selected_mod_list <- cort_list[selected_mods]

cort_plots <- list()
sub_plots <- list()

for (model_name in names(selected_mod_list)) {
  title_ <- sub("model_", "", model_name)  # Remove "model_" from the title
  
  # Create variable names
  fig_t_name <- paste0(title_, "_fig_t")
  sig_t_name <- paste0(title_, "_sig_t")
  fig_sub_t_name <- paste0(title_, "_fig_sub_t")
  sig_sub_t_name <- paste0(title_, "_sig_sub_t")
  
  # Append the plots to the list
  cort_plots[[fig_t_name]] <- create_ggseg_plot(cort_list, model_name, is_cortical = TRUE)
  cort_plots[[sig_t_name]] <- create_ggseg_plot(cort_list_sig, model_name, is_cortical = TRUE)
  sub_plots[[fig_sub_t_name]] <- create_ggseg_plot(sub_list, model_name, is_cortical = FALSE)
  sub_plots[[sig_sub_t_name]] <- create_ggseg_plot(sub_list_sig, model_name, is_cortical = FALSE)
}

age_row=plot_grid(cort_plots[["age_fig_t"]], sub_plots[["age_fig_sub_t"]],  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") ) 
sex_row=plot_grid(cort_plots[["sex_fig_t"]], sub_plots[['sex_fig_sub_t']],  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
ab_row=plot_grid(cort_plots[["ab_fig_t"]], sub_plots[['ab_fig_sub_t']],  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
race_row=plot_grid(cort_plots[["race_fig_t"]], sub_plots[['race_fig_sub_t']],  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
int_row=plot_grid(cort_plots[["int_fig_t"]], sub_plots[["int_fig_sub_t"]],  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
merged_brain=plot_grid(age_row, sex_row, ab_row, race_row, int_row, ncol=1)+theme(plot.background = element_rect(fill = "black"), panel.border = element_blank())

tiff('/Users/hobbsd/Documents/GitHub/Brian_Requests/forBrian/SexAgeTau/AAIC_2024/PVC_v_nonPVC_NO-DIAN/ADNI_PVC/Supplemental_Figure1_Unthresholded.tiff', units="in", width=7, height=9, res=200,type="cairo" )
print(merged_brain)
dev.off()

age_row_sig=plot_grid(cort_plots[["age_sig_t"]], sub_plots[["age_sig_sub_t"]],  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") ) 
sex_row_sig=plot_grid(cort_plots[["sex_sig_t"]], sub_plots[["sex_sig_sub_t"]],  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
int_row_sig=plot_grid(cort_plots[["int_sig_t"]], sub_plots[["int_sig_sub_t"]],  align="hv", axis="tb", rel_widths = c(1, .3))+theme(plot.background = element_rect(fill = "black") )  
merged_brain_sig=plot_grid(age_row_sig, sex_row_sig, ncol=1)+theme(plot.background = element_rect(fill = "black"), panel.border = element_blank())

tiff('/Users/hobbsd/Documents/GitHub/Brian_Requests/forBrian/SexAgeTau/AAIC_2024/PVC_v_nonPVC_NO-DIAN/ADNI_PVC/Figure2_Thresholded.tiff', units="in", width=7, height=4, res=200,type="cairo" )
print(merged_brain_sig)
dev.off()

##### Make Age Scatterplots #####
scatter_theme = list(
  geom_point(size=2),  geom_smooth(method=lm, se=FALSE, fullrange=TRUE, aes(colour = factor(sex))),  scale_color_manual(values=c('#3182BD','#CC6666')), #gold #E69F00
  theme_bw(), theme(text = element_text(size=16)),  theme(legend.position="none"), theme(legend.title=element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()), scale_y_continuous(
    labels = scales::number_format(accuracy = 0.1))
)
# theme(legend.position="bottom"

### Most significant cortical and subcortical regions
## AGE
# Cortical: medorbfrn, suppar
(cort_list_sig$model_age[order(cort_list_sig$model_age[,5]),])[1:2,]
age_medorbfrn=ggplot(tauage, aes(x=age, y=MEDORBFRN, linetype=sex )) +
  xlab("Age") + ylab("Medial Orbitofrontal  Tau")+  scatter_theme
age_suppar=ggplot(tauage, aes(x=age, y=SUPERPRTL, linetype=sex )) +
  xlab("Age") + ylab("Superior Parietal Tau")+  scatter_theme


# Subcortical: putamen, caudate
(sub_list_sig$model_age[order(sub_list_sig$model_age[,5]),])[1:2,]

age_putamen=ggplot(tauage, aes(x=age, y=PUTAMEN, linetype=sex)) +
  xlab("Age") + ylab("Putamen Tau")+scatter_theme+ theme(axis.title.x = element_blank())
age_caudate=ggplot(tauage, aes(x=age, y=CAUD, linetype=sex )) +
  xlab("Age") + ylab("Caudate Tau")+  scatter_theme + theme(axis.title.x = element_blank())


## INTERACTION: NO significant interactions 
(cort_list_sig$model_int[order(cort_list_sig$model_int[,5]),])[1:2,]
(sub_list_sig$model_int[order(sub_list_sig$model_int[,5]),])[1:2,]

##### Merge panels into figures for the paper using cowplot
scatter_age=plot_grid(age_putamen, age_caudate, age_medorbfrn, age_suppar, ncol=2, labels = c("A)", "B)", "C)","D)"), align = "v", axis = "bt")+theme(plot.background = element_rect(fill = "white"), panel.border = element_blank())
tiff('/Users/hobbsd/Documents/GitHub/Brian_Requests/forBrian/SexAgeTau/AAIC_2024/PVC_v_nonPVC_NO-DIAN/ADNI_PVC/Scatter_age.tiff', units="in", width=8.5, height=8.5, res=200,type="cairo" )
print(scatter_age)
dev.off()


##### Make Sex Violin Plots #####
violin_theme=list(
  geom_violin(trim=FALSE, alpha=1),  geom_boxplot(width=0.1, fill="white", outlier.shape = NA),scale_fill_manual(values=c('#3182BD','#CC6666'), labels=c("Men","Women")), #gold #E69F00
  theme_classic(), theme(text = element_text(size=16)),  theme(legend.position="none"), theme(legend.title=element_blank()), scale_x_discrete(breaks=c("male","female"), labels=c("Men","Women")), theme(axis.title.x = element_blank()), scale_y_continuous(
    labels = scales::number_format(accuracy = 0.1))
)


### Most significant cortical and subcortical regions
## SEX 
# Cortical: pericalc, cuneus, lingual
(cort_list_sig$model_sex[order(cort_list_sig$model_sex[,5]),])[1:3,]

sex_pericalc=ggplot(tauage, aes(x=sex, y=PERICLCRN, fill=sex))+ ylab("Pericalcarine Tau")+  violin_theme  # If you want figure title labs(title="Rostral Middle Frontal Cortex",x="Sex", y = "Rostral Middle Frontal Tau")
sex_cuneus=ggplot(tauage, aes(x=sex, y=CUNEUS, fill=sex))+  ylab("Cuneus Tau")+  violin_theme  #xlab("Group") + ylab("Temporal Pole Tau") #If you don't want a figure title
sex_lingual=ggplot(tauage, aes(x=sex, y=LINGUAL, fill=sex))+ ylab("Lingual Tau")+  violin_theme  

# Subcortical:thalamus 
(sub_list_sig$model_sex[order(sub_list_sig$model_sex[,5]),])[1:2,]
sex_thalamus=ggplot(tauage, aes(x=sex, y=THALAMUS_PROPER, fill=sex))+ ylab("Thalamus Tau")+  violin_theme  


##### Merge panels into figures for the paper using cowplot
violin_sex=plot_grid(sex_pericalc, sex_cuneus, sex_lingual, sex_thalamus, ncol=2, labels = c("E)", "F)", "G)","H)"), align = "v", axis = "bt")+theme(plot.background = element_rect(fill = NA), panel.border = element_blank())
tiff('/Users/hobbsd/Documents/GitHub/Brian_Requests/forBrian/SexAgeTau/AAIC_2024/PVC_v_nonPVC_NO-DIAN/ADNI_PVC/Violin_sex.tiff', units="in", width=8.5, height=8.5, res=200,type="cairo" )
print(violin_sex)
dev.off()


combined_sex_age=plot_grid(age_putamen, age_caudate, NULL , sex_pericalc, sex_cuneus, age_medorbfrn, age_suppar, NULL , sex_lingual, sex_thalamus, nrow=2, labels = c("A)", "B)", "" ,"E)","F)", "C)", "D)", "" , "G)","H)"),  align = "v", axis = "bt", rel_widths = c(1.5, 1.5, .1, 1, 1))
# combined_sex_age=plot_grid(age_putamen, age_caudate, NULL , sex_rmf, sex_front_pole, nrow=1, labels = c("A)", "B)", "" ,"E)", "F)"),  align = "hv", rel_widths = c(1, 1, .2, 1, 1), axis = "bt") 
tiff('/Users/hobbsd/Documents/GitHub/Brian_Requests/forBrian/SexAgeTau/AAIC_2024/PVC_v_nonPVC_NO-DIAN/ADNI_PVC/Combined_sex_age.tiff', units="in", width=12, height=6.75, res=200,type="cairo" )
print(combined_sex_age)
dev.off()


all_data <- list()
all_data_sig <- list()
all_data_comb <- list()
for (model_name in c(names(cort_list), names(cort_list_sig))){
  if (model_name %in% names(cort_list)) {
    model_data <- cort_list[[model_name]]
    model_sub_data <- sub_list[[model_name]]
    all_data[[model_name]] <- rbind(model_data, model_sub_data)
  }
  if (model_name %in% names(cort_list_sig)) {
    model_data_sig <- cort_list_sig[[model_name]]
    model_sub_data_sig <- sub_list_sig[[model_name]]
    all_data_sig[[model_name]] <- rbind(model_data_sig, model_sub_data_sig)
  }
}

for (model_name in names(all_data)){
  if (!grepl("_rob", model_name)) {
    model_reg_name <- model_name
    model_rob_name <- paste0(model_name, "_rob")
    
    if (model_rob_name %in% names(all_data)) {
      model_reg <- all_data[[model_reg_name]]
      model_rob <- all_data[[model_rob_name]][, 2:5]
      all_data_comb[[model_reg_name]] <- cbind.data.frame(model_reg, model_rob)
    }
  }                                             
}
for (model_name in names(all_data_comb)) {
  cols19 <- c(1:9);  colnames(all_data_comb[[model_name]])[cols19] <- c("Region", "B", "t", "p", "p (cor)", "B", "t", "p", "p (adj)")
  cols26 <- c(2,6); cols37 <- c(3,7); cols4589 <- c(4,5,8,9)
  all_data_comb[[model_name]][, cols26] <- lapply(all_data_comb[[model_name]][, cols26], function(x) round(x, digits = 3))
  all_data_comb[[model_name]][, cols37] <- lapply(all_data_comb[[model_name]][, cols37], function(x) round(x, digits = 2))
  all_data_comb[[model_name]][, cols4589] <- lapply(all_data_comb[[model_name]][, cols4589], function(x) round(x, digits = 10))
}

cols15 <- c(1,5)
colnames(all_data_sig[["model_age"]])[cols15] <- c("Region", "p (adj)")
colnames(all_data_sig[["model_sex"]])[cols15] <- c("Region", "p (adj)")

all_data_sig[["model_age"]][2] <- lapply(all_data_sig[['model_age']][2], function(x) round(x, digits = 3))
all_data_sig[['model_age']][3] <- lapply(all_data_sig[['model_age']][3], function(x) round(x, digits = 2))
all_data_sig[['model_age']][c(4,5)] <- lapply(all_data_sig[['model_age']][c(4,5)], function(x) round(x, digits = 10))

all_data_sig[['model_sex']][2] <- lapply(all_data_sig[["model_sex"]][2], function(x) round(x, digits = 3))
all_data_sig[['model_sex']][3] <- lapply(all_data_sig[["model_sex"]][3], function(x) round(x, digits = 2))
all_data_sig[['model_sex']][c(4,5)] <- lapply(all_data_sig[["model_sex"]][c(4,5)], function(x) round(x, digits = 10))




pacman::p_load(gridExtra, grid, gtable)

create_and_save_table_sig <- function(data, filename) {
  g <- tableGrob(data, rows = NULL)
  g <- gtable_add_grob(g, grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)), t = 2, b = nrow(g), l = 1, r = ncol(g))
  g <- gtable_add_grob(g, grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)), t = 1, l = 1, r = ncol(g))
  
  tiff(filename, units = "in", width = 6, height = 7, res = 200, type = "cairo")
  grid.draw(g)
  dev.off()
}
create_and_save_table <- function(data, filename) {
  g <- tableGrob(data, rows = NULL)
  g <- gtable_add_grob(g, grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)), t = 2, b = nrow(g), l = 1, r = ncol(g))
  g <- gtable_add_grob(g, grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)), t = 1, l = 1, r = ncol(g))
  g <- gtable_add_grob(g, grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)), t = 1, b = nrow(g), l = 2, r = 5)
  
  tiff(filename, units = "in", width = 9, height = 12, res = 200, type = "cairo")
  grid.draw(g)
  dev.off()
}
# Just significant tables for the main paper
create_and_save_table_sig(all_data_sig[["model_age"]], '/Users/hobbsd/Documents/GitHub/Brian_Requests/forBrian/SexAgeTau/AAIC_2024/PVC_v_nonPVC_NO-DIAN/ADNI_PVC/age_table1.tiff')
create_and_save_table_sig(all_data_sig[["model_sex"]], '/Users/hobbsd/Documents/GitHub/Brian_Requests/forBrian/SexAgeTau/AAIC_2024/PVC_v_nonPVC_NO-DIAN/ADNI_PVC/sex_table2.tiff')

# Supplemental tables
create_and_save_table(all_data_comb[["model_age"]], '/Users/hobbsd/Documents/GitHub/Brian_Requests/forBrian/SexAgeTau/AAIC_2024/PVC_v_nonPVC_NO-DIAN/ADNI_PVC/Supplemental1_age.tiff')
create_and_save_table(all_data_comb[["model_sex"]], '/Users/hobbsd/Documents/GitHub/Brian_Requests/forBrian/SexAgeTau/AAIC_2024/PVC_v_nonPVC_NO-DIAN/ADNI_PVC/Supplemental2_sex.tiff')
create_and_save_table(all_data_comb[["model_ab"]], '/Users/hobbsd/Documents/GitHub/Brian_Requests/forBrian/SexAgeTau/AAIC_2024/PVC_v_nonPVC_NO-DIAN/ADNI_PVC/Supplemental3_ab.tiff')
create_and_save_table(all_data_comb[["model_race"]], '/Users/hobbsd/Documents/GitHub/Brian_Requests/forBrian/SexAgeTau/AAIC_2024/PVC_v_nonPVC_NO-DIAN/ADNI_PVC/Supplemental4_race.tiff')
create_and_save_table(all_data_comb[["model_int"]], '/Users/hobbsd/Documents/GitHub/Brian_Requests/forBrian/SexAgeTau/AAIC_2024/PVC_v_nonPVC_NO-DIAN/ADNI_PVC/Supplemental5_interaction.tiff')

age.all <- all_data[["model_age"]] %>% rename(age.t = "t") %>% dplyr::select("region", "age.t")
sex.all <- all_data[["model_sex"]] %>% rename(sex.t = "t") %>% dplyr::select("region", "sex.t")
race.all <- all_data[["model_race"]] %>% rename(race.t = "t") %>% dplyr::select("region", "race.t")
ab.all <- all_data[["model_ab"]] %>% rename(centiloids.t = "t") %>% dplyr::select("region", "centiloids.t")
int.all <- all_data[["model_int"]] %>% rename(int.t = "t") %>% dplyr::select("region", "int.t")

df1 <- merge(age.all, sex.all, by='region', sort=F)
df2 <- merge(race.all, ab.all, by='region', sort=F)
df <- merge(df1, df2, by='region', sort=F)

write.csv(age.all, "/Users/hobbsd/Documents/GitHub/Brian_Requests/forBrian/SexAgeTau/AAIC_2024/PVC_v_nonPVC_NO-DIAN/ADNI_PVC/ADNI_PVC_age.all.csv")
write.csv(sex.all, "/Users/hobbsd/Documents/GitHub/Brian_Requests/forBrian/SexAgeTau/AAIC_2024/PVC_v_nonPVC_NO-DIAN/ADNI_PVC/ADNI_PVC_sex.all.csv")
write.csv(race.all, "/Users/hobbsd/Documents/GitHub/Brian_Requests/forBrian/SexAgeTau/AAIC_2024/PVC_v_nonPVC_NO-DIAN/ADNI_PVC/ADNI_PVC_race.all.csv")
write.csv(ab.all, "/Users/hobbsd/Documents/GitHub/Brian_Requests/forBrian/SexAgeTau/AAIC_2024/PVC_v_nonPVC_NO-DIAN/ADNI_PVC/ADNI_PVC_centiloids.all.csv")
write.csv(df, "/Users/hobbsd/Documents/GitHub/Brian_Requests/forBrian/SexAgeTau/AAIC_2024/PVC_v_nonPVC_NO-DIAN/ADNI_PVC/ADNI_PVC_t.all.csv")
