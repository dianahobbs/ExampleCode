# Title: ADRC - Make merged brain image for AAIC abstract 
# Author: Diana Hobbs
# Date: April 2022

#collapse on CVD 
source("2.1_Rubin_Replication_Analyses_ANOVA_MR.R")
source("2.1_Rubin_Replication_Analyses_ANOVA_Tau.R")
source("2.1_Rubin_Replication_Analyses_ANOVA_PIB.R")

merged_brain=plot_grid(Colsp_MRUnthresh_cvd_row, Colsp_TauUnthresh_cvd_row, Colsp_PiBUnthresh_cvd_row, nrow=1)+theme(plot.background = element_rect(fill = "black"), panel.border = element_blank())
tiff('./Rabin_Replication_Output/CVD_GGSEG_Unthreshold.tiff', units="in", width=7, height=9, res=200,type="cairo")
print(merged_brain)
dev.off()

merged_brain.sig=plot_grid(Colsp_MRThresh_cvd_row, Colsp_TauThresh_cvd_row, Colsp_PiBThresh_cvd_row, nrow=1)+theme(plot.background = element_rect(fill = "black"), panel.border = element_blank())
tiff('./Rabin_Replication_Output/CVD_GGSEG_Thresholded.tiff', units="in", width=7, height=9, res=200,type="cairo")
print(merged_brain.sig)
dev.off()

merged_brain_SubCort=plot_grid(Colsp_MRUnthresh_cvd_row_SubCort, Colsp_TauUnthresh_cvd_row_SubCort, Colsp_PiBUnthresh_cvd_row_SubCort, nrow=1)+theme(plot.background = element_rect(fill = "black"), panel.border = element_blank())
tiff('./Rabin_Replication_Output/CVD_SubCort_GGSEG_Unthreshold.tiff', units="in", width=7, height=9, res=200,type="cairo")
print(merged_brain_SubCort)
dev.off()

merged_brain_SubCort.sig=plot_grid(Colsp_MRThresh_cvd_row_SubCort, Colsp_TauThresh_cvd_row_SubCort, Colsp_PiBThresh_cvd_row_SubCort, nrow=1)+theme(plot.background = element_rect(fill = "black"), panel.border = element_blank())
tiff('./Rabin_Replication_Output/CVD_SubCort_GGSEG_Thresholded.tiff', units="in", width=7, height=9, res=200,type="cairo")
print(merged_brain_SubCort.sig)
dev.off()


Test=plot_grid(merged_brain, merged_brain.sig, ncol=1)+labs(title = "CVD Risk")+theme(plot.background = element_rect(fill = "black"), panel.border = element_blank())
tiff('./Rabin_Replication_Output/CVD_GGSEG_MergedMRTau.tiff', units="in", width=7, height=2, res=200,type="cairo")
print(Test)
dev.off()

merged_brain_vert=plot_grid(Colsp_MR_row, Colsp_Tau_row,Colsp_PiB_row, ncol = 1)
tiff('./Rabin_Replication_Output/CVD_GGSEG_MergedMRTau_Vertical.tiff', units="in", width=5, height=3, res=200,type="cairo")
print(merged_brain_vert)
dev.off()

Test_SubCort=plot_grid(merged_brain_SubCort, merged_brain_SubCort.sig, ncol=1)+labs(title = "CVD Risk")+theme(plot.background = element_rect(fill = "black"), panel.border = element_blank())
tiff('./Rabin_Replication_Output/CVD_SubCort_GGSEG_MergedMRTau.tiff', units="in", width=7, height=2, res=200,type="cairo")
print(Test_SubCort)
dev.off()

merged_brain_vert_SubCort=plot_grid(Colsp_MR_row_SubCort, Colsp_Tau_row_SubCort,Colsp_PiB_row_SubCort, ncol = 1)
tiff('./Rabin_Replication_Output/CVD_SubCort_GGSEG_MergedMRTau_Vertical.tiff', units="in", width=5, height=3, res=200,type="cairo")
print(merged_brain_vert_SubCort)
dev.off()

