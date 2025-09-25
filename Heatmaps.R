library(gplots)
library(DESeq2)
library(dplyr)

#### Brain Microplastics heatmaps ####

# sig genes from Br_res_s_MP_DMSO
#Rows = genes
#Columns = samples

#Log the data

BrMP_Norm <- counts(MP_Br_conds_dds, normalized = TRUE)

BrMP_DMSO <- intersect(rownames(Br_res_s_MP_DMSO), rownames(BrMP_Norm))
BrMP_DMSO_Norm_sig <- BrMP_Norm[BrMP_DMSO, ]

logBrMP <- log2(BrMP_DMSO_Norm_sig + 0.5)
sc_logBrMP <- t(scale(t(logBrMP)))
sc_logBrMP[sc_logBrMP > 3] <- 3
sc_logBrMP[sc_logBrMP < -3] <- -3
heatmap.2(sc_logBrMP, trace ='none', scale = 'none', col=bluered(50), labRow = rownames(sc_logBrMP), labCol = colnames(sc_logBrMP), cexRow = 0.7, cexCol = 0.7, breaks = seq(-3,3, l=51))


# sig genes from Br_res_s_MPD_DMSO

BrMPD_DMSO<- intersect(rownames(Br_res_s_MPD_DMSO), rownames(BrMP_Norm))
BrMPD_DMSO_Norm_sig <- BrMP_Norm[BrMPD_DMSO, ]

logBrMPD <- log2(BrMPD_DMSO_Norm_sig + 0.5)
sc_logBrMPD <- t(scale(t(logBrMPD)))
sc_logBrMPD[sc_logBrMPD > 3] <- 3
sc_logBrMPD[sc_logBrMPD < -3] <- -3
heatmap.2(sc_logBrMPD, trace ='none', scale = 'none', col=bluered(50), labRow = rownames(sc_logBrMPD), labCol = colnames(sc_logBrMPD), cexRow = 0.7, cexCol = 0.7, breaks = seq(-3,3, l=51))

#### Brain F1 ####
#sig genes from Br_res_s_F1

BrF1_Norm <- counts(F1_Br_conds_dds, normalized = TRUE)

BrF1 <- intersect(rownames(Br_res_s_F1), rownames(BrF1_Norm))
BrF1_Norm_sig <- BrF1_Norm[BrF1, ]

logBrF1 <- log2(BrF1_Norm_sig + 0.5)
sc_logBrF1 <- t(scale(t(logBrF1)))
sc_logBrF1[sc_logBrF1 > 3] <- 3
sc_logBrF1[sc_logBrF1 < -3] <- -3
heatmap.2(sc_logBrF1, trace ='none', scale = 'none', col=bluered(50), labRow = "", labCol = colnames(sc_logBrF1), cexRow = 0.7, cexCol = 0.7, breaks = seq(-3,3, l=51))

#### Gonad High-Low-DMSO ####
#sig genes from Go_res_s_HD

GoHLD_Norm <- counts(HLD_Go_conds_dds, normalized = TRUE)

Gonad_sample_conditions <- colData(HLD_Go_conds_dds)$conds
Gonad_keep_samples <- colnames(HLD_Go_conds_dds)[Gonad_sample_conditions %in% c("High", "DMSO")]
GoHLD_Norm_sub <- GoHLD_Norm[, Gonad_keep_samples]

GoHD <- intersect(rownames(Go_res_s_HD), rownames(GoHLD_Norm_sub))
GoHD_Norm_sig <- GoHLD_Norm_sub[GoHD, ]

logGoHD <- log2(GoHD_Norm_sig + 0.5)
sc_logGoHD <- t(scale(t(logGoHD)))
sc_logGoHD[sc_logGoHD > 3] <- 3
sc_logGoHD[sc_logGoHD < -3] <- -3
heatmap.2(sc_logGoHD, trace ='none', scale = 'none', col=bluered(50), labRow = "", labCol = colnames(sc_logGoHD), cexRow = 0.7, cexCol = 0.7, breaks = seq(-3,3, l=51))

#sig genes from Go_res_s_HL

GoHL <- intersect(rownames(Go_res_s_HL), rownames(GoHLD_Norm))
GoHL_Norm_sig <- GoHLD_Norm[GoHL, ]

logGoHL <- log2(GoHD_Norm_sig + 0.5)
sc_logGoHL <- t(scale(t(logGoHL)))
sc_logGoHL[sc_logGoHL > 3] <- 3
sc_logGoHL[sc_logGoHL < -3] <- -3
heatmap.2(sc_logGoHL, trace ='none', scale = 'none', col=bluered(50), labRow = "", labCol = colnames(sc_logGoHL), cexRow = 0.7, cexCol = 0.7, breaks = seq(-3,3, l=51))

#sig genes from Go_res_s_LD

GoLD <- intersect(rownames(Go_res_s_LD), rownames(GoHLD_Norm))
GoLD_Norm_sig <- GoHLD_Norm[GoLD, ]

logGoLD <- log2(GoLD_Norm_sig + 0.5)
sc_logGoLD <- t(scale(t(logGoLD)))
sc_logGoLD[sc_logGoLD > 3] <- 3
sc_logGoLD[sc_logGoLD < -3] <- -3
heatmap.2(sc_logGoLD, trace ='none', scale = 'none', col=bluered(50), labRow = rownames(sc_logGoLD), labCol = colnames(sc_logGoLD), cexRow = 0.7, cexCol = 0.7, breaks = seq(-3,3, l=51))

#### Gonad Microplastics ####

# sig genes from Go_res_s_MP_MPD

GoMP_Norm <- counts(MP_Go_conds_dds, normalized = TRUE)

GoMP_MPD <- intersect(rownames(Go_res_s_MP_MPD), rownames(GoMP_Norm))
GoMP_MPD_Norm_sig <- GoMP_Norm[GoMP_MPD, ]

logGoMP_MPD <- log2(GoMP_MPD_Norm_sig + 0.5)
sc_logMP_MPD <- t(scale(t(logGoMP_MPD)))
sc_logMP_MPD[sc_logMP_MPD > 3] <- 3
sc_logMP_MPD[sc_logMP_MPD < -3] <- -3
heatmap.2(sc_logMP_MPD, trace ='none', scale = 'none', col=bluered(50), labRow = "", labCol = colnames(sc_logMP_MPD), cexRow = 0.7, cexCol = 0.7, breaks = seq(-3,3, l=51))

# sig genes from Go_res_s_MP_DMSO

MP_Gonad_sample_conditions <- colData(MP_Go_conds_dds)$conds
MP_Gonad_keep_samples <- colnames(MP_Go_conds_dds)[MP_Gonad_sample_conditions %in% c("MP", "DMSO")]
GoMP_Norm_sub <- GoMP_Norm[, MP_Gonad_keep_samples]

GoMP_DMSO <- intersect(rownames(Go_res_s_MP_DMSO), rownames(GoMP_Norm_sub))
GoMP_DMSO_Norm_sig <- GoMP_Norm_sub[GoMP_DMSO, ]

logGoMP_DMSO <- log2(GoMP_DMSO_Norm_sig + 0.5)
sc_logMP_DMSO <- t(scale(t(logGoMP_DMSO)))
sc_logMP_DMSO[sc_logMP_DMSO > 3] <- 3
sc_logMP_DMSO[sc_logMP_DMSO < -3] <- -3
heatmap.2(sc_logMP_DMSO, trace ='none', scale = 'none', col=bluered(50), labRow = "", labCol = colnames(sc_logMP_DMSO), cexRow = 0.7, cexCol = 0.7, breaks = seq(-3,3, l=51))

#sig genes from Go_res_s_MPD_DMSO

GoMPD_DMSO <- intersect(rownames(Go_res_s_MPD_DMSO), rownames(GoMP_Norm))
GoMPD_DMSO_Norm_sig <- GoMP_Norm[GoMPD_DMSO, ]

logGoMPD_DMSO <- log2(GoMPD_DMSO_Norm_sig + 0.5)
sc_logMPD_DMSO <- t(scale(t(logGoMPD_DMSO)))
sc_logMPD_DMSO[sc_logMPD_DMSO > 3] <- 3
sc_logMPD_DMSO[sc_logMPD_DMSO < -3] <- -3
heatmap.2(sc_logMPD_DMSO, trace ='none', scale = 'none', col=bluered(50), labRow = rownames(sc_logMPD_DMSO), labCol = colnames(sc_logMPD_DMSO), cexRow = 0.7, cexCol = 0.7, breaks = seq(-3,3, l=51))

#### Gonad F1 ####

GoF1_Norm <- counts(F1_Go_conds_dds, normalized = TRUE)

GoF1 <- intersect(rownames(Go_res_s_F1), rownames(GoF1_Norm))
GoF1_Norm_sig <- GoF1_Norm[GoF1, ]

logGoF1 <- log2(GoF1_Norm_sig + 0.5)
sc_logF1 <- t(scale(t(logGoF1)))
sc_logF1[sc_logF1 > 3] <- 3
sc_logF1[sc_logF1 < -3] <- -3
heatmap.2(sc_logF1, trace ='none', scale = 'none', col=bluered(50), labRow = "", labCol = colnames(sc_logF1), cexRow = 0.7, cexCol = 0.7, breaks = seq(-3,3, l=51))

#### Liver HLD ####

#sig genes from Li_res_s_LD
LiHLD_Norm <- counts(HLD_Li_conds_dds, normalized = TRUE)
Liver_sample_conditions <- colData(HLD_Li_conds_dds)$conds
Liver_keep_samples <- colnames(HLD_Li_conds_dds)[Liver_sample_conditions %in% c("Low", "DMSO")]
LiHLD_Norm_sub <- LiHLD_Norm[, Liver_keep_samples]

LiLD <- intersect(rownames(Li_res_s_LD), rownames(LiHLD_Norm_sub))
LiLD_Norm_sig <- LiHLD_Norm_sub[LiLD, ]

logLiLD <- log2(LiLD_Norm_sig + 0.5)
sc_logLiLD <- t(scale(t(logLiLD)))
sc_logLiLD[sc_logLiLD > 3] <- 3
sc_logLiLD[sc_logLiLD < -3] <- -3
heatmap.2(sc_logLiLD, trace ='none', scale = 'none', col=bluered(50), labRow = "", labCol = colnames(sc_logLiLD), cexRow = 0.7, cexCol = 0.7, breaks = seq(-3,3, l=51))

#sig genes from Li_res_s_HD
LiHD <- intersect(rownames(Li_res_s_HD), rownames(LiHLD_Norm))
LiHD_Norm_sig <- LiHLD_Norm[LiHD, ]

logLiHD <- log2(LiHD_Norm_sig + 0.5)
sc_logLiHD <- t(scale(t(logLiHD)))
sc_logLiHD[sc_logLiHD > 3] <- 3
sc_logLiHD[sc_logLiHD < -3] <- -3
heatmap.2(sc_logLiHD, trace ='none', scale = 'none', col=bluered(50), labRow = rownames(sc_logLiHD), labCol = colnames(sc_logLiHD), cexRow = 0.7, cexCol = 0.7, breaks = seq(-3,3, l=51))


#sig genes from Li_res_s_LD
LiLD <- intersect(rownames(Li_res_s_LD), rownames(LiHLD_Norm))
LiLD_Norm_sig <- LiHLD_Norm[LiLD, ]

logLiLD <- log2(LiLD_Norm_sig + 0.5)
sc_logLiLD <- t(scale(t(logLiLD)))
sc_logLiLD[sc_logLiLD > 3] <- 3
sc_logLiLD[sc_logLiLD < -3] <- -3
heatmap.2(sc_logLiLD, trace ='none', scale = 'none', col=bluered(50), labRow = "", labCol = colnames(sc_logLiLD), cexRow = 0.7, cexCol = 0.7, breaks = seq(-3,3, l=51))

#### F1 Liver ####
LiF1_Norm <- counts(F1_Li_conds_dds, normalized = TRUE)

LiF1 <- intersect(rownames(Li_res_s_F1), rownames(LiF1_Norm))
LiF1_Norm_sig <- LiF1_Norm[LiF1, ]

logLiF1 <- log2(LiF1_Norm_sig + 0.5)
sc_logLiF1 <- t(scale(t(logLiF1)))
sc_logLiF1[sc_logLiF1 > 3] <- 3
sc_logLiF1[sc_logLiF1 < -3] <- -3
heatmap.2(sc_logLiF1, trace ='none', scale = 'none', col=bluered(50), labRow = "", labCol = colnames(sc_logLiF1), cexRow = 0.7, cexCol = 0.7, breaks = seq(-3,3, l=51))


#Save
#saveRDS(MP_Br_conds_dds, file = "DDS objects/MP_Br_conds_dds.rds")
#saveRDS(F1_Br_conds_dds, file = "DDS objects/F1_Br_conds_dds.rds")
#saveRDS(HLD_Go_conds_dds, file = "DDS objects/HLD_Go_conds_dds.rds")
#saveRDS(MP_Go_conds_dds, file = "DDS objects/MP_Go_conds_dds.rds")
#saveRDS(F1_Go_conds_dds, file = "DDS objects/F1_Go_conds_dds.rds")
#saveRDS(HLD_Li_conds_dds, file = "DDS objects/HLD_Li_conds_dds.rds")
#saveRDS(F1_Li_conds_dds, file = "DDS objects/F1_Li_conds_dds.rds")
