library(gplots)
library(DESeq2)

#### Brain Microplastics heatmaps ####

# sig genes from Br.res_s_MP.DMSO
#Rows = genes
#Columns = samples

#Log the data

BrMP_Norm <- counts(MP.Br_conds_dds, normalized = TRUE)

BrMP.DMSO <- intersect(rownames(Br.res_s_MP.DMSO), rownames(BrMP_Norm))
BrMP.DMSO_Norm_sig <- BrMP_Norm[BrMP.DMSO, ]

logBrMP <- log2(BrMP.DMSO_Norm_sig + 0.5)
sc_logBrMP <- t(scale(t(logBrMP)))
sc_logBrMP[sc_logBrMP > 3] <- 3
sc_logBrMP[sc_logBrMP < -3] <- -3
heatmap.2(sc_logBrMP, trace ='none', scale = 'none', col=bluered(50), labRow = rownames(sc_logBrMP), labCol = colnames(sc_logBrMP), cexRow = 0.7, cexCol = 0.7, breaks = seq(-3,3, l=51))


# sig genes from Br.res_s_MPD.DMSO

BrMPD.DMSO<- intersect(rownames(Br.res_s_MPD.DMSO), rownames(BrMP_Norm))
BrMPD.DMSO_Norm_sig <- BrMP_Norm[BrMPD.DMSO, ]

logBrMPD <- log2(BrMPD.DMSO_Norm_sig + 0.5)
sc_logBrMPD <- t(scale(t(logBrMPD)))
sc_logBrMPD[sc_logBrMPD > 3] <- 3
sc_logBrMPD[sc_logBrMPD < -3] <- -3
heatmap.2(sc_logBrMPD, trace ='none', scale = 'none', col=bluered(50), labRow = rownames(sc_logBrMPD), labCol = colnames(sc_logBrMPD), cexRow = 0.7, cexCol = 0.7, breaks = seq(-3,3, l=51))

#### Brain F1 ####
#sig genes from Br.res_s_F1

BrF1_Norm <- counts(F1.Br_conds_dds, normalized = TRUE)

BrF1 <- intersect(rownames(Br.res_s_F1), rownames(BrF1_Norm))
BrF1_Norm_sig <- BrF1_Norm[BrF1, ]

logBrF1 <- log2(BrF1_Norm_sig + 0.5)
sc_logBrF1 <- t(scale(t(logBrF1)))
sc_logBrF1[sc_logBrF1 > 3] <- 3
sc_logBrF1[sc_logBrF1 < -3] <- -3
heatmap.2(sc_logBrF1, trace ='none', scale = 'none', col=bluered(50), labRow = "", labCol = colnames(sc_logBrF1), cexRow = 0.7, cexCol = 0.7, breaks = seq(-3,3, l=51))

#### Gonad High-Low-DMSO ####
#sig genes from Go.res_s_HL

GoHLD_Norm <- counts(HLD.Go_conds_dds, normalized = TRUE)

GoHL <- intersect(rownames(Go.res_s_HL), rownames(GoHLD_Norm))
GoHL_Norm_sig <- GoHLD_Norm[GoHL, ]

logGoHL <- log2(GoHL_Norm_sig + 0.5)
sc_logGoHL <- t(scale(t(logGoHL)))
sc_logGoHL[sc_logGoHL > 3] <- 3
sc_logGoHL[sc_logGoHL < -3] <- -3
heatmap.2(sc_logGoHL, trace ='none', scale = 'none', col=bluered(50), labRow = "", labCol = colnames(sc_logGoHL), cexRow = 0.7, cexCol = 0.7, breaks = seq(-3,3, l=51))

#sig genes from Go.res_s_HD

GoHD <- intersect(rownames(Go.res_s_HD), rownames(GoHLD_Norm))
GoHD_Norm_sig <- GoHLD_Norm[GoHD, ]

logGoHD <- log2(GoHD_Norm_sig + 0.5)
sc_logGoHD <- t(scale(t(logGoHD)))
sc_logGoHD[sc_logGoHD > 3] <- 3
sc_logGoHD[sc_logGoHD < -3] <- -3
heatmap.2(sc_logGoHD, trace ='none', scale = 'none', col=bluered(50), labRow = "", labCol = colnames(sc_logGoHD), cexRow = 0.7, cexCol = 0.7, breaks = seq(-3,3, l=51))

#sig genes from Go.res_s_LD

GoLD <- intersect(rownames(Go.res_s_LD), rownames(GoHLD_Norm))
GoLD_Norm_sig <- GoHLD_Norm[GoLD, ]

logGoLD <- log2(GoLD_Norm_sig + 0.5)
sc_logGoLD <- t(scale(t(logGoLD)))
sc_logGoLD[sc_logGoLD > 3] <- 3
sc_logGoLD[sc_logGoLD < -3] <- -3
heatmap.2(sc_logGoLD, trace ='none', scale = 'none', col=bluered(50), labRow = rownames(sc_logGoLD), labCol = colnames(sc_logGoLD), cexRow = 0.7, cexCol = 0.7, breaks = seq(-3,3, l=51))

#### Gonad Microplastics ####

# sig genes from Go.res_s_MP.MPD

GoMP_Norm <- counts(MP.Go_conds_dds, normalized = TRUE)

GoMP.MPD <- intersect(rownames(Go.res_s_MP.MPD), rownames(GoMP_Norm))
GoMP.MPD_Norm_sig <- GoMP_Norm[GoMP.MPD, ]

logGoMP.MPD <- log2(GoMP.MPD_Norm_sig + 0.5)
sc_logMP.MPD <- t(scale(t(logGoMP.MPD)))
sc_logMP.MPD[sc_logMP.MPD > 3] <- 3
sc_logMP.MPD[sc_logMP.MPD < -3] <- -3
heatmap.2(sc_logMP.MPD, trace ='none', scale = 'none', col=bluered(50), labRow = "", labCol = colnames(sc_logMP.MPD), cexRow = 0.7, cexCol = 0.7, breaks = seq(-3,3, l=51))

# sig genes from Go.res_s_MP.DMSO

GoMP.DMSO <- intersect(rownames(Go.res_s_MP.DMSO), rownames(GoMP_Norm))
GoMP.DMSO_Norm_sig <- GoMP_Norm[GoMP.DMSO, ]

logGoMP.DMSO <- log2(GoMP.DMSO_Norm_sig + 0.5)
sc_logMP.DMSO <- t(scale(t(logGoMP.DMSO)))
sc_logMP.DMSO[sc_logMP.DMSO > 3] <- 3
sc_logMP.DMSO[sc_logMP.DMSO < -3] <- -3
heatmap.2(sc_logMP.DMSO, trace ='none', scale = 'none', col=bluered(50), labRow = "", labCol = colnames(sc_logMP.DMSO), cexRow = 0.7, cexCol = 0.7, breaks = seq(-3,3, l=51))

#sig genes from Go.res_s_MPD.DMSO

GoMPD.DMSO <- intersect(rownames(Go.res_s_MPD.DMSO), rownames(GoMP_Norm))
GoMPD.DMSO_Norm_sig <- GoMP_Norm[GoMPD.DMSO, ]

logGoMPD.DMSO <- log2(GoMPD.DMSO_Norm_sig + 0.5)
sc_logMPD.DMSO <- t(scale(t(logGoMPD.DMSO)))
sc_logMPD.DMSO[sc_logMPD.DMSO > 3] <- 3
sc_logMPD.DMSO[sc_logMPD.DMSO < -3] <- -3
heatmap.2(sc_logMPD.DMSO, trace ='none', scale = 'none', col=bluered(50), labRow = rownames(sc_logMPD.DMSO), labCol = colnames(sc_logMPD.DMSO), cexRow = 0.7, cexCol = 0.7, breaks = seq(-3,3, l=51))

#### Gonad F1 ####

GoF1_Norm <- counts(F1.Go_conds_dds, normalized = TRUE)

GoF1 <- intersect(rownames(Go.res_s_F1), rownames(GoF1_Norm))
GoF1_Norm_sig <- GoF1_Norm[GoF1, ]

logGoF1 <- log2(GoF1_Norm_sig + 0.5)
sc_logF1 <- t(scale(t(logGoF1)))
sc_logF1[sc_logF1 > 3] <- 3
sc_logF1[sc_logF1 < -3] <- -3
heatmap.2(sc_logF1, trace ='none', scale = 'none', col=bluered(50), labRow = "", labCol = colnames(sc_logF1), cexRow = 0.7, cexCol = 0.7, breaks = seq(-3,3, l=51))

#### Liver HLD ####

#sig genes from Li.res_s_HL
LiHLD_Norm <- counts(HLD.Li_conds_dds, normalized = TRUE)

LiHL <- intersect(rownames(Li.res_s_HL), rownames(LiHLD_Norm))
LiHL_Norm_sig <- LiHLD_Norm[LiHL, ]

logLiHL <- log2(LiHL_Norm_sig + 0.5)
sc_logLiHL <- t(scale(t(logLiHL)))
sc_logLiHL[sc_logLiHL > 3] <- 3
sc_logLiHL[sc_logLiHL < -3] <- -3
heatmap.2(sc_logLiHL, trace ='none', scale = 'none', col=bluered(50), labRow = rownames(sc_logLiHL), labCol = colnames(sc_logLiHL), cexRow = 0.7, cexCol = 0.7, breaks = seq(-3,3, l=51))

#sig genes from Li.res_s_HD
LiHD <- intersect(rownames(Li.res_s_HD), rownames(LiHLD_Norm))
LiHD_Norm_sig <- LiHLD_Norm[LiHD, ]

logLiHD <- log2(LiHD_Norm_sig + 0.5)
sc_logLiHD <- t(scale(t(logLiHD)))
sc_logLiHD[sc_logLiHD > 3] <- 3
sc_logLiHD[sc_logLiHD < -3] <- -3
heatmap.2(sc_logLiHD, trace ='none', scale = 'none', col=bluered(50), labRow = rownames(sc_logLiHD), labCol = colnames(sc_logLiHD), cexRow = 0.7, cexCol = 0.7, breaks = seq(-3,3, l=51))


#sig genes from Li.res_s_LD
LiLD <- intersect(rownames(Li.res_s_LD), rownames(LiHLD_Norm))
LiLD_Norm_sig <- LiHLD_Norm[LiLD, ]

logLiLD <- log2(LiLD_Norm_sig + 0.5)
sc_logLiLD <- t(scale(t(logLiLD)))
sc_logLiLD[sc_logLiLD > 3] <- 3
sc_logLiLD[sc_logLiLD < -3] <- -3
heatmap.2(sc_logLiLD, trace ='none', scale = 'none', col=bluered(50), labRow = "", labCol = colnames(sc_logLiLD), cexRow = 0.7, cexCol = 0.7, breaks = seq(-3,3, l=51))

#### F1 Liver ####
LiF1_Norm <- counts(F1.Li_conds_dds, normalized = TRUE)

LiF1 <- intersect(rownames(Li.res_s_F1), rownames(LiF1_Norm))
LiF1_Norm_sig <- LiF1_Norm[LiF1, ]

logLiF1 <- log2(LiF1_Norm_sig + 0.5)
sc_logLiF1 <- t(scale(t(logLiF1)))
sc_logLiF1[sc_logLiF1 > 3] <- 3
sc_logLiF1[sc_logLiF1 < -3] <- -3
heatmap.2(sc_logLiF1, trace ='none', scale = 'none', col=bluered(50), labRow = "", labCol = colnames(sc_logLiF1), cexRow = 0.7, cexCol = 0.7, breaks = seq(-3,3, l=51))


#Save
saveRDS(MP.Br_conds_dds, file = "DDS objects/MP.Br_conds_dds.rds")
saveRDS(F1.Br_conds_dds, file = "DDS objects/F1.Br_conds_dds.rds")
saveRDS(HLD.Go_conds_dds, file = "DDS objects/HLD.Go_conds_dds.rds")
saveRDS(MP.Go_conds_dds, file = "DDS objects/MP.Go_conds_dds.rds")
saveRDS(F1.Go_conds_dds, file = "DDS objects/F1.Go_conds_dds.rds")
saveRDS(HLD.Li_conds_dds, file = "DDS objects/HLD.Li_conds_dds.rds")
saveRDS(F1.Li_conds_dds, file = "DDS objects/F1.Li_conds_dds.rds")
