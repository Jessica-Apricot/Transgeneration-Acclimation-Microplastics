#Brain

Brain_cols <- c("MP_F0_M1_B", "D_F0_M1_B", "L_F0_M1_B", "H_F1_M3_B", "H_F1_M3_B.2.", "D_F0_M2_B", "D_F1_M3_B", "MPD_F0_M1_B", "D_F0_M3_B", "L_F0_M2_B", "L_F0_M3_B", "H_F0_M2_B", "H_F0_M1_B", "H_F0_M3_B", "MP_F0_M3_B", "MP_F0_M2_B", "MPD_F0_M2_B", "MPD_F0_M3_B", "D_F1_M1_B", "D_F1_M1_B.2.", "H_F1_M2_B", "MP_F0_M4_B", "H_F0_M4_B", "D_F0_M4_B")
Brain <- Brain_Counts[, Brain_cols]
rownames(Brain) <- rownames(Brain_Counts)


Br_conds <- c("MP", "DMSO", "Low", "High", "High", "DMSO", "DMSO", "MPD", "DMSO", "Low", "Low", "High", "High", "High", "MP", "MP", "MPD", "MPD", "DMSO", "DMSO", "High", "MP", "High", "DMSO" )

Br_metadata <- data.frame(row.names = colnames(Brain),
                              conds = factor(Br_conds))
Br_conds_dds <- DESeqDataSetFromMatrix(countData = as.matrix(Brain),
                                           colData = Br_metadata,
                                           design = ~ conds)

keep <- rowSums(cpm(counts(Br_conds_dds)) >= 1) >= 6
Br_conds_dds <- Br_conds_dds[keep,]

Br_conds_dds <- DESeq(Br_conds_dds)

Br_PCA <- rlog(Br_conds_dds, blind = TRUE) 
Br_rv <- rowVars(assay(Br_PCA)) 
Br_top_500 <- order(Br_rv, decreasing = TRUE)[1:500] 

Br_top_matrix <- assay(Br_PCA)[Br_top_500, ] 

Br_pca_res <- prcomp(t(Br_top_matrix), scale. = TRUE)

Br_metadata <- as.data.frame(colData(Br_conds_dds)) 

Br_pca_df <- as.data.frame(Br_pca_res$x)

Br_pca_df$conds <- Br_metadata$conds 

percentVar <- round(100 * (Br_pca_res$sdev)^2 / sum(Br_pca_res$sdev^2), 1)

#with jitter 
set.seed(9)
Br_pca_df$jitter_PC1 <- Br_pca_df$PC1 + runif(nrow(Br_pca_df), -0.3, 0.3)
Br_pca_df$jitter_PC2 <- Br_pca_df$PC2 + runif(nrow(Br_pca_df), -0.3, 0.3)

ggplot(Br_pca_df, aes(x = jitter_PC1, y = jitter_PC2, color = conds)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  labs(title = "",
       x = paste0("PC1: ", percentVar[1], "%"),
       y = paste0("PC2: ", percentVar[2], "%"))



#Without jitter
ggplot(Br_pca_df, aes(x = PC1, y = PC2, color = conds)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  labs(title = "",
       x = paste0("PC1: ", percentVar[1], "%"),
       y = paste0("PC2: ", percentVar[2], "%")) 


#Gonad PCA

Gonad_cols <- c("D_F0_M2_G", "L_F0_M2_G", "L_F0_M3_G", "H_F0_M2_G","D_F1_M1_G", "D_F1_M1_G.2.", "D_F1_M3_G", "H_F1_M3_G", "H_F1_M3_G.2.", "H_F1_M2_G", "MP_F0_M1_G", "MPD_F0_M1_G", "H_F0_M1_G", "H_F0_M3_G", "MP_F0_M3_G", "MP_F0_M2_G", "MPD_F0_M2_G", "MPD_F0_M3_G", "D_F0_M3_G", "D_F0_M1_G", "L_F0_M1_G", "MP_F0_M4_G", "H_F0_M4_G", "D_F0_M4_G")
Gonad <- Gonad_Counts[, Gonad_cols]
rownames(Gonad) <- rownames(Gonad_Counts)


Go_conds <- c("DMSO", "Low", "Low", "High", "DMSO", "DMSO", "DMSO", "High", "High", "High", "MP", "MPD", "High", "High", "MP", "MP", "MPD", "MPD", "DMSO", "DMSO", "Low", "MP", "High", "DMSO" )

Go_metadata <- data.frame(row.names = colnames(Gonad),
                          conds = factor(Go_conds))
Go_conds_dds <- DESeqDataSetFromMatrix(countData = as.matrix(Gonad),
                                       colData = Go_metadata,
                                       design = ~ conds)

keep <- rowSums(cpm(counts(Go_conds_dds)) >= 1) >= 6
Go_conds_dds <- Go_conds_dds[keep,]

Go_conds_dds <- DESeq(Go_conds_dds)

Go_PCA <- rlog(Go_conds_dds, blind = TRUE) 
Go_rv <- rowVars(assay(Go_PCA)) 
Go_top_500 <- order(Go_rv, decreasing = TRUE)[1:500] 

Go_top_matrix <- assay(Go_PCA)[Go_top_500, ] 

Go_pca_res <- prcomp(t(Go_top_matrix), scale. = TRUE)

Go_metadata <- as.data.frame(colData(Go_conds_dds)) 

Go_pca_df <- as.data.frame(Go_pca_res$x)

Go_pca_df$conds <- Go_metadata$conds 

percentVar <- round(100 * (Go_pca_res$sdev)^2 / sum(Go_pca_res$sdev^2), 1)

#with jitter 
set.seed(9)
Go_pca_df$jitter_PC1 <- Go_pca_df$PC1 + runif(nrow(Go_pca_df), -0.3, 0.3)
Go_pca_df$jitter_PC2 <- Go_pca_df$PC2 + runif(nrow(Go_pca_df), -0.3, 0.3)

ggplot(Go_pca_df, aes(x = jitter_PC1, y = jitter_PC2, color = conds)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  labs(title = "",
       x = paste0("PC1: ", percentVar[1], "%"),
       y = paste0("PC2: ", percentVar[2], "%"))



#Without jitter
ggplot(Go_pca_df, aes(x = PC1, y = PC2, color = conds)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  labs(title = "",
       x = paste0("PC1: ", percentVar[1], "%"),
       y = paste0("PC2: ", percentVar[2], "%")) 

#Liver

Liver_cols <- c("L_F0_M1_L", "L_F0_M2_L", "L_F0_M3_L", "H_F0_M1_L", "MP_F0_M2_L", "D_F1_M3_L", "H_F1_M3_L", "H_F1_M3_L.2.", "H_F1_M2_L", "MP_F0_M1_L", "MPD_F0_M1_L", "D_F0_M3_L", "D_F0_M2_L", "D_F0_M1_L", "H_F0_M2_L", "H_F0_M3_L", "MP_F0_M3_L", "MPD_F0_M2_L", "MPD_F0_M3_L", "D_F1_M1_L", "D_F1_M1_L.2.", "MP_F0_M4_L", "H_F0_M4_L", "D_F0_M4_L")
Liver <- Liver_Counts[, Liver_cols]
rownames(Liver) <- rownames(Liver_Counts)


Li_conds <- c("Low", "Low", "Low", "High", "MP", "DMSO", "High", "High", "High", "MP", "MPD", "DMSO", "DMSO", "DMSO", "High", "High", "MP", "MPD", "MPD", "DMSO", "DMSO", "MP", "High", "DMSO" )

Li_metadata <- data.frame(row.names = colnames(Liver),
                          conds = factor(Li_conds))
Li_conds_dds <- DESeqDataSetFromMatrix(countData = as.matrix(Liver),
                                       colData = Li_metadata,
                                       design = ~ conds)

keep <- rowSums(cpm(counts(Li_conds_dds)) >= 1) >= 6
Li_conds_dds <- Li_conds_dds[keep,]

Li_conds_dds <- DESeq(Li_conds_dds)

Li_PCA <- rlog(Li_conds_dds, blind = TRUE) 
Li_rv <- rowVars(assay(Li_PCA)) 
Li_top_500 <- order(Li_rv, decreasing = TRUE)[1:500] 

Li_top_matrix <- assay(Li_PCA)[Li_top_500, ] 

Li_pca_res <- prcomp(t(Li_top_matrix), scale. = TRUE)

Li_metadata <- as.data.frame(colData(Li_conds_dds)) 

Li_pca_df <- as.data.frame(Li_pca_res$x)

Li_pca_df$conds <- Li_metadata$conds 

percentVar <- round(100 * (Li_pca_res$sdev)^2 / sum(Li_pca_res$sdev^2), 1)

#with jitter 
set.seed(9)
Li_pca_df$jitter_PC1 <- Li_pca_df$PC1 + runif(nrow(Li_pca_df), -0.3, 0.3)
Li_pca_df$jitter_PC2 <- Li_pca_df$PC2 + runif(nrow(Li_pca_df), -0.3, 0.3)

ggplot(Li_pca_df, aes(x = jitter_PC1, y = jitter_PC2, color = conds)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  labs(title = "",
       x = paste0("PC1: ", percentVar[1], "%"),
       y = paste0("PC2: ", percentVar[2], "%"))



#Without jitter
ggplot(Li_pca_df, aes(x = PC1, y = PC2, color = conds)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  labs(title = "",
       x = paste0("PC1: ", percentVar[1], "%"),
       y = paste0("PC2: ", percentVar[2], "%")) 

