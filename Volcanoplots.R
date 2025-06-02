#Volcano plots
library(ggplot2)
#BiocManager::install("ggiraph")
library(ggiraph)
#BiocManager::install("plotly")
library(plotly)

#Define significance thresholds
pval_thres <- 0.05
fc_thres <- 1 #log2 fold

### Brain High-Low ####
#plotting df
volc_Br.HL <- data.frame(
  logFC = Br.res_HL$log2FoldChange,
  negLogPval = -log10(Br.res_HL$padj),
  adj.P.Val = Br.res_HL$padj,
  ID = rownames(Br.res_HL))

#Column to catagorise genes, if above fc threshold labelled upregulated (greater than 1), if below threshold (-1), downregulated
volc_Br.HL$category <- ifelse(Br.res_HL$padj <= pval_thres,
                              ifelse(volc_Br.HL$logFC >= fc_thres, "Upregulated", 
                                     ifelse(volc_Br.HL$logFC <= -fc_thres, "Downregulated", "Not Significant")),
                              "Not Significant")

# Create the volcano plot using ggplot2
ggplot(volc_Br.HL, aes(x = logFC, y = negLogPval, color = category)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  geom_vline(xintercept = c(-fc_thres, fc_thres), linetype = "dashed") +
  geom_hline(yintercept = -log10(pval_thres), linetype = "dashed") +
  labs(
    title = "Differential Gene Expression in High DEHP vs Low DEHP in Brain",
    subtitle = paste("Thresholds: |log2FC| >", fc_thres, "and adjusted p-value <", pval_thres),
    x = "log2 Fold Change",
    y = "-log10(p-value)",
    color = "Differential Expression"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12)
  )




### Brain High-DMSO #####
#plotting df
volc_Br.HD <- data.frame(
  logFC = Br.res_HD$log2FoldChange,
  negLogPval = -log10(Br.res_HD$padj),
  adj.P.Val = Br.res_HD$padj,
  ID = rownames(Br.res_HD))

#Column to catagorise genes, if above fc threshold labelled upregulated (greater than 1), if below threshold (-1), downregulated
volc_Br.HD$category <- ifelse(Br.res_HD$padj <= pval_thres,
                              ifelse(volc_Br.HD$logFC >= fc_thres, "Upregulated", 
                                     ifelse(volc_Br.HD$logFC <= -fc_thres, "Downregulated", "Not Significant")),
                              "Not Significant")

# Create the volcano plot using ggplot2
ggplot(volc_Br.HD, aes(x = logFC, y = negLogPval, color = category)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  geom_vline(xintercept = c(-fc_thres, fc_thres), linetype = "dashed") +
  geom_hline(yintercept = -log10(pval_thres), linetype = "dashed") +
  labs(
    title = "Differential Gene Expression in High DEHP vs DMSO control in Brain",
    subtitle = paste("Thresholds: |log2FC| >", fc_thres, "and adjusted p-value <", pval_thres),
    x = "log2 Fold Change",
    y = "-log10(p-value)",
    color = "Differential Expression"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12)
  )


### Brain Low-DMSO ####
#plotting df
volc_Br.LD <- data.frame(
  logFC = Br.res_LD$log2FoldChange,
  negLogPval = -log10(Br.res_LD$padj),
  adj.P.Val = Br.res_LD$padj,
  ID = rownames(Br.res_LD))

#Column to catagorise genes, if above fc threshold labelled upregulated (greater than 1), if below threshold (-1), downregulated
volc_Br.LD$category <- ifelse(Br.res_LD$padj <= pval_thres,
                              ifelse(volc_Br.LD$logFC >= fc_thres, "Upregulated", 
                                     ifelse(volc_Br.LD$logFC <= -fc_thres, "Downregulated", "Not Significant")),
                              "Not Significant")

# Create the volcano plot using ggplot2
ggplot(volc_Br.LD, aes(x = logFC, y = negLogPval, color = category)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  geom_vline(xintercept = c(-fc_thres, fc_thres), linetype = "dashed") +
  geom_hline(yintercept = -log10(pval_thres), linetype = "dashed") +
  labs(
    title = "Differential Gene Expression in Low DEHP vs DMSO control in Brain",
    subtitle = paste("Thresholds: |log2FC| >", fc_thres, "and adjusted p-value <", pval_thres),
    x = "log2 Fold Change",
    y = "-log10(p-value)",
    color = "Differential Expression"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12)
  )



### Brain MP- MPD ####
#plotting df
volc_Br.MP.MPD <- data.frame(
  logFC = Br.res_MP.MPD$log2FoldChange,
  negLogPval = -log10(Br.res_MP.MPD$padj),
  adj.P.Val = Br.res_MP.MPD$padj,
  ID = rownames(Br.res_MP.MPD))

#Column to catagorise genes, if above fc threshold labelled upregulated (greater than 1), if below threshold (-1), downregulated
volc_Br.MP.MPD$category <- ifelse(Br.res_MP.MPD$padj <= pval_thres,
                              ifelse(volc_Br.MP.MPD$logFC >= fc_thres, "Upregulated", 
                                     ifelse(volc_Br.MP.MPD$logFC <= -fc_thres, "Downregulated", "Not Significant")),
                              "Not Significant")

# Create the volcano plot using ggplot2
ggplot(volc_Br.MP.MPD, aes(x = logFC, y = negLogPval, color = category)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  geom_vline(xintercept = c(-fc_thres, fc_thres), linetype = "dashed") +
  geom_hline(yintercept = -log10(pval_thres), linetype = "dashed") +
  labs(
    title = "Differential Gene Expression in Microplastics control vs Microplastics+DEHP in Brain",
    subtitle = paste("Thresholds: |log2FC| >", fc_thres, "and adjusted p-value <", pval_thres),
    x = "log2 Fold Change",
    y = "-log10(p-value)",
    color = "Differential Expression"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12) )


### Brain MP-DMSO ####
#plotting df
volc_Br.MP.DMSO <- data.frame(
  logFC = Br.res_MP.DMSO$log2FoldChange,
  negLogPval = -log10(Br.res_MP.DMSO$padj),
  adj.P.Val = Br.res_MP.DMSO$padj,
  ID = rownames(Br.res_MP.DMSO))

#Column to catagorise genes, if above fc threshold labelled upregulated (greater than 1), if below threshold (-1), downregulated
volc_Br.MP.DMSO$category <- ifelse(Br.res_MP.DMSO$padj <= pval_thres,
                                  ifelse(volc_Br.MP.DMSO$logFC >= fc_thres, "Upregulated", 
                                         ifelse(volc_Br.MP.DMSO$logFC <= -fc_thres, "Downregulated", "Not Significant")),
                                  "Not Significant")

# Create the volcano plot using ggplot2
ggplot(volc_Br.MP.DMSO, aes(x = logFC, y = negLogPval, color = category)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  geom_vline(xintercept = c(-fc_thres, fc_thres), linetype = "dashed") +
  geom_hline(yintercept = -log10(pval_thres), linetype = "dashed") +
  labs(
    title = "Differential Gene Expression in Microplastics control vs DMSO in Brain",
    subtitle = paste("Thresholds: |log2FC| >", fc_thres, "and adjusted p-value <", pval_thres),
    x = "log2 Fold Change",
    y = "-log10(p-value)",
    color = "Differential Expression"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12) )



### Brain MPD-DMSO ####
#One of the sig genes are overlapping
#plotting df
volc_Br.MPD.DMSO <- data.frame(
  logFC = Br.res_MPD.DMSO$log2FoldChange,
  negLogPval = -log10(Br.res_MPD.DMSO$padj),
  adj.P.Val = Br.res_MPD.DMSO$padj,
  ID = rownames(Br.res_MPD.DMSO))

#Column to catagorise genes, if above fc threshold labelled upregulated (greater than 1), if below threshold (-1), downregulated
volc_Br.MPD.DMSO$category <- ifelse(Br.res_MPD.DMSO$padj <= pval_thres,
                                   ifelse(volc_Br.MPD.DMSO$logFC >= fc_thres, "Upregulated", 
                                          ifelse(volc_Br.MPD.DMSO$logFC <= -fc_thres, "Downregulated", "Not Significant")),
                                   "Not Significant")

# Create the volcano plot using ggplot2
ggplot(volc_Br.MPD.DMSO, aes(x = logFC, y = negLogPval, color = category)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  geom_vline(xintercept = c(-fc_thres, fc_thres), linetype = "dashed") +
  geom_hline(yintercept = -log10(pval_thres), linetype = "dashed") +
  labs(
    title = "Differential Gene Expression in Microplastics + DEHP vs DMSO in Brain",
    subtitle = paste("Thresholds: |log2FC| >", fc_thres, "and adjusted p-value <", pval_thres),
    x = "log2 Fold Change",
    y = "-log10(p-value)",
    color = "Differential Expression"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12) )


### Brain F1 ####
#plotting df
volc_Br.F1 <- data.frame(
  logFC = Br.res_F1$log2FoldChange,
  negLogPval = -log10(Br.res_F1$padj),
  adj.P.Val = Br.res_F1$padj,
  ID = rownames(Br.res_F1))

#Column to catagorise genes, if above fc threshold labelled upregulated (greater than 1), if below threshold (-1), downregulated
volc_Br.F1$category <- ifelse(Br.res_F1$padj <= pval_thres,
                                   ifelse(volc_Br.F1$logFC >= fc_thres, "Upregulated", 
                                          ifelse(volc_Br.F1$logFC <= -fc_thres, "Downregulated", "Not Significant")),
                                   "Not Significant")

# Create the volcano plot using ggplot2
ggplot(volc_Br.F1, aes(x = logFC, y = negLogPval, color = category)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  geom_vline(xintercept = c(-fc_thres, fc_thres), linetype = "dashed") +
  geom_hline(yintercept = -log10(pval_thres), linetype = "dashed") +
  labs(
    title = "Differential Gene Expression in Brain F1 High vs DMSO control",
    subtitle = paste("Thresholds: |log2FC| >", fc_thres, "and adjusted p-value <", pval_thres),
    x = "log2 Fold Change",
    y = "-log10(p-value)",
    color = "Differential Expression"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12) )


### Gonad High-Low ####
#plotting df
volc_Go.HL <- data.frame(
  logFC = Go.res_HL$log2FoldChange,
  negLogPval = -log10(Go.res_HL$padj),
  adj.P.Val = Go.res_HL$padj,
  ID = rownames(Go.res_HL))

#Column to catagorise genes, if above fc threshold labelled upregulated (greater than 1), if below threshold (-1), downregulated
volc_Go.HL$category <- ifelse(Go.res_HL$padj <= pval_thres,
                              ifelse(volc_Go.HL$logFC >= fc_thres, "Upregulated", 
                                     ifelse(volc_Go.HL$logFC <= -fc_thres, "Downregulated", "Not Significant")),
                              "Not Significant")

# Create the volcano plot using ggplot2
ggplot(volc_Go.HL, aes(x = logFC, y = negLogPval, color = category)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  geom_vline(xintercept = c(-fc_thres, fc_thres), linetype = "dashed") +
  geom_hline(yintercept = -log10(pval_thres), linetype = "dashed") +
  labs(
    title = "Differential Gene Expression in High DEHP vs Low DEHP in Gonads",
    subtitle = paste("Thresholds: |log2FC| >", fc_thres, "and adjusted p-value <", pval_thres),
    x = "log2 Fold Change",
    y = "-log10(p-value)",
    color = "Differential Expression"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12)
  )

### Gonad High-DMSO ####
#plotting df
volc_Go.HD <- data.frame(
  logFC = Go.res_HD$log2FoldChange,
  negLogPval = -log10(Go.res_HD$padj),
  adj.P.Val = Go.res_HD$padj,
  ID = rownames(Go.res_HD))

#Column to catagorise genes, if above fc threshold labelled upregulated (greater than 1), if below threshold (-1), downregulated
volc_Go.HD$category <- ifelse(Go.res_HD$padj <= pval_thres,
                              ifelse(volc_Go.HD$logFC >= fc_thres, "Upregulated", 
                                     ifelse(volc_Go.HD$logFC <= -fc_thres, "Downregulated", "Not Significant")),
                              "Not Significant")

# Create the volcano plot using ggplot2
ggplot(volc_Go.HD, aes(x = logFC, y = negLogPval, color = category)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  geom_vline(xintercept = c(-fc_thres, fc_thres), linetype = "dashed") +
  geom_hline(yintercept = -log10(pval_thres), linetype = "dashed") +
  labs(
    title = "Differential Gene Expression in High DEHP vs DMSO control in Gonad",
    subtitle = paste("Thresholds: |log2FC| >", fc_thres, "and adjusted p-value <", pval_thres),
    x = "log2 Fold Change",
    y = "-log10(p-value)",
    color = "Differential Expression"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12)
  )

### Gonad Low-DMSO ####
#plotting df
volc_Go.LD <- data.frame(
  logFC = Go.res_LD$log2FoldChange,
  negLogPval = -log10(Go.res_LD$padj),
  adj.P.Val = Go.res_LD$padj,
  ID = rownames(Go.res_LD))

#Column to catagorise genes, if above fc threshold labelled upregulated (greater than 1), if below threshold (-1), downregulated
volc_Go.LD$category <- ifelse(Go.res_LD$padj <= pval_thres,
                              ifelse(volc_Go.LD$logFC >= fc_thres, "Upregulated", 
                                     ifelse(volc_Go.LD$logFC <= -fc_thres, "Downregulated", "Not Significant")),
                              "Not Significant")

# Create the volcano plot using ggplot2
ggplot(volc_Go.LD, aes(x = logFC, y = negLogPval, color = category)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  geom_vline(xintercept = c(-fc_thres, fc_thres), linetype = "dashed") +
  geom_hline(yintercept = -log10(pval_thres), linetype = "dashed") +
  labs(
    title = "Differential Gene Expression in Low DEHP vs DMSO control in Gonad",
    subtitle = paste("Thresholds: |log2FC| >", fc_thres, "and adjusted p-value <", pval_thres),
    x = "log2 Fold Change",
    y = "-log10(p-value)",
    color = "Differential Expression"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12)
  )

### Gonad MP- MPD ####
volc_Go.MP.MPD <- data.frame(
  logFC = Go.res_MP.MPD$log2FoldChange,
  negLogPval = -log10(Go.res_MP.MPD$padj),
  adj.P.Val = Go.res_MP.MPD$padj,
  ID = rownames(Go.res_MP.MPD))

#Column to catagorise genes, if above fc threshold labelled upregulated (greater than 1), if below threshold (-1), downregulated
volc_Go.MP.MPD$category <- ifelse(Go.res_MP.MPD$padj <= pval_thres,
                              ifelse(volc_Go.MP.MPD$logFC >= fc_thres, "Upregulated", 
                                     ifelse(volc_Go.MP.MPD$logFC <= -fc_thres, "Downregulated", "Not Significant")),
                              "Not Significant")

# Create the volcano plot using ggplot2
ggplot(volc_Go.LD, aes(x = logFC, y = negLogPval, color = category)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  geom_vline(xintercept = c(-fc_thres, fc_thres), linetype = "dashed") +
  geom_hline(yintercept = -log10(pval_thres), linetype = "dashed") +
  labs(
    title = "Differential Gene Expression in Microplastic Control vs Microplastics+DEHP in Gonad",
    subtitle = paste("Thresholds: |log2FC| >", fc_thres, "and adjusted p-value <", pval_thres),
    x = "log2 Fold Change",
    y = "-log10(p-value)",
    color = "Differential Expression"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12)
  )


### Gonad MP- DMSO ####
volc_Go.MP.DMSO <- data.frame(
  logFC = Go.res_MP.DMSO$log2FoldChange,
  negLogPval = -log10(Go.res_MP.DMSO$padj),
  adj.P.Val = Go.res_MP.DMSO$padj,
  ID = rownames(Go.res_MP.DMSO))

#Column to catagorise genes, if above fc threshold labelled upregulated (greater than 1), if below threshold (-1), downregulated
volc_Go.MP.DMSO$category <- ifelse(Go.res_MP.DMSO$padj <= pval_thres,
                                  ifelse(volc_Go.MP.DMSO$logFC >= fc_thres, "Upregulated", 
                                         ifelse(volc_Go.MP.DMSO$logFC <= -fc_thres, "Downregulated", "Not Significant")),
                                  "Not Significant")

# Create the volcano plot using ggplot2
ggplot(volc_Go.MP.DMSO, aes(x = logFC, y = negLogPval, color = category)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  geom_vline(xintercept = c(-fc_thres, fc_thres), linetype = "dashed") +
  geom_hline(yintercept = -log10(pval_thres), linetype = "dashed") +
  labs(
    title = "Differential Gene Expression in Microplastic Control vs DMSO control in Gonad",
    subtitle = paste("Thresholds: |log2FC| >", fc_thres, "and adjusted p-value <", pval_thres),
    x = "log2 Fold Change",
    y = "-log10(p-value)",
    color = "Differential Expression"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12)
  )


### Gonad MPD - DMSO ####
volc_Go.MPD.DMSO <- data.frame(
  logFC = Go.res_MPD.DMSO$log2FoldChange,
  negLogPval = -log10(Go.res_MPD.DMSO$padj),
  adj.P.Val = Go.res_MPD.DMSO$padj,
  ID = rownames(Go.res_MPD.DMSO))

#Column to catagorise genes, if above fc threshold labelled upregulated (greater than 1), if below threshold (-1), downregulated
volc_Go.MPD.DMSO$category <- ifelse(Go.res_MPD.DMSO$padj <= pval_thres,
                                   ifelse(volc_Go.MPD.DMSO$logFC >= fc_thres, "Upregulated", 
                                          ifelse(volc_Go.MPD.DMSO$logFC <= -fc_thres, "Downregulated", "Not Significant")),
                                   "Not Significant")

# Create the volcano plot using ggplot2
ggplot(volc_Go.MPD.DMSO, aes(x = logFC, y = negLogPval, color = category)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  geom_vline(xintercept = c(-fc_thres, fc_thres), linetype = "dashed") +
  geom_hline(yintercept = -log10(pval_thres), linetype = "dashed") +
  labs(
    title = "Differential Gene Expression in Microplastics+DEHP vs DMSO control in Gonad",
    subtitle = paste("Thresholds: |log2FC| >", fc_thres, "and adjusted p-value <", pval_thres),
    x = "log2 Fold Change",
    y = "-log10(p-value)",
    color = "Differential Expression"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12)
  )

### Gonad F1 ####
volc_Go.F1 <- data.frame(
  logFC = Go.res_F1$log2FoldChange,
  negLogPval = -log10(Go.res_F1$padj),
  adj.P.Val = Go.res_F1$padj,
  ID = rownames(Go.res_F1))

#Column to catagorise genes, if above fc threshold labelled upregulated (greater than 1), if below threshold (-1), downregulated
volc_Go.F1$category <- ifelse(Go.res_F1$padj <= pval_thres,
                                    ifelse(volc_Go.F1$logFC >= fc_thres, "Upregulated", 
                                           ifelse(volc_Go.F1$logFC <= -fc_thres, "Downregulated", "Not Significant")),
                                    "Not Significant")

# Create the volcano plot using ggplot2
ggplot(volc_Go.F1, aes(x = logFC, y = negLogPval, color = category)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  geom_vline(xintercept = c(-fc_thres, fc_thres), linetype = "dashed") +
  geom_hline(yintercept = -log10(pval_thres), linetype = "dashed") +
  labs(
    title = "Differential Gene Expression in Gonad F1 High vs DMSO control",
    subtitle = paste("Thresholds: |log2FC| >", fc_thres, "and adjusted p-value <", pval_thres),
    x = "log2 Fold Change",
    y = "-log10(p-value)",
    color = "Differential Expression"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12)
  )

### Liver High-Low ####
#plotting df
volc_Li.HL <- data.frame(
  logFC = Li.res_HL$log2FoldChange,
  negLogPval = -log10(Li.res_HL$padj),
  adj.P.Val = Li.res_HL$padj,
  ID = rownames(Li.res_HL))

#Column to catagorise genes, if above fc threshold labelled upregulated (greater than 1), if below threshold (-1), downregulated
volc_Li.HL$category <- ifelse(Li.res_HL$padj <= pval_thres,
                              ifelse(volc_Li.HL$logFC >= fc_thres, "Upregulated", 
                                     ifelse(volc_Li.HL$logFC <= -fc_thres, "Downregulated", "Not Significant")),
                              "Not Significant")

# Create the volcano plot using ggplot2
ggplot(volc_Li.HL, aes(x = logFC, y = negLogPval, color = category)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  geom_vline(xintercept = c(-fc_thres, fc_thres), linetype = "dashed") +
  geom_hline(yintercept = -log10(pval_thres), linetype = "dashed") +
  labs(
    title = "Differential Gene Expression in High DEHP vs Low DEHP in Liver",
    subtitle = paste("Thresholds: |log2FC| >", fc_thres, "and adjusted p-value <", pval_thres),
    x = "log2 Fold Change",
    y = "-log10(p-value)",
    color = "Differential Expression"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12)
  )


### Liver High-DMSO ####
#plotting df
volc_Li.HD <- data.frame(
  logFC = Li.res_HD$log2FoldChange,
  negLogPval = -log10(Li.res_HD$padj),
  adj.P.Val = Li.res_HD$padj,
  ID = rownames(Li.res_HD))

#Column to catagorise genes, if above fc threshold labelled upregulated (greater than 1), if below threshold (-1), downregulated
volc_Li.HD$category <- ifelse(Li.res_HD$padj <= pval_thres,
                              ifelse(volc_Li.HD$logFC >= fc_thres, "Upregulated", 
                                     ifelse(volc_Li.HD$logFC <= -fc_thres, "Downregulated", "Not Significant")),
                              "Not Significant")

# Create the volcano plot using ggplot2
ggplot(volc_Li.HD, aes(x = logFC, y = negLogPval, color = category)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  geom_vline(xintercept = c(-fc_thres, fc_thres), linetype = "dashed") +
  geom_hline(yintercept = -log10(pval_thres), linetype = "dashed") +
  labs(
    title = "Differential Gene Expression in High DEHP vs DMSO control in Liver",
    subtitle = paste("Thresholds: |log2FC| >", fc_thres, "and adjusted p-value <", pval_thres),
    x = "log2 Fold Change",
    y = "-log10(p-value)",
    color = "Differential Expression"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12)
  )

### Liver Low-DMSO ####
#plotting df
volc_Li.LD <- data.frame(
  logFC = Li.res_LD$log2FoldChange,
  negLogPval = -log10(Li.res_LD$padj),
  adj.P.Val = Li.res_LD$padj,
  ID = rownames(Li.res_LD))

#Column to catagorise genes, if above fc threshold labelled upregulated (greater than 1), if below threshold (-1), downregulated
volc_Li.LD$category <- ifelse(Li.res_LD$padj <= pval_thres,
                              ifelse(volc_Li.LD$logFC >= fc_thres, "Upregulated", 
                                     ifelse(volc_Li.LD$logFC <= -fc_thres, "Downregulated", "Not Significant")),
                              "Not Significant")

# Create the volcano plot using ggplot2
ggplot(volc_Li.LD, aes(x = logFC, y = negLogPval, color = category)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  geom_vline(xintercept = c(-fc_thres, fc_thres), linetype = "dashed") +
  geom_hline(yintercept = -log10(pval_thres), linetype = "dashed") +
  labs(
    title = "Differential Gene Expression in Low DEHP vs DMSO control in Liver",
    subtitle = paste("Thresholds: |log2FC| >", fc_thres, "and adjusted p-value <", pval_thres),
    x = "log2 Fold Change",
    y = "-log10(p-value)",
    color = "Differential Expression"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12)
  )

### Liver MP-MPD ####
#plotting df
volc_Li.MP.MPD <- data.frame(
  logFC = Li.res_MP.MPD$log2FoldChange,
  negLogPval = -log10(Li.res_MP.MPD$padj),
  adj.P.Val = Li.res_MP.MPD$padj,
  ID = rownames(Li.res_MP.MPD))

#Column to catagorise genes, if above fc threshold labelled upregulated (greater than 1), if below threshold (-1), downregulated
volc_Li.MP.MPD$category <- ifelse(Li.res_MP.MPD$padj <= pval_thres,
                              ifelse(volc_Li.MP.MPD$logFC >= fc_thres, "Upregulated", 
                                     ifelse(volc_Li.MP.MPD$logFC <= -fc_thres, "Downregulated", "Not Significant")),
                              "Not Significant")

# Create the volcano plot using ggplot2
ggplot(volc_Li.MP.MPD, aes(x = logFC, y = negLogPval, color = category)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  geom_vline(xintercept = c(-fc_thres, fc_thres), linetype = "dashed") +
  geom_hline(yintercept = -log10(pval_thres), linetype = "dashed") +
  labs(
    title = "Differential Gene Expression in Microplastics control vs Microplastics+DEHP in Liver",
    subtitle = paste("Thresholds: |log2FC| >", fc_thres, "and adjusted p-value <", pval_thres),
    x = "log2 Fold Change",
    y = "-log10(p-value)",
    color = "Differential Expression"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12)
  )

### Liver MP-DMSO ####
#plotting df
volc_Li.MP.DMSO <- data.frame(
  logFC = Li.res_MP.DMSO$log2FoldChange,
  negLogPval = -log10(Li.res_MP.DMSO$padj),
  adj.P.Val = Li.res_MP.DMSO$padj,
  ID = rownames(Li.res_MP.DMSO))

#Column to catagorise genes, if above fc threshold labelled upregulated (greater than 1), if below threshold (-1), downregulated
volc_Li.MP.DMSO$category <- ifelse(Li.res_MP.DMSO$padj <= pval_thres,
                                  ifelse(volc_Li.MP.DMSO$logFC >= fc_thres, "Upregulated", 
                                         ifelse(volc_Li.MP.DMSO$logFC <= -fc_thres, "Downregulated", "Not Significant")),
                                  "Not Significant")

# Create the volcano plot using ggplot2
ggplot(volc_Li.MP.DMSO, aes(x = logFC, y = negLogPval, color = category)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  geom_vline(xintercept = c(-fc_thres, fc_thres), linetype = "dashed") +
  geom_hline(yintercept = -log10(pval_thres), linetype = "dashed") +
  labs(
    title = "Differential Gene Expression in Microplastics control vs DMSO control in Liver",
    subtitle = paste("Thresholds: |log2FC| >", fc_thres, "and adjusted p-value <", pval_thres),
    x = "log2 Fold Change",
    y = "-log10(p-value)",
    color = "Differential Expression"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12)
  )

### Liver MPD-DMSO ####
#plotting df
volc_Li.MPD.DMSO <- data.frame(
  logFC = Li.res_MPD.DMSO$log2FoldChange,
  negLogPval = -log10(Li.res_MPD.DMSO$padj),
  adj.P.Val = Li.res_MPD.DMSO$padj,
  ID = rownames(Li.res_MPD.DMSO))

#Column to catagorise genes, if above fc threshold labelled upregulated (greater than 1), if below threshold (-1), downregulated
volc_Li.MPD.DMSO$category <- ifelse(Li.res_MPD.DMSO$padj <= pval_thres,
                                   ifelse(volc_Li.MPD.DMSO$logFC >= fc_thres, "Upregulated", 
                                          ifelse(volc_Li.MPD.DMSO$logFC <= -fc_thres, "Downregulated", "Not Significant")),
                                   "Not Significant")

# Create the volcano plot using ggplot2
ggplot(volc_Li.MPD.DMSO, aes(x = logFC, y = negLogPval, color = category)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  geom_vline(xintercept = c(-fc_thres, fc_thres), linetype = "dashed") +
  geom_hline(yintercept = -log10(pval_thres), linetype = "dashed") +
  labs(
    title = "Differential Gene Expression in Microplastics+DEHP vs DMSO control in Liver",
    subtitle = paste("Thresholds: |log2FC| >", fc_thres, "and adjusted p-value <", pval_thres),
    x = "log2 Fold Change",
    y = "-log10(p-value)",
    color = "Differential Expression"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12)
  )

### Liver F1 ####
#plotting df
volc_Li.F1 <- data.frame(
  logFC = Li.res_F1$log2FoldChange,
  negLogPval = -log10(Li.res_F1$padj),
  adj.P.Val = Li.res_F1$padj,
  ID = rownames(Li.res_F1))

#Column to catagorise genes, if above fc threshold labelled upregulated (greater than 1), if below threshold (-1), downregulated
volc_Li.F1$category <- ifelse(Li.res_F1$padj <= pval_thres,
                                    ifelse(volc_Li.F1$logFC >= fc_thres, "Upregulated", 
                                           ifelse(volc_Li.F1$logFC <= -fc_thres, "Downregulated", "Not Significant")),
                                    "Not Significant")

# Create the volcano plot using ggplot2
ggplot(volc_Li.F1, aes(x = logFC, y = negLogPval, color = category)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  geom_vline(xintercept = c(-fc_thres, fc_thres), linetype = "dashed") +
  geom_hline(yintercept = -log10(pval_thres), linetype = "dashed") +
  labs(
    title = "Differential Gene Expression in Liver F1 High vs DMSO control",
    subtitle = paste("Thresholds: |log2FC| >", fc_thres, "and adjusted p-value <", pval_thres),
    x = "log2 Fold Change",
    y = "-log10(p-value)",
    color = "Differential Expression"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12)
  )
