#Volcano plots
#install.packages(c("systemfonts", "ragg","ggrepel", "plotly", "ggiraph"))
library(ggplot2)
library(ggiraph)
library(plotly)
library(ggrepel)
library(systemfonts)
library(ragg)

#install.packages(c("ggplot2", "grDevices", "Cairo"))

#Define significance thresholds
pval_thres <- 0.05
fc_thres <- 1 #log2 fold

### Brain High-Low ####
#plotting df
volc_Br_HL <- data.frame(
  logFC = Br_res_HL$log2FoldChange,
  negLogPval = -log10(Br_res_HL$padj),
  adj.P.Val = Br_res_HL$padj,
  ID = rownames(Br_res_HL))

#Column to catagorise genes, if above fc threshold labelled upregulated (greater than 1), if below threshold (-1), downregulated
volc_Br_HL$category <- ifelse(Br_res_HL$padj <= pval_thres,
                              ifelse(volc_Br_HL$logFC >= fc_thres, "Upregulated", 
                                     ifelse(volc_Br_HL$logFC <= -fc_thres, "Downregulated", "Not Significant")),
                              "Not Significant")

# Create the volcano plot using ggplot2
ggplot(volc_Br_HL, aes(x = logFC, y = negLogPval, color = category)) +
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
volc_Br_HD <- data.frame(
  logFC = Br_res_HD$log2FoldChange,
  negLogPval = -log10(Br_res_HD$padj),
  adj.P.Val = Br_res_HD$padj,
  ID = rownames(Br_res_HD))

#Column to catagorise genes, if above fc threshold labelled upregulated (greater than 1), if below threshold (-1), downregulated
volc_Br_HD$category <- ifelse(Br_res_HD$padj <= pval_thres,
                              ifelse(volc_Br_HD$logFC >= fc_thres, "Upregulated", 
                                     ifelse(volc_Br_HD$logFC <= -fc_thres, "Downregulated", "Not Significant")),
                              "Not Significant")

# Create the volcano plot using ggplot2
ggplot(volc_Br_HD, aes(x = logFC, y = negLogPval, color = category)) +
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
volc_Br_LD <- data.frame(
  logFC = Br_res_LD$log2FoldChange,
  negLogPval = -log10(Br_res_LD$padj),
  adj.P.Val = Br_res_LD$padj,
  ID = rownames(Br_res_LD))

#Column to catagorise genes, if above fc threshold labelled upregulated (greater than 1), if below threshold (-1), downregulated
volc_Br_LD$category <- ifelse(Br_res_LD$padj <= pval_thres,
                              ifelse(volc_Br_LD$logFC >= fc_thres, "Upregulated", 
                                     ifelse(volc_Br_LD$logFC <= -fc_thres, "Downregulated", "Not Significant")),
                              "Not Significant")

# Create the volcano plot using ggplot2
ggplot(volc_Br_LD, aes(x = logFC, y = negLogPval, color = category)) +
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
volc_Br_MP_MPD <- data.frame(
  logFC = Br_res_MP_MPD$log2FoldChange,
  negLogPval = -log10(Br_res_MP_MPD$padj),
  adj.P.Val = Br_res_MP_MPD$padj,
  ID = rownames(Br_res_MP_MPD))

volc_Br_MP_MPD$category <- ifelse(Br_res_MP_MPD$padj <= pval_thres,
                              ifelse(volc_Br_MP_MPD$logFC >= fc_thres, "Upregulated", 
                                     ifelse(volc_Br_MP_MPD$logFC <= -fc_thres, "Downregulated", "Not Significant")),
                              "Not Significant")

# Create the volcano plot using ggplot2
ggplot(volc_Br_MP_MPD, aes(x = logFC, y = negLogPval, color = category)) +
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
    plot.title = element_text(hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 12) )


### Brain MP-DMSO ####
#plotting df
volc_Br_MP_DMSO <- data.frame(
  logFC = Br_res_MP_DMSO$log2FoldChange,
  negLogPval = -log10(Br_res_MP_DMSO$padj),
  adj.P.Val = Br_res_MP_DMSO$padj,
  ID = rownames(Br_res_MP_DMSO))

#Column to catagorise genes, if above fc threshold labelled upregulated (greater than 1), if below threshold (-1), downregulated
volc_Br_MP_DMSO$category <- ifelse(Br_res_MP_DMSO$padj <= pval_thres,
                                  ifelse(volc_Br_MP_DMSO$logFC >= fc_thres, "Upregulated", 
                                         ifelse(volc_Br_MP_DMSO$logFC <= -fc_thres, "Downregulated", "Not Significant")),
                                  "Not Significant")

#Labelling
volc_Br_MP_DMSO$label <- ifelse(
  volc_Br_MP_DMSO$category != "Not Significant" &
    abs(volc_Br_MP_DMSO$logFC) > 1, volc_Br_MP_DMSO$ID,
  NA)


# Create the volcano plot using ggplot2
ggplot(volc_Br_MP_DMSO, aes(x = logFC, y = negLogPval, color = category)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  geom_vline(xintercept = c(-fc_thres, fc_thres), linetype = "dashed") +
  geom_hline(yintercept = -log10(pval_thres), linetype = "dashed")+
  geom_text_repel(
    data = subset(volc_Br_MP_DMSO, !is.na(label)),
    aes(x = logFC, y = negLogPval, label = label),  # re-specify x and y
    inherit.aes = FALSE,  # avoid pulling global aesthetics like color
    size = 3,
    max.overlaps = 20,
    box.padding = 0.3,
    point.padding = 0.2,
    segment.color = 'grey50'
  ) +
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
    plot.title = element_text(hjust = 0.5, size = 13),
    plot.subtitle = element_text(hjust = 0.5, size = 12) )



### Brain MPD-DMSO ####
#One of the sig genes are overlapping
#plotting df
volc_Br_MPD_DMSO <- data.frame(
  logFC = Br_res_MPD_DMSO$log2FoldChange,
  negLogPval = -log10(Br_res_MPD_DMSO$padj),
  adj.P.Val = Br_res_MPD_DMSO$padj,
  ID = rownames(Br_res_MPD_DMSO))

#Column to catagorise genes, if above fc threshold labelled upregulated (greater than 1), if below threshold (-1), downregulated
volc_Br_MPD_DMSO$category <- ifelse(Br_res_MPD_DMSO$padj <= pval_thres,
                                   ifelse(volc_Br_MPD_DMSO$logFC >= fc_thres, "Upregulated", 
                                          ifelse(volc_Br_MPD_DMSO$logFC <= -fc_thres, "Downregulated", "Not Significant")),
                                   "Not Significant")

# Create the volcano plot using ggplot2
ggplot(volc_Br_MPD_DMSO, aes(x = logFC, y = negLogPval, color = category)) +
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
    plot.title = element_text(hjust = 0.5, size = 13),
    plot.subtitle = element_text(hjust = 0.5, size = 12) )


### Brain F1 ####
#plotting df
volc_Br_F1 <- data.frame(
  logFC = Br_res_F1$log2FoldChange,
  negLogPval = -log10(Br_res_F1$padj),
  adj.P.Val = Br_res_F1$padj,
  ID = rownames(Br_res_F1))


volc_Br_F1$category <- ifelse(Br_res_F1$padj <= pval_thres,
                              ifelse(volc_Br_F1$logFC >= fc_thres, "Upregulated", 
                                     ifelse(volc_Br_F1$logFC <= -fc_thres, "Downregulated", "Not Significant")),
                              "Not Significant")


top_n <- 10
BrF1_sig <- volc_Br_F1[
  abs(volc_Br_F1$logFC) > fc_thres & volc_Br_F1$adj.P.Val < pval_thres,
]
BrF1_top <- BrF1_sig[order(-abs(BrF1_sig$logFC)), ][1:min(top_n, nrow(BrF1_sig)), ]

# Volcano plot
ggplot(volc_Br_F1, aes(x = logFC, y = negLogPval, color = category)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_text_repel(data = BrF1_top, aes(label = ID), size = 3, max.overlaps = Inf) +
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
  scale_x_continuous(
    breaks = seq(-6, 6, by = 1) 
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    plot.title = element_blank(),
    plot.subtitle = element_blank()
  )


### Gonad High-Low ####
#plotting df
volc_Go_HL <- data.frame(
  logFC = Go_res_HL$log2FoldChange,
  negLogPval = -log10(Go_res_HL$padj),
  adj.P.Val = Go_res_HL$padj,
  ID = rownames(Go_res_HL))

#Column to catagorise genes, if above fc threshold labelled upregulated (greater than 1), if below threshold (-1), downregulated
volc_Go_HL$category <- ifelse(Go_res_HL$padj <= pval_thres,
                              ifelse(volc_Go_HL$logFC >= fc_thres, "Upregulated", 
                                     ifelse(volc_Go_HL$logFC <= -fc_thres, "Downregulated", "Not Significant")),
                              "Not Significant")

# Create the volcano plot using ggplot2
ggplot(volc_Go_HL, aes(x = logFC, y = negLogPval, color = category)) +
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
    plot.title = element_text(hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 12)
  )

### Gonad High-DMSO ####
#plotting df
volc_Go_HD <- data.frame(
  logFC = Go_res_HD$log2FoldChange,
  negLogPval = -log10(Go_res_HD$padj),
  adj.P.Val = Go_res_HD$padj,
  ID = rownames(Go_res_HD))


#Column to catagorise genes, if above fc threshold labelled upregulated (greater than 1), if below threshold (-1), downregulated
volc_Go_HD$category <- ifelse(Go_res_HD$padj <= pval_thres,
                              ifelse(volc_Go_HD$logFC >= fc_thres, "Upregulated", 
                                     ifelse(volc_Go_HD$logFC <= -fc_thres, "Downregulated", "Not Significant")),
                              "Not Significant")

# Create the volcano plot using ggplot2
ggplot(volc_Go_HD, aes(x = logFC, y = negLogPval, color = category)) +
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
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    plot.title = element_blank(),
    plot.subtitle = element_blank()
  )

### Gonad Low-DMSO ####
#plotting df
volc_Go_LD <- data.frame(
  logFC = Go_res_LD$log2FoldChange,
  negLogPval = -log10(Go_res_LD$padj),
  adj.P.Val = Go_res_LD$padj,
  ID = rownames(Go_res_LD))

#Column to catagorise genes, if above fc threshold labelled upregulated (greater than 1), if below threshold (-1), downregulated
volc_Go_LD$category <- ifelse(Go_res_LD$padj <= pval_thres,
                              ifelse(volc_Go_LD$logFC >= fc_thres, "Upregulated", 
                                     ifelse(volc_Go_LD$logFC <= -fc_thres, "Downregulated", "Not Significant")),
                              "Not Significant")

# Create the volcano plot using ggplot2
ggplot(volc_Go_LD, aes(x = logFC, y = negLogPval, color = category)) +
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
volc_Go_MP_MPD <- data.frame(
  logFC = Go_res_MP_MPD$log2FoldChange,
  negLogPval = -log10(Go_res_MP_MPD$padj),
  adj.P.Val = Go_res_MP_MPD$padj,
  ID = rownames(Go_res_MP_MPD))

#Column to catagorise genes, if above fc threshold labelled upregulated (greater than 1), if below threshold (-1), downregulated
volc_Go_MP_MPD$category <- ifelse(Go_res_MP_MPD$padj <= pval_thres,
                              ifelse(volc_Go_MP_MPD$logFC >= fc_thres, "Upregulated", 
                                     ifelse(volc_Go_MP_MPD$logFC <= -fc_thres, "Downregulated", "Not Significant")),
                              "Not Significant")

# Create the volcano plot using ggplot2
ggplot(volc_Go_LD, aes(x = logFC, y = negLogPval, color = category)) +
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
volc_Go_MP_DMSO <- data.frame(
  logFC = Go_res_MP_DMSO$log2FoldChange,
  negLogPval = -log10(Go_res_MP_DMSO$padj),
  adj.P.Val = Go_res_MP_DMSO$padj,
  ID = rownames(Go_res_MP_DMSO))

#Column to catagorise genes, if above fc threshold labelled upregulated (greater than 1), if below threshold (-1), downregulated
volc_Go_MP_DMSO$category <- ifelse(Go_res_MP_DMSO$padj <= pval_thres,
                                  ifelse(volc_Go_MP_DMSO$logFC >= fc_thres, "Upregulated", 
                                         ifelse(volc_Go_MP_DMSO$logFC <= -fc_thres, "Downregulated", "Not Significant")),
                                  "Not Significant")

# Create the volcano plot using ggplot2
ggplot(volc_Go_MP_DMSO, aes(x = logFC, y = negLogPval, color = category)) +
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
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    plot.title = element_blank(),
    plot.subtitle = element_blank()
  )


### Gonad MPD - DMSO ####
volc_Go_MPD_DMSO <- data.frame(
  logFC = Go_res_MPD_DMSO$log2FoldChange,
  negLogPval = -log10(Go_res_MPD_DMSO$padj),
  adj.P.Val = Go_res_MPD_DMSO$padj,
  ID = rownames(Go_res_MPD_DMSO))

#Column to catagorise genes, if above fc threshold labelled upregulated (greater than 1), if below threshold (-1), downregulated
volc_Go_MPD_DMSO$category <- ifelse(Go_res_MPD_DMSO$padj <= pval_thres,
                                   ifelse(volc_Go_MPD_DMSO$logFC >= fc_thres, "Upregulated", 
                                          ifelse(volc_Go_MPD_DMSO$logFC <= -fc_thres, "Downregulated", "Not Significant")),
                                   "Not Significant")

# Create the volcano plot using ggplot2
ggplot(volc_Go_MPD_DMSO, aes(x = logFC, y = negLogPval, color = category)) +
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
volc_Go_F1 <- data.frame(
  logFC = Go_res_F1$log2FoldChange,
  negLogPval = -log10(Go_res_F1$padj),
  adj.P.Val = Go_res_F1$padj,
  ID = rownames(Go_res_F1))


#Column to catagorise genes, if above fc threshold labelled upregulated (greater than 1), if below threshold (-1), downregulated
volc_Go_F1$category <- ifelse(Go_res_F1$padj <= pval_thres,
                                    ifelse(volc_Go_F1$logFC >= fc_thres, "Upregulated", 
                                           ifelse(volc_Go_F1$logFC <= -fc_thres, "Downregulated", "Not Significant")),
                                    "Not Significant")

GoF1_sig <- volc_Go_F1[
  abs(volc_Go_F1$logFC) > fc_thres & volc_Go_F1$adj.P.Val < pval_thres,
]
GoF1_top <- GoF1_sig[order(-abs(GoF1_sig$logFC)), ][1:min(top_n, nrow(GoF1_sig)), ]


# Create the volcano plot using ggplot2
ggplot(volc_Go_F1, aes(x = logFC, y = negLogPval, color = category)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_text_repel(data = GoF1_top, aes(label = ID), size = 3, max.overlaps = Inf)+
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
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    plot.title = element_blank(),
    plot.subtitle = element_blank()
  )

### Liver High-Low ####
#plotting df
volc_Li_HL <- data.frame(
  logFC = Li_res_HL$log2FoldChange,
  negLogPval = -log10(Li_res_HL$padj),
  adj.P.Val = Li_res_HL$padj,
  ID = rownames(Li_res_HL))

#Column to catagorise genes, if above fc threshold labelled upregulated (greater than 1), if below threshold (-1), downregulated
volc_Li_HL$category <- ifelse(Li_res_HL$padj <= pval_thres,
                              ifelse(volc_Li_HL$logFC >= fc_thres, "Upregulated", 
                                     ifelse(volc_Li_HL$logFC <= -fc_thres, "Downregulated", "Not Significant")),
                              "Not Significant")

# Create the volcano plot using ggplot2
ggplot(volc_Li_HL, aes(x = logFC, y = negLogPval, color = category)) +
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
volc_Li_HD <- data.frame(
  logFC = Li_res_HD$log2FoldChange,
  negLogPval = -log10(Li_res_HD$padj),
  adj.P.Val = Li_res_HD$padj,
  ID = rownames(Li_res_HD))

#Column to catagorise genes, if above fc threshold labelled upregulated (greater than 1), if below threshold (-1), downregulated
volc_Li_HD$category <- ifelse(Li_res_HD$padj <= pval_thres,
                              ifelse(volc_Li_HD$logFC >= fc_thres, "Upregulated", 
                                     ifelse(volc_Li_HD$logFC <= -fc_thres, "Downregulated", "Not Significant")),
                              "Not Significant")

# Create the volcano plot using ggplot2
ggplot(volc_Li_HD, aes(x = logFC, y = negLogPval, color = category)) +
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
volc_Li_LD <- data.frame(
  logFC = Li_res_LD$log2FoldChange,
  negLogPval = -log10(Li_res_LD$padj),
  adj.P.Val = Li_res_LD$padj,
  ID = rownames(Li_res_LD))

#Column to catagorise genes, if above fc threshold labelled upregulated (greater than 1), if below threshold (-1), downregulated
volc_Li_LD$category <- ifelse(Li_res_LD$padj <= pval_thres,
                              ifelse(volc_Li_LD$logFC >= fc_thres, "Upregulated", 
                                     ifelse(volc_Li_LD$logFC <= -fc_thres, "Downregulated", "Not Significant")),
                              "Not Significant")

# Create the volcano plot using ggplot2
ggplot(volc_Li_LD, aes(x = logFC, y = negLogPval, color = category)) +
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
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    plot.title = element_blank(),
    plot.subtitle = element_blank()
  )

### Liver MP-MPD ####
#plotting df
volc_Li_MP_MPD <- data.frame(
  logFC = Li_res_MP_MPD$log2FoldChange,
  negLogPval = -log10(Li_res_MP_MPD$padj),
  adj.P.Val = Li_res_MP_MPD$padj,
  ID = rownames(Li_res_MP_MPD))

#Column to catagorise genes, if above fc threshold labelled upregulated (greater than 1), if below threshold (-1), downregulated
volc_Li_MP_MPD$category <- ifelse(Li_res_MP_MPD$padj <= pval_thres,
                              ifelse(volc_Li_MP_MPD$logFC >= fc_thres, "Upregulated", 
                                     ifelse(volc_Li_MP_MPD$logFC <= -fc_thres, "Downregulated", "Not Significant")),
                              "Not Significant")

# Create the volcano plot using ggplot2
ggplot(volc_Li_MP_MPD, aes(x = logFC, y = negLogPval, color = category)) +
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
volc_Li_MP_DMSO <- data.frame(
  logFC = Li_res_MP_DMSO$log2FoldChange,
  negLogPval = -log10(Li_res_MP_DMSO$padj),
  adj.P.Val = Li_res_MP_DMSO$padj,
  ID = rownames(Li_res_MP_DMSO))

#Column to catagorise genes, if above fc threshold labelled upregulated (greater than 1), if below threshold (-1), downregulated
volc_Li_MP_DMSO$category <- ifelse(Li_res_MP_DMSO$padj <= pval_thres,
                                  ifelse(volc_Li_MP_DMSO$logFC >= fc_thres, "Upregulated", 
                                         ifelse(volc_Li_MP_DMSO$logFC <= -fc_thres, "Downregulated", "Not Significant")),
                                  "Not Significant")

# Create the volcano plot using ggplot2
ggplot(volc_Li_MP_DMSO, aes(x = logFC, y = negLogPval, color = category)) +
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
volc_Li_MPD_DMSO <- data.frame(
  logFC = Li_res_MPD_DMSO$log2FoldChange,
  negLogPval = -log10(Li_res_MPD_DMSO$padj),
  adj.P.Val = Li_res_MPD_DMSO$padj,
  ID = rownames(Li_res_MPD_DMSO))

#Column to catagorise genes, if above fc threshold labelled upregulated (greater than 1), if below threshold (-1), downregulated
volc_Li_MPD_DMSO$category <- ifelse(Li_res_MPD_DMSO$padj <= pval_thres,
                                   ifelse(volc_Li_MPD_DMSO$logFC >= fc_thres, "Upregulated", 
                                          ifelse(volc_Li_MPD_DMSO$logFC <= -fc_thres, "Downregulated", "Not Significant")),
                                   "Not Significant")

# Create the volcano plot using ggplot2
ggplot(volc_Li_MPD_DMSO, aes(x = logFC, y = negLogPval, color = category)) +
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
volc_Li_F1 <- data.frame(
  logFC = Li_res_F1$log2FoldChange,
  negLogPval = -log10(Li_res_F1$padj),
  adj.P.Val = Li_res_F1$padj,
  ID = rownames(Li_res_F1))



#Column to catagorise genes, if above fc threshold labelled upregulated (greater than 1), if below threshold (-1), downregulated
volc_Li_F1$category <- ifelse(Li_res_F1$padj <= pval_thres,
                                    ifelse(volc_Li_F1$logFC >= fc_thres, "Upregulated", 
                                           ifelse(volc_Li_F1$logFC <= -fc_thres, "Downregulated", "Not Significant")),
                                    "Not Significant")

LiF1_sig <- volc_Li_F1[
  abs(volc_Li_F1$logFC) > fc_thres & volc_Li_F1$adj.P.Val < pval_thres,
]
LiF1_top <- LiF1_sig[order(-abs(LiF1_sig$logFC)), ][1:min(top_n, nrow(LiF1_sig)), ]

# Create the volcano plot using ggplot2
ggplot(volc_Li_F1, aes(x = logFC, y = negLogPval, color = category)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_text_repel(data = LiF1_top, aes(label = ID), size = 3, max.overlaps = Inf) +
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
  ) + xlim(-4, 4) + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    plot.title = element_blank(),
    plot.subtitle = element_blank()
  )
