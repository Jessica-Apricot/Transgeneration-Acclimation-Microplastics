## Goseq
#BiocManager::install("goseq")
#BiocManager::install("org.Dr.eg.db")
library(goseq)
library(dplyr)
library(org.Dr.eg.db)
#For now im going to use danRer6 as it has both gene length and go annotations, will revisit if need danRer11

#We will do Brain F1, Gonad HLD, Gonad Microplastics, Gonad F1 and Liver HLD

gene_info <- read.delim("gene_info.tsv", header = TRUE, sep = "\t")

#### Brain F1 Go analysis #####

length_vector <- gene_info$Length
names(length_vector) <- gene_info$GeneID

BrainF1_genes <- ifelse(Br_res_F1$padj < 0.05, 1, 0)
names(BrainF1_genes) <- rownames(Br_res_F1)

BrF1_valid_genes <- intersect(names(BrainF1_genes), names(length_vector))

BrainF1_genes <- BrainF1_genes[BrF1_valid_genes]
BrF1_length_vector <- length_vector[BrF1_valid_genes]


#Calculate gene weights
BrF1_pwf <- nullp(BrainF1_genes, bias.data = BrF1_length_vector)
par(mfrow=c(1,2))
hist(BrF1_pwf$bias.data, 30)
hist(BrF1_pwf$pwf, 30)

Br_F1_GO <- goseq(BrF1_pwf, "danRer6", "geneSymbol")
Br_F1_GO_padj <- p.adjust(Br_F1_GO$over_represented_pvalue, method = "fdr")
sum(Br_F1_GO_padj < 0.05)
# No significant GO terms, Not enough genes?
View(Br_F1_GO_padj)


#### Gonad HLD #####

##### Gonad High-Low #####

GoHL_genes <- ifelse(Go_res_HL$padj < 0.05, 1, 0)
names(GoHL_genes) <- rownames(Go_res_HL)


GoHL_valid_genes <- intersect(names(GoHL_genes), names(length_vector))

GoHL_genes <- GoHL_genes[GoHL_valid_genes]
GOHL_length_vector <- length_vector[GoHL_valid_genes]

#Calculate gene weight values
GoHL_pwf <- nullp(GoHL_genes, bias.data = GOHL_length_vector)
hist(GoHL_pwf$bias.data, 30)
hist(GoHL_pwf$pwf, 30)

GoHL_GO <- goseq(GoHL_pwf, "danRer6", "geneSymbol")
GoHL_GO_padj <- p.adjust(GoHL_GO$over_represented_pvalue, method = "fdr")
sum(GoHL_GO_padj < 0.05)
#225 significant GO terms

GoHL_GO_sig <- GoHL_GO$category[GoHL_GO_padj < 0.05]
length(GoHL_GO_sig)
head(GoHL_GO_sig)
View(GoHL_GO_sig)

GoHL_GO_df <- data.frame(GO_ID = GoHL_GO$category[GoHL_GO_padj < 0.05],
                         padj = GoHL_GO_padj[GoHL_GO_padj < 0.05])
write.table(GoHL_GO_df, file = "GoHL_GO_sig_with_padj.txt", quote = FALSE, row.names = FALSE, sep = "\t")


#First 5 GO Terms
#GO:0045202 - Synapse
#GO:0030054 - Cell junction
#GO:0099537 - Trans synaptic signalling



##### Gonad High-DMSO #####

GoHD_genes <- ifelse(Go_res_HD$padj < 0.05, 1, 0)
names(GoHD_genes) <- rownames(Go_res_HD)

GoHD_valid_genes <- intersect(names(GoHD_genes), names(length_vector))

GoHD_genes <- GoHD_genes[GoHD_valid_genes]
GoHD_length_vector <- length_vector[GoHD_valid_genes]

#Calculate gene weight values
GoHD_pwf <- nullp(GoHD_genes, bias.data = GoHD_length_vector)
hist(GoHD_pwf$bias.data, 30)
hist(GoHD_pwf$pwf, 30)

GoHD_GO <- goseq(GoHD_pwf, "danRer6", "geneSymbol")
GoHD_GO_padj <- p.adjust(GoHD_GO$over_represented_pvalue, method = "fdr")
sum(GoHD_GO_padj < 0.05)

GoHD_GO_sig <- GoHD_GO$category[GoHD_GO_padj < 0.05]
length(GoHD_GO_sig)
head(GoHD_GO_sig)

#278 significant GO terms


#### Gonad Microplastics #####
##### Gonad MP vs MPD ######

GoMP_MPD_genes <- ifelse(Go_res_MP_MPD$padj < 0.05, 1, 0)
names(GoMP_MPD_genes) <- rownames(Go_res_MP_MPD)

GoMP_MPD_valid_genes <- intersect(names(GoMP_MPD_genes), names(length_vector))

GoMP_MPD_genes <- GoMP_MPD_genes[GoMP_MPD_valid_genes]
GoMP_MPD_length_vector <- length_vector[GoMP_MPD_valid_genes]

#Calculate gene weight values
GoMP_MPD_pwf <- nullp(GoMP_MPD_genes, bias.data = GoMP_MPD_length_vector)
hist(GoMP_MPD_pwf$bias.data, 30)
hist(GoMP_MPD_pwf$pwf, 30)

GoMP_MPD_GO <- goseq(GoMP_MPD_pwf, "danRer6", "geneSymbol")
GoMP_MPD_GO_padj <- p.adjust(GoHD_GO$over_represented_pvalue, method = "fdr")
sum(GoMP_MPD_GO_padj < 0.05)
#278 significant GO terms

GoMP_MPD_GO_sig <- GoMP_MPD_GO$category[GoMP_MPD_GO_padj < 0.05]
length(GoMP_MPD_GO_sig)
head(GoMP_MPD_GO_sig)
head(GoMP_MPD_GO)

#First Five significant genes
#GO:0045202 - synapse
#GO:0030054 - cell junction
#GO:0007267 - cell-cell signaling
#GO:0099537 - trans-synaptic signaling 
#GO:0099536 - synaptic signaling

##### Gonad MP vs DMSO #####
GoMP_DMSO_genes <- ifelse(Go_res_MP_DMSO$padj < 0.05, 1, 0)
names(GoMP_DMSO_genes) <- rownames(Go_res_MP_DMSO)

GoMP_DMSO_valid_genes <- intersect(names(GoMP_DMSO_genes), names(length_vector))

GoMP_DMSO_genes <- GoMP_DMSO_genes[GoMP_DMSO_valid_genes]
GoMP_DMSO_length_vector <- length_vector[GoMP_DMSO_valid_genes]

#Calculate gene weight values
GoMP_DMSO_pwf <- nullp(GoMP_DMSO_genes, bias.data = GoMP_DMSO_length_vector)
hist(GoMP_DMSO_pwf$bias.data, 30)
hist(GoMP_DMSO_pwf$pwf, 30)

GoMP_DMSO_GO <- goseq(GoMP_DMSO_pwf, "danRer6", "geneSymbol")
GoMP_DMSO_GO_padj <- p.adjust(GoMP_DMSO_GO$over_represented_pvalue, method = "fdr")
sum(GoMP_DMSO_GO_padj < 0.05)
#57 significant GO terms

GoMP_DMSO_GO_sig <- GoMP_DMSO_GO$category[GoMP_DMSO_GO_padj < 0.05]
length(GoMP_DMSO_GO_sig)
head(GoMP_DMSO_GO_sig)
head(GoMP_DMSO_GO)

#First five significant GO Terms
#GO:0045202 - synapse
#GO:0030054 - cell junction
#GO:0023052 - signaling       
#GO:0016020 - membrane
#GO:0007154 - cell communication

#### Gonad F1 #####

GoF1_genes <- ifelse(Go_res_F1$padj < 0.05, 1, 0)
names(GoF1_genes) <- rownames(Go_res_F1)

GoF1_valid_genes <- intersect(names(GoF1_genes), names(length_vector))

GoF1_genes <- GoF1_genes[GoF1_valid_genes]
GoF1_length_vector <- length_vector[GoF1_valid_genes]

#Calculate gene weight values
GoF1_pwf <- nullp(GoF1_genes, bias.data = GoF1_length_vector)
hist(GoF1_pwf$bias.data, 30)
hist(GoF1_pwf$pwf, 30)

GoF1_GO <- goseq(GoF1_pwf, "danRer6", "geneSymbol")
GoF1_GO_padj <- p.adjust(GoF1_GO$over_represented_pvalue, method = "fdr")
sum(GoF1_GO_padj < 0.05)
#0 significant GO Terms

#### Liver HLD ####
###### Liver Low-DMSO ####

LiLD_genes <- ifelse(Li_res_LD$padj < 0.05, 1, 0)
names(LiLD_genes) <- rownames(Li_res_LD)

LiLD_valid_genes <- intersect(names(LiLD_genes), names(length_vector))

LiLD_genes <- LiLD_genes[LiLD_valid_genes]
LiLD_length_vector <- length_vector[LiLD_valid_genes]

#Calculate gene weight values
LiLD_pwf <- nullp(LiLD_genes, bias.data = LiLD_length_vector)
hist(LiLD_pwf$bias.data, 30)
hist(LiLD_pwf$pwf, 30)

LiLD_GO <- goseq(LiLD_pwf, "danRer6", "geneSymbol")
LiLD_GO_padj <- p.adjust(LiLD_GO$over_represented_pvalue, method = "fdr")
sum(LiLD_GO_padj < 0.05)
#0 significant GO terms




# packages
#install.packages("tidyverse")
#install.packages("forcats")
library(tidyverse)
library(forcats)


Up_Fun_GoHD <- read_csv("Up_Function_Go_HD.csv")   

# --- ORDER TERMS WITHIN EACH PANEL BY SIGNIFICANCE ---------------------------


# Order terms by significance

Up_Fun_GoHD_ord <- Up_Fun_GoHD %>%
  filter(!is.na(padj)) %>%  
  mutate(term = fct_reorder(term, padj, .desc = TRUE),
         neg_log10_padj = -log10(padj))


ggplot(Up_Fun_GoHD_ord, aes(neg_log10_padj, term, color = category)) +
  geom_point(size = 3) +
  scale_color_manual(values = c(Function = "orange")) +
  labs(x = expression(-log[10]~"(Adjusted P-value)"), y = "GO Terms", color = "Category") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.title.y = element_text(margin = margin(r = 10))
  )





# --- ORDER TERMS WITHIN EACH PANEL BY SIGNIFICANCE ---------------------------
# smaller padj = more significant; we’ll show the most significant near x=0
Up_Fun_Pro_GoHD <- read_csv("Up_Fun_Pro_GoHD.csv")  
Up_Fun_Pro_GoHD_ord <- Up_Fun_Pro_GoHD %>%
  group_by(category) %>%
  mutate(term = fct_reorder(term, padj, .desc = TRUE), neg_log10_padj = -log10(padj)) %>%  # order within panel
  ungroup()

# --- PLOT --------------------------------------------------------------------
ggplot(Up_Fun_Pro_GoHD_ord, aes(neg_log10_padj, term, color = category)) +
  geom_point(size = 3) +
  facet_wrap(~ category, ncol = 1, scales = "free_y") +
  scale_color_manual(
    values = c(Process = "pink",
               Function = "orange")  
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  labs(x = expression(-log[10]~"(Adjusted P-value)"), y = "GO Terms", color = "Category") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.title.y = element_text(margin = margin(r = 10))
  )






Up_Fun_MP_DMSO <- read_csv("Up_Fuction_MP_DMSO.csv")   

# --- ORDER TERMS WITHIN EACH PANEL BY SIGNIFICANCE ---------------------------


# Order terms by significance

Up_Fun_MP_DMSO_ord <- Up_Fun_MP_DMSO %>%
  filter(!is.na(padj)) %>%  
  mutate(term = fct_reorder(term, padj, .desc = TRUE),
         neg_log10_padj = -log10(padj))


ggplot(Up_Fun_MP_DMSO_ord, aes(neg_log10_padj, term, color = category)) +
  geom_point(size = 3) +
  scale_color_manual(values = c(Function = "orange")) +
  labs(x = expression(-log[10]~"(Adjusted P-value)"), y = "GO Terms", color = "Category") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.title.y = element_text(margin = margin(r = 10))
  )


# --- ORDER TERMS WITHIN EACH PANEL BY SIGNIFICANCE ---------------------------
# smaller padj = more significant; we’ll show the most significant near x=0
Up_Fun_Pro_Mp_DMSO <- read_csv("Up_Fun_Pro_Mp_DMSO.csv")  
Up_Fun_Pro_Mp_DMSO_ord <- Up_Fun_Pro_Mp_DMSO %>%
  group_by(category) %>%
  mutate(term = fct_reorder(term, padj, .desc = TRUE), neg_log10_padj = -log10(padj)) %>%  # order within panel
  ungroup()

# --- PLOT --------------------------------------------------------------------
ggplot(Up_Fun_Pro_Mp_DMSO_ord, aes(neg_log10_padj, term, color = category)) +
  geom_point(size = 3) +
  facet_wrap(~ category, ncol = 1, scales = "free_y") +
  scale_color_manual(
    values = c(Process = "pink",
               Function = "orange")  
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  labs(x = expression(-log[10]~"(Adjusted P-value)"), y = "GO Terms", color = "Category") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.title.y = element_text(margin = margin(r = 10))
  )

