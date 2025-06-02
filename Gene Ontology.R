## Goseq
#BiocManager::install("goseq")
#BiocManager::install("org.Dr.eg.db")
BiocManager::install("topGO")
library(goseq)
library(dplyr)
library(org.Dr.eg.db)
library(topGO)
#For now im going to use danRer6 as it has both gene length and go annotations, will revisit if need danRer11

#We will do Brain F1, Gonad HLD, Gonad Microplastics, Gonad F1 and Liver HLD

gene_info <- read.delim("gene_info.tsv", header = TRUE, sep = "\t")

#### Brain F1 Go analysis #####

length_vector <- gene_info$Length
names(length_vector) <- gene_info$GeneID

BrainF1_genes <- ifelse(Br.res_F1$padj < 0.05, 1, 0)
names(BrainF1_genes) <- rownames(Br.res_F1)

BrF1_valid_genes <- intersect(names(BrainF1_genes), names(length_vector))

BrainF1_genes <- BrainF1_genes[BrF1_valid_genes]
length_vector <- length_vector[BrF1_valid_genes]


#Calculate gene weights
BrF1_pwf <- nullp(BrainF1_genes, bias.data = length_vector)
par(mfrow=c(1,2))
hist(BrF1_pwf$bias.data, 30)
hist(BrF1_pwf$pwf, 30)

Br_F1_GO <- goseq(BrF1_pwf, "danRer6", "geneSymbol")
Br_F1_GO.padj <- p.adjust(Br_F1_GO$over_represented_pvalue, method = "fdr")
sum(Br_F1_GO.padj < 0.05)
# No significant GO terms, Not enough genes?


#### Gonad HLD #####

##### Gonad High-Low #####

GoHL_genes <- ifelse(Go.res_HL$padj < 0.05, 1, 0)
names(GoHL_genes) <- rownames(Go.res_HL)


GoHL_valid_genes <- intersect(names(GoHL_genes), names(length_vector))

GoHL_genes <- GoHL_genes[GoHL_valid_genes]
length_vector <- length_vector[GoHL_valid_genes]

#Calculate gene weight values
GoHL_pwf <- nullp(GoHL_genes, bias.data = length_vector)
hist(GoHL_pwf$bias.data, 30)
hist(GoHL_pwf$pwf, 30)

GoHL_GO <- goseq(GoHL_pwf, "danRer6", "geneSymbol")
GoHL_GO.padj <- p.adjust(GoHL_GO$over_represented_pvalue, method = "fdr")
sum(GoHL_GO.padj < 0.05)
#257 significant GO terms

GoHL_GO.sig <- GoHL_GO$category[GoHL_GO.padj < 0.05]
length(GoHL_GO.sig)
head(GoHL_GO.sig)
View(GoHL_GO.sig)

showSigOfNodes(GoHL_GO.sig, score(sigGenes), firstSigNodes = 5, useInfo = "all")
#First 5 GO Terms
#GO:0045202 - Synapse
#GO:0030054 - Cell junction
#GO:0007267 - cell - cell signalling
#GO:0099537 - Trans synaptic signalling
#GO:0005886 - plasma membrane

##### Gonad High-DMSO #####

GoHD_genes <- ifelse(Go.res_HD$padj < 0.05, 1, 0)
names(GoHD_genes) <- rownames(Go.res_HD)

GoHD_valid_genes <- intersect(names(GoHD_genes), names(length_vector))

GoHD_genes <- GoHD_genes[GoHD_valid_genes]
length_vector <- length_vector[GoHD_valid_genes]

#Calculate gene weight values
GoHD_pwf <- nullp(GoHD_genes, bias.data = length_vector)
hist(GoHD_pwf$bias.data, 30)
hist(GoHD_pwf$pwf, 30)

GoHD_GO <- goseq(GoHD_pwf, "danRer6", "geneSymbol")
GoHD_GO.padj <- p.adjust(GoHD_GO$over_represented_pvalue, method = "fdr")
sum(GoHD_GO.padj < 0.05)

GoHD_GO.sig <- GoHD_GO$category[GoHD_GO.padj < 0.05]
length(GoHD_GO.sig)
head(GoHD_GO.sig)
#282 significant GO terms

#First 5 significant GO terms
#GO:0045202 - Synapse
#GO:0030054 - Cell Junction
#GO:0043005 - Neuron Projection
#

#### Gonad Microplastics #####
##### Gonad MP vs MPD ######

GoMP.MPD_genes <- ifelse(Go.res_MP.MPD$padj < 0.05, 1, 0)
names(GoMP.MPD_genes) <- rownames(Go.res_MP.MPD)

GoMP.MPD_valid_genes <- intersect(names(GoMP.MPD_genes), names(length_vector))

GoMP.MPD_genes <- GoMP.MPD_genes[GoMP.MPD_valid_genes]
length_vector <- length_vector[GoMP.MPD_valid_genes]

#Calculate gene weight values
GoMP.MPD_pwf <- nullp(GoMP.MPD_genes, bias.data = length_vector)
hist(GoMP.MPD_pwf$bias.data, 30)
hist(GoMP.MPD_pwf$pwf, 30)

GoMP.MPD_GO <- goseq(GoMP.MPD_pwf, "danRer6", "geneSymbol")
GoMP.MPD_GO.padj <- p.adjust(GoHD_GO$over_represented_pvalue, method = "fdr")
sum(GoMP.MPD_GO.padj < 0.05)

GoMP.MPD_GO.sig <- GoMP.MPD_GO$category[GoMP.MPD_GO.padj < 0.05]
length(GoMP.MPD_GO.sig)
head(GoMP.MPD_GO.sig)


##### Gonad MP vs DMSO #####
GoMP.DMSO_genes <- ifelse(Go.res_MP.DMSO$padj < 0.05, 1, 0)
names(GoMP.DMSO_genes) <- rownames(Go.res_MP.DMSO)

GoMP.DMSO_valid_genes <- intersect(names(GoMP.DMSO_genes), names(length_vector))

GoMP.DMSO_genes <- GoMP.DMSO_genes[GoMP.DMSO_valid_genes]
length_vector <- length_vector[GoMP.DMSO_valid_genes]

#Calculate gene weight values
GoMP.DMSO_pwf <- nullp(GoMP.DMSO_genes, bias.data = length_vector)
hist(GoMP.DMSO_pwf$bias.data, 30)
hist(GoMP.DMSO_pwf$pwf, 30)

GoMP.DMSO_GO <- goseq(GoMP.DMSO_pwf, "danRer6", "geneSymbol")
GoMP.DMSO_GO.padj <- p.adjust(GoMP.DMSO_GO$over_represented_pvalue, method = "fdr")
sum(GoMP.DMSO_GO.padj < 0.05)

GoMP.DMSO_GO.sig <- GoMP.DMSO_GO$category[GoMP.DMSO_GO.padj < 0.05]
length(GoMP.DMSO_GO.sig)
head(GoMP.DMSO_GO.sig)

#### Gonad F1 #####

GoF1_genes <- ifelse(Go.res_F1$padj < 0.05, 1, 0)
names(GoF1_genes) <- rownames(Go.res_F1)

GoF1_valid_genes <- intersect(names(GoF1_genes), names(length_vector))

GoF1_genes <- GoF1_genes[GoF1_valid_genes]
length_vector <- length_vector[GoF1_valid_genes]

#Calculate gene weight values
GoF1_pwf <- nullp(GoF1_genes, bias.data = length_vector)
hist(GoF1_pwf$bias.data, 30)
hist(GoF1_pwf$pwf, 30)

GoF1_GO <- goseq(GoF1_pwf, "danRer6", "geneSymbol")
GoF1_GO.padj <- p.adjust(GoF1_GO$over_represented_pvalue, method = "fdr")
sum(GoF1_GO.padj < 0.05)

#### Liver HLD ####
###### Liver Low-DMSO ####

LiLD_genes <- ifelse(Li.res_LD$padj < 0.05, 1, 0)
names(LiLD_genes) <- rownames(Li.res_LD)

LiLD_valid_genes <- intersect(names(LiLD_genes), names(length_vector))

LiLD_genes <- LiLD_genes[LiLD_valid_genes]
length_vector <- length_vector[LiLD_valid_genes]

#Calculate gene weight values
LiLD_pwf <- nullp(LiLD_genes, bias.data = length_vector)
hist(LiLD_pwf$bias.data, 30)
hist(LiLD_pwf$pwf, 30)

LiLD_GO <- goseq(LiLD_pwf, "danRer6", "geneSymbol")
LiLD_GO.padj <- p.adjust(LiLD_GO$over_represented_pvalue, method = "fdr")
sum(LiLD_GO.padj < 0.05)



