#install.packages("BiocManager")
#BiocManager::install("Rsubread")
library(dplyr)
library(DESeq2)
library(edgeR)
library(limma)


zero <- rownames(All_Counts)[rowSums(All_Counts == 0) == ncol(All_Counts)]
length(zero)
#WE have 2484 with zero expression across all samples. 
colSums(All_Counts)


colSums(All_Counts) %>% barplot(., las = 3, ylab = "Reads Mapped Per Sample")

boxplot(as.matrix(All_Counts) ~ col(All_Counts), 
        ylab = "Counts", xlab = "Samples", names = colnames(All_Counts))


#Brain Graphs
colSums(Brain_Counts)
colSums(Brain_Counts) %>% barplot(., las = 3, ylab = "Reads Mapped Per Sample", cex.names = 0.7, main = "Brain Samples")

#Plus 0.5 to avoid log(0) errors
logBrain = log2(as.matrix(Brain_Counts) + 0.5)
boxplot(as.matrix(logBrain) ~ col(Brain_Counts), las=3, cex.axis = 0.4,
        ylab = "Counts", xlab = "Samples", names = colnames(Brain_Counts), main = "Log counts of Brain samples")

#Gonad Graphs
colSums(Gonad_Counts)
colSums(Gonad_Counts) %>% barplot(., las = 3, ylab = "Reads Mapped Per Sample", cex.names = 0.7, main = "Gonad Samples")

logGonad = log2(as.matrix(Gonad_Counts) + 0.5)
boxplot(as.matrix(logGonad) ~ col(Gonad_Counts), las=3, cex.axis = 0.4,
        ylab = "Counts", xlab = "Samples", names = colnames(Gonad_Counts), main = "Log counts of Gonad samples")

#Liver Graphs
colSums(Liver_Counts)
colSums(Liver_Counts) %>% barplot(., las = 3, ylab = "Reads Mapped Per Sample", cex.names = 0.7, main = "Liver Samples")

logLiver = log2(as.matrix(Liver_Counts) + 0.5)
boxplot(as.matrix(logLiver) ~ col(Liver_Counts), las=3, cex.axis = 0.4,
        ylab = "Counts", xlab = "Samples", names = colnames(Liver_Counts), main = "Log counts of Liver samples")


fc_results <- read.table("fc_results.txt", sep= "\t", header=TRUE)
All_Counts <- read.table("All_Counts.txt", sep="\t", header = TRUE,  row.names = 1)
Brain_Counts <- read.table("Brain_Counts.txt", sep="\t", header=TRUE, row.names =1)
Gonad_Counts <- read.table("Gonad_Counts.txt", sep="\t", header=TRUE, row.names = 1)
Liver_Counts <- read.table("Liver_Counts.txt", sep="\t", header=TRUE, row.names = 1)

### DESeq Analysis ####

#### Separate into conditions ####

#HLD Brain
HLD.Brain.cols <- c("D_F0_M1_B", "L_F0_M1_B", "D_F0_M2_B", "D_F0_M3_B", 
                  "L_F0_M2_B", "L_F0_M3_B", "H_F0_M2_B", "H_F0_M1_B", 
                  "H_F0_M3_B", "H_F0_M4_B", "D_F0_M4_B")
H.L.D.Brain <- Brain_Counts[, HLD.Brain.cols]
rownames(H.L.D.Brain) <- rownames(Brain_Counts)

#Mciroplastics Brain
MP.MPD.D.Brain.cols <- c("MP_F0_M1_B", "D_F0_M1_B", "D_F0_M2_B", "MPD_F0_M1_B", "D_F0_M3_B", "MP_F0_M3_B", "MP_F0_M2_B", "MPD_F0_M2_B", "MPD_F0_M3_B", "MP_F0_M4_B", "D_F0_M4_B")
MP.MPD.D.Brain <- Brain_Counts[, MP.MPD.D.Brain.cols]
rownames(MP.MPD.D.Brain) <- rownames(Brain_Counts)

#Brain F1- High and DMSO
F1.Brain.cols <- c("H_F1_M3_B", "H_F1_M3_B.2.", "D_F1_M3_B", "D_F1_M1_B", "D_F1_M1_B.2.", "H_F1_M2_B")
F1.Brain <- Brain_Counts[,F1.Brain.cols]
rownames(F1.Brain) <-rownames(Brain_Counts)

#HLD Gonad
H.L.D.Gonad.cols <- c("D_F0_M2_G", "L_F0_M2_G", "L_F0_M3_G", "H_F0_M2_G", "H_F0_M1_G", "H_F0_M3_G", "D_F0_M3_G", "D_F0_M1_G", "L_F0_M1_G", "H_F0_M4_G", "D_F0_M4_G")
H.L.D.Gonad <-Gonad_Counts[,H.L.D.Gonad.cols]
rownames(H.L.D.Gonad) <- rownames(Gonad_Counts)

#Microplastics Gonad
MP.MPD.D.Gonad.cols <- c("D_F0_M2_G", "MP_F0_M1_G", "MPD_F0_M1_G", "MP_F0_M3_G", "MP_F0_M2_G", "MPD_F0_M2_G", "MPD_F0_M3_G", "D_F0_M3_G", "D_F0_M1_G", "MP_F0_M4_G", "D_F0_M4_G")
MP.MPD.D.Gonad <- Gonad_Counts[, MP.MPD.D.Gonad.cols]
rownames(MP.MPD.D.Gonad) <-rownames(Gonad_Counts)

#F1 Gonads 
F1.Gonad.cols <- c("D_F1_M1_G", "D_F1_M1_G.2.", "D_F1_M3_G", "H_F1_M3_G", "H_F1_M3_G.2.", "H_F1_M2_G")
F1.Gonad <- Gonad_Counts[, F1.Gonad.cols]
rownames(F1.Gonad) <- rownames(Gonad_Counts)

#Liver
#HLD Liver
H.L.D.Liver.cols <- c("L_F0_M1_L", "L_F0_M2_L", "L_F0_M3_L", "H_F0_M1_L", "D_F0_M3_L", "D_F0_M2_L", "D_F0_M1_L", "H_F0_M2_L", "H_F0_M3_L", "H_F0_M4_L", "D_F0_M4_L")
H.L.D.Liver <- Liver_Counts[, H.L.D.Liver.cols]
rownames(H.L.D.Liver) <- rownames(Liver_Counts)

#Microplastics Liver
MP.MPD.D.Liver.cols <- c("MP_F0_M2_L", "MP_F0_M1_L", "MPD_F0_M1_L", "D_F0_M3_L", "D_F0_M2_L", "D_F0_M1_L", "MP_F0_M3_L", "MPD_F0_M2_L", "MPD_F0_M3_L", "MP_F0_M4_L", "D_F0_M4_L")
MP.MPD.D.Liver <- Liver_Counts[, MP.MPD.D.Liver.cols]
rownames(MP.MPD.D.Liver) <- rownames(Liver_Counts)

#F1 Liver
F1.Liver.cols <- c("D_F1_M3_L", "H_F1_M3_L", "H_F1_M3_L.2.", "H_F1_M2_L", "D_F1_M1_L", "D_F1_M1_L.2.")
F1.Liver <- Liver_Counts[, F1.Liver.cols]
rownames(F1.Liver) <- rownames(Liver_Counts)

#### Brain High-Low- DMSO ####

#install.packages("BiocManager")  # if not already installed
#BiocManager::install(version = "3.21")  # or latest version for your R
#BiocManager::install()  # update all Bioconductor packages
#BiocManager::install("DESeq2")
library(DESeq2)

#HLD Brain 
HLD.Br_conds <- c("DMSO", "Low", "DMSO", "DMSO", "Low", "Low", "High", "High", "High", "High", "DMSO")

HLD_metadata <- data.frame(row.names = colnames(H.L.D.Brain),
                           conds = factor(HLD.Br_conds))
HLD.Br_conds_dds <- DESeqDataSetFromMatrix(countData = as.matrix(H.L.D.Brain),
                                           colData = HLD_metadata,
                                           design = ~ conds)

keep <- rowSums(cpm(counts(HLD.Br_conds_dds)) >= 1) >= 6
HLD.Br_conds_dds <- HLD.Br_conds_dds[keep,]

HLD.Br_conds_dds <- DESeq(HLD.Br_conds_dds)

#Pairwise
# High vs Low
Br.res_HL <- results(HLD.Br_conds_dds, contrast = c("conds", "High", "Low"))
# High vs DMSO
Br.res_HD <- results(HLD.Br_conds_dds, contrast = c("conds", "High", "DMSO"))
# Low vs DMSO
Br.res_LD <- results(HLD.Br_conds_dds, contrast = c("conds", "Low", "DMSO"))

#remove na values
Br.res_HL <- na.omit(Br.res_HL)
Br.res_HD <- na.omit(Br.res_HD)
Br.res_LD <- na.omit(Br.res_LD)
#significant p-values
Br.res_s_HL <- Br.res_HL[Br.res_HL$padj <= 0.05 ,]
dim(Br.res_s_HL) #0 significant
Br.res_s_HD <- Br.res_HD[Br.res_HD$padj <= 0.05 ,]
dim(Br.res_s_HD) #1 significant
Br.res_s_LD <- Br.res_LD[Br.res_LD$padj <= 0.05 ,]
dim(Br.res_s_LD) # 1 significant

#PCA
plotPCA(rlog(HLD.Br_conds_dds, blind=TRUE), intgroup="conds")


#### Brain Microplastics ####

#conditions
MP.Br_conds <- c("MP", "DMSO", "DMSO", "MPDMSO", "DMSO", "MP", "MP", "MPDMSO", "MPDMSO", "MP", "DMSO")

MP.Br_metadata <- data.frame(row.names = colnames(MP.MPD.D.Brain),
                           conds = factor(MP.Br_conds))
MP.Br_conds_dds <- DESeqDataSetFromMatrix(countData = as.matrix(MP.MPD.D.Brain),
                                           colData = MP.Br_metadata,
                                           design = ~ conds)

keep <- rowSums(cpm(counts(MP.Br_conds_dds)) >= 1) >= 6
MP.Br_conds_dds <- MP.Br_conds_dds[keep,]

MP.Br_conds_dds <- DESeq(MP.Br_conds_dds)

#Pairwise
#MP vs MPDMSO
Br.res_MP.MPD <- results(MP.Br_conds_dds, contrast = c("conds", "MP", "MPDMSO"))
# MP vs DMSO
Br.res_MP.DMSO <- results(MP.Br_conds_dds, contrast = c("conds", "MP", "DMSO"))
# MPDMSO vs DMSO
Br.res_MPD.DMSO <- results(MP.Br_conds_dds, contrast = c("conds", "MPDMSO", "DMSO"))

#remove na values
Br.res_MP.MPD <- na.omit(Br.res_MP.MPD)
Br.res_MP.DMSO <- na.omit(Br.res_MP.DMSO)
Br.res_MPD.DMSO <- na.omit(Br.res_MPD.DMSO)
#Significant p-values
Br.res_s_MP.MPD <- Br.res_MP.MPD[Br.res_MP.MPD$padj <= 0.05 ,]
dim(Br.res_s_MP.MPD) #0 significant genes
Br.res_s_MP.DMSO <-Br.res_MP.DMSO[Br.res_MP.DMSO$padj <= 0.05 ,]
dim(Br.res_s_MP.DMSO) #2 significant genes
Br.res_s_MPD.DMSO <- Br.res_MPD.DMSO[Br.res_MPD.DMSO$padj <= 0.05 ,]
dim(Br.res_s_MPD.DMSO) # 3 significant gene

#PCA
plotPCA(rlog(MP.Br_conds_dds, blind=TRUE), intgroup="conds")

#### Brain F1 ####

#conditions
F1.Br_conds <- c("High", "High", "DMSO", "DMSO", "DMSO", "High")

F1.Br_metadata <- data.frame(row.names = colnames(F1.Brain),
                             conds = factor(F1.Br_conds))
F1.Br_conds_dds <- DESeqDataSetFromMatrix(countData = as.matrix(F1.Brain),
                                          colData = F1.Br_metadata,
                                          design = ~ conds)

keep <- rowSums(cpm(counts(F1.Br_conds_dds)) >= 1) >= 6
F1.Br_conds_dds <- F1.Br_conds_dds[keep,]

F1.Br_conds_dds <- DESeq(F1.Br_conds_dds)

#results
Br.res_F1 <- results(F1.Br_conds_dds, contrast = c("conds", "High", "DMSO"))
#remove na values
Br.res_F1 <- na.omit(Br.res_F1)
#Significant p-values
Br.res_s_F1 <- Br.res_F1[Br.res_F1$padj <= 0.05 ,]
dim(Br.res_s_F1) #106 significant genes

#PCA
plotPCA(rlog(F1.Br_conds_dds, blind=TRUE), intgroup="conds")

#### Gonad High-Low-DMSO ####

HLD.Go_conds <- c("DMSO", "Low", "Low", "High", "High", "High", "DMSO", "DMSO", "Low", "High", "DMSO")

HLD.Go_metadata <- data.frame(row.names = colnames(H.L.D.Gonad),
                             conds = factor(HLD.Go_conds))
HLD.Go_conds_dds <- DESeqDataSetFromMatrix(countData = as.matrix(H.L.D.Gonad),
                                          colData = HLD.Go_metadata,
                                          design = ~ conds)

keep <- rowSums(cpm(counts(HLD.Go_conds_dds)) >= 1) >= 6
HLD.Go_conds_dds <- HLD.Go_conds_dds[keep,]

HLD.Go_conds_dds <- DESeq(HLD.Go_conds_dds)

#Pairwise
#High vs Low
Go.res_HL <- results(HLD.Go_conds_dds, contrast = c("conds", "High", "Low"))
# High vs DMSO
Go.res_HD <- results(HLD.Go_conds_dds, contrast = c("conds", "High", "DMSO"))
# Low vs DMSO
Go.res_LD <- results(HLD.Go_conds_dds, contrast = c("conds", "Low", "DMSO"))

#remove na values
Go.res_HL <- na.omit(Go.res_HL)
Go.res_HD <- na.omit(Go.res_HD)
Go.res_LD <- na.omit(Go.res_LD)
#Significant p-values
Go.res_s_HL <- Go.res_HL[Go.res_HL$padj <= 0.05 ,]
dim(Go.res_s_HL) #5032 significant genes
Go.res_s_HD <-Go.res_HD[Go.res_HD$padj <= 0.05 ,]
dim(Go.res_s_HD) #5184 significant genes
Go.res_s_LD <- Go.res_LD[Go.res_LD$padj <= 0.05 ,]
dim(Go.res_s_LD) # 3 significant gene

#PCA
plotPCA(rlog(HLD.Go_conds_dds, blind=TRUE), intgroup="conds")

#### Gonad Microplastics ####

#conditions
MP.Go_conds <- c("DMSO", "MP", "MPDMSO", "MP", "MP", "MPDMSO", "MPDMSO", "DMSO", "DMSO", "MP", "DMSO")

MP.Go_metadata <- data.frame(row.names = colnames(MP.MPD.D.Gonad),
                             conds = factor(MP.Go_conds))
MP.Go_conds_dds <- DESeqDataSetFromMatrix(countData = as.matrix(MP.MPD.D.Gonad),
                                          colData = MP.Go_metadata,
                                          design = ~ conds)

keep <- rowSums(cpm(counts(MP.Go_conds_dds)) >= 1) >= 6
MP.Go_conds_dds <- MP.Go_conds_dds[keep,]

MP.Go_conds_dds <- DESeq(MP.Go_conds_dds)

#Pairwise
#MP vs MPDMSO
Go.res_MP.MPD <- results(MP.Go_conds_dds, contrast = c("conds", "MP", "MPDMSO"))
# MP vs DMSO
Go.res_MP.DMSO <- results(MP.Go_conds_dds, contrast = c("conds", "MP", "DMSO"))
# MPDMSO vs DMSO
Go.res_MPD.DMSO <- results(MP.Go_conds_dds, contrast = c("conds", "MPDMSO", "DMSO"))

#remove na values
Go.res_MP.MPD <- na.omit(Go.res_MP.MPD)
Go.res_MP.DMSO <- na.omit(Go.res_MP.DMSO)
Go.res_MPD.DMSO <- na.omit(Go.res_MPD.DMSO)
#Significant p-values
Go.res_s_MP.MPD <- Go.res_MP.MPD[Go.res_MP.MPD$padj <= 0.05 ,]
dim(Go.res_s_MP.MPD) #1517 significant genes
Go.res_s_MP.DMSO <-Go.res_MP.DMSO[Go.res_MP.DMSO$padj <= 0.05 ,]
dim(Go.res_s_MP.DMSO) #5758 significant genes
Go.res_s_MPD.DMSO <- Go.res_MPD.DMSO[Go.res_MPD.DMSO$padj <= 0.05 ,]
dim(Go.res_s_MPD.DMSO) # 9 significant gene

#PCA
plotPCA(rlog(MP.Go_conds_dds, blind=TRUE), intgroup="conds")

#### Gonad F1 ####

#conditions
F1.Go_conds <- c("DMSO", "DMSO", "DMSO", "High", "High", "High")

F1.Go_metadata <- data.frame(row.names = colnames(F1.Gonad),
                             conds = factor(F1.Go_conds))
F1.Go_conds_dds <- DESeqDataSetFromMatrix(countData = as.matrix(F1.Gonad),
                                          colData = F1.Go_metadata,
                                          design = ~ conds)

keep <- rowSums(cpm(counts(F1.Go_conds_dds)) >= 1) >= 6
F1.Go_conds_dds <- F1.Go_conds_dds[keep,]

F1.Go_conds_dds <- DESeq(F1.Go_conds_dds)



#results
Go.res_F1 <- results(F1.Go_conds_dds, contrast = c("conds", "High", "DMSO"))
#remove na values
Go.res_F1 <- na.omit(Go.res_F1)
#Significant p-values
Go.res_s_F1 <- Go.res_F1[Go.res_F1$padj <= 0.05 ,]
dim(Go.res_s_F1) #71 significant genes

#PCA
plotPCA(rlog(F1.Go_conds_dds, blind=TRUE), intgroup="conds")

#### Liver High-Low-DMSO ####

HLD.Li_conds <- c("Low", "Low", "Low", "High", "DMSO", "DMSO", "DMSO", "High", "High", "High", "DMSO")

HLD.Li_metadata <- data.frame(row.names = colnames(H.L.D.Liver),
                              conds = factor(HLD.Li_conds))
HLD.Li_conds_dds <- DESeqDataSetFromMatrix(countData = as.matrix(H.L.D.Liver),
                                           colData = HLD.Li_metadata,
                                           design = ~ conds)

keep <- rowSums(cpm(counts(HLD.Li_conds_dds)) >= 1) >= 6
HLD.Li_conds_dds <- HLD.Li_conds_dds[keep,]

HLD.Li_conds_dds <- DESeq(HLD.Li_conds_dds)

#Pairwise
#High vs Low
Li.res_HL <- results(HLD.Li_conds_dds, contrast = c("conds", "High", "Low"))
# High vs DMSO
Li.res_HD <- results(HLD.Li_conds_dds, contrast = c("conds", "High", "DMSO"))
# Low vs DMSO
Li.res_LD <- results(HLD.Li_conds_dds, contrast = c("conds", "Low", "DMSO"))

#remove na values
Li.res_HL <- na.omit(Li.res_HL)
Li.res_HD <- na.omit(Li.res_HD)
Li.res_LD <- na.omit(Li.res_LD)
#Significant p-values
Li.res_s_HL <- Li.res_HL[Li.res_HL$padj <= 0.05 ,]
dim(Li.res_s_HL) #4 significant genes
Li.res_s_HD <-Li.res_HD[Li.res_HD$padj <= 0.05 ,]
dim(Li.res_s_HD) #2 significant genes
Li.res_s_LD <- Li.res_LD[Li.res_LD$padj <= 0.05 ,]
dim(Li.res_s_LD) # 415 significant gene

#PCA
plotPCA(rlog(HLD.Li_conds_dds, blind=TRUE), intgroup="conds")

#### Liver Microplastics ####

MP.Li_conds <- c("MP", "MP", "MPDMSO", "DMSO", "DMSO", "DMSO", "MP", "MPDMSO", "MPDMSO", "MP", "DMSO")

MP.Li_metadata <- data.frame(row.names = colnames(MP.MPD.D.Liver),
                             conds = factor(MP.Li_conds))
MP.Li_conds_dds <- DESeqDataSetFromMatrix(countData = as.matrix(MP.MPD.D.Liver),
                                          colData = MP.Li_metadata,
                                          design = ~ conds)

keep <- rowSums(cpm(counts(MP.Li_conds_dds)) >= 1) >= 6
MP.Li_conds_dds <- MP.Li_conds_dds[keep,]

MP.Li_conds_dds <- DESeq(MP.Li_conds_dds)

#Pairwise
#MP vs MPDMSO
Li.res_MP.MPD <- results(MP.Li_conds_dds, contrast = c("conds", "MP", "MPDMSO"))
# MP vs DMSO
Li.res_MP.DMSO <- results(MP.Li_conds_dds, contrast = c("conds", "MP", "DMSO"))
# MPDMSO vs DMSO
Li.res_MPD.DMSO <- results(MP.Li_conds_dds, contrast = c("conds", "MPDMSO", "DMSO"))

#remove na values
Li.res_MP.MPD <- na.omit(Li.res_MP.MPD)
Li.res_MP.DMSO <- na.omit(Li.res_MP.DMSO)
Li.res_MPD.DMSO <- na.omit(Li.res_MPD.DMSO)
#Significant p-values
Li.res_s_MP.MPD <- Li.res_MP.MPD[Li.res_MP.MPD$padj <= 0.05 ,]
dim(Li.res_s_MP.MPD) #0 significant genes
Li.res_s_MP.DMSO <-Li.res_MP.DMSO[Li.res_MP.DMSO$padj <= 0.05 ,]
dim(Li.res_s_MP.DMSO) #0 significant genes
Li.res_s_MPD.DMSO <- Li.res_MPD.DMSO[Li.res_MPD.DMSO$padj <= 0.05 ,]
dim(Li.res_s_MPD.DMSO) # 0 significant gene

#PCA
plotPCA(rlog(MP.Li_conds_dds, blind=TRUE), intgroup="conds")

#### F1 Liver ####

#conditions

F1.Li_conds <- c("DMSO", "High", "High", "High", "DMSO", "DMSO")

F1.Li_metadata <- data.frame(row.names = colnames(F1.Liver),
                             conds = factor(F1.Li_conds))
F1.Li_conds_dds <- DESeqDataSetFromMatrix(countData = as.matrix(F1.Liver),
                                          colData = F1.Li_metadata,
                                          design = ~ conds)

keep <- rowSums(cpm(counts(F1.Li_conds_dds)) >= 1) >= 6
F1.Li_conds_dds <- F1.Li_conds_dds[keep,]

F1.Li_conds_dds <- DESeq(F1.Li_conds_dds)

#results
Li.res_F1 <- results(F1.Li_conds_dds, contrast = c("conds", "High", "DMSO"))
#remove na values
Li.res_F1 <- na.omit(Li.res_F1)
#Significant p-values
Li.res_s_F1 <- Li.res_F1[Li.res_F1$padj <= 0.05 ,]
dim(Li.res_s_F1) #52 significant genes

#PCA
plotPCA(rlog(F1.Li_conds_dds, blind=TRUE), intgroup="conds")

#Remove metadata
rm(HLD.Li_metadata, HLD_metadata, HLD.Go_metadata, MP.Br_metadata, MP.Go_metadata, MP.Li_metadata, F1.Br_metadata, F1.Go_metadata, F1.Li_metadata)

#What are our genes of interest
head(Br.res_s_HD)
head(Br.res_s_LD)

head(Br.res_s_MP.DMSO)
head(Br.res_s_MPD.DMSO)

head(Br.res_s_F1, 10)


### Only keep results - Free space
#rm(HLD.Br_conds_dds, MP.Br_conds_dds, F1.Br_conds_dds, HLD.Go_conds_dds,  MP.Go_conds_dds, F1.Go_conds_dds, HLD.Li_conds_dds, MP.Li_conds_dds, F1.Li_conds_dds, keep)

### Save the significant genes 
write.csv(as.data.frame(Br.res_s_HD), file = "Br_res_s_HD.csv")
write.csv(as.data.frame(Br.res_s_LD), file = "Br_res_s_LD.csv")
write.csv(as.data.frame(Br.res_s_HL), file = "Br_res_s_HL.csv")
write.csv(as.data.frame(Br.res_s_MP.DMSO), file = "Br.res_s_MP.DMSO.csv")
write.csv(as.data.frame(Br.res_s_MP.MPD), file = "Br.res_s_MP.MPD.csv")
write.csv(as.data.frame(Br.res_s_MPD.DMSO), file = "Br.res_s_MPD.DMSO.csv")
write.csv(as.data.frame(Br.res_s_F1), file = "Br.res_s_F1.csv")

write.csv(as.data.frame(Go.res_s_HD), file = "Go.res_s_HD.csv")
write.csv(as.data.frame(Go.res_s_LD), file = "Go.res_s_LD.csv")
write.csv(as.data.frame(Go.res_s_HL), file = "Go.res_s_HL.csv")
write.csv(as.data.frame(Go.res_s_MP.DMSO), file = "Go.res_s_MP.DMSO.csv")
write.csv(as.data.frame(Go.res_s_MP.MPD), file = "Go.res_s_MP.MPD.csv")
write.csv(as.data.frame(Go.res_s_MPD.DMSO), file = "Go.res_s_MPD.DMSO.csv")
write.csv(as.data.frame(Go.res_s_F1), file = "Go.res_s_F1.csv")

write.csv(as.data.frame(Li.res_s_HD), file = "Li.res_s_HD.csv")
write.csv(as.data.frame(Li.res_s_LD), file = "Li.res_s_LD.csv")
write.csv(as.data.frame(Li.res_s_HL), file = "Li.res_s_HL.csv")
write.csv(as.data.frame(Li.res_s_MP.DMSO), file = "Li.res_s_MP.DMSO.csv")
write.csv(as.data.frame(Li.res_s_MP.MPD), file = "Li.res_s_MP.MPD.csv")
write.csv(as.data.frame(Li.res_s_MPD.DMSO), file = "Li.res_s_MPD.DMSO.csv")
write.csv(as.data.frame(Li.res_s_F1), file = "Li.res_s_F1.csv")

