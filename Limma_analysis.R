
#Requires the conditions from the DESeq analysis



### Limma Analysis ####
library(limma)
#BiocManager::install("edgeR")
library(edgeR)
#### Brain High-Low-DMSO ####
#Design Matrix - conditions = HLD.Br_conds
#~0+treatments makes it a means model rather than a means-reference model = can be used for pairwise comparison
HLDBr.design = model.matrix(~0+HLD.Br_conds)
Lm_HLD.Br = DGEList(counts=H.L.D.Brain)
Lm_HLD.Br = calcNormFactors(Lm_HLD.Br)
#In case of plotting - might not use
#Lm_HLD.BrlogCPM = cpm(Lm_HLD.Br, log=TRUE, prior.count=3)

keep <- filterByExpr(Lm_HLD.Br, HLDBr.design)
Lm_HLD.Br <- Lm_HLD.Br[keep, , keep.lib.sizes=FALSE]

#voom transformation
v.HLDBr <- voom(Lm_HLD.Br, HLDBr.design, plot = TRUE) 
fit.HLDBr <- lmFit(v.HLDBr, HLDBr.design)

#Defining Pairwise
HLD.Brcontrast <- makeContrasts(
  Br.DMSO_vs_High = HLD.Br_condsDMSO - HLD.Br_condsHigh,
  Br.DMSO_vs_Low = HLD.Br_condsDMSO - HLD.Br_condsLow,
  Br.High_vs_Low = HLD.Br_condsHigh - HLD.Br_condsLow,
  levels = HLDBr.design
)

fit2HLD.Br <- contrasts.fit(fit.HLDBr, HLD.Brcontrast)
fit2HLD.Br <- eBayes(fit2HLD.Br)

#Extracting results
#DMSO vs High
results_Br.DMSO_vs_High <- topTable(fit2HLD.Br, coef = "Br.DMSO_vs_High", number = Inf, adjust.method = "BH")

#DMSO vs Low
results_Br.DMSO_vs_Low <- topTable(fit2HLD.Br, coef = "Br.DMSO_vs_Low", number = Inf, adjust.method = "BH")

#High vs Low
results_Br.High_vs_Low <- topTable(fit2HLD.Br, coef = "Br.High_vs_Low", number = Inf, adjust.method = "BH")

#Significant p-values
results_Br.DMSO_vs_High <- results_Br.DMSO_vs_High[results_Br.DMSO_vs_High$adj.P.Val <= 0.05, ]
dim(results_Br.DMSO_vs_High) #0 significant genes

results_Br.DMSO_vs_Low <- results_Br.DMSO_vs_Low[results_Br.DMSO_vs_Low$adj.P.Val <= 0.05, ]
dim(results_Br.DMSO_vs_Low) # 0 significant genes

results_Br.High_vs_Low <- results_Br.High_vs_Low[results_Br.High_vs_Low$adj.P.Val <= 0.05, ]
dim(results_Br.DMSO_vs_Low) #0 significant genes

#### Brain Microplastics ####
MP.Br.design = model.matrix(~0+MP.Br_conds)
Lm_MP.Br = DGEList(counts=MP.MPD.D.Brain)
Lm_MP.Br = calcNormFactors(Lm_MP.Br)
keep <- filterByExpr(Lm_MP.Br, MP.Br.design)
Lm_MP.Br <- Lm_MP.Br[keep, , keep.lib.sizes=FALSE]

#voom transformation
v.MPBr <- voom(Lm_MP.Br, MP.Br.design, plot = TRUE) 
fit.MPBr <- lmFit(v.MPBr, MP.Br.design)

#Defining Pairwise
MP.Brcontrast <- makeContrasts(
  Br.DMSO_vs_MP = MP.Br_condsDMSO - MP.Br_condsMP,
  Br.DMSO_vs_MPDMSO = MP.Br_condsDMSO - MP.Br_condsMPDMSO,
  Br.MP_vs_MPDMSO = MP.Br_condsMP - MP.Br_condsMPDMSO,
  levels = MP.Br.design
)

fit2MP.Br <- contrasts.fit(fit.MPBr, MP.Brcontrast)
fit2MP.Br <- eBayes(fit2MP.Br)

#Extracting results
#DMSO vs MP
results_Br.DMSO_vs_MP <- topTable(fit2MP.Br, coef = "Br.DMSO_vs_MP", number = Inf, adjust.method = "BH")

#DMSO vs MPDMSO
results_Br.DMSO_vs_MPDMSO <- topTable(fit2MP.Br, coef = "Br.DMSO_vs_MPDMSO", number = Inf, adjust.method = "BH")

#MP vs MPDMSO
results_Br.MP_vs_MPDMSO <- topTable(fit2MP.Br, coef = "Br.MP_vs_MPDMSO", number = Inf, adjust.method = "BH")

#Significant p-values
results_Br.DMSO_vs_MP <- results_Br.DMSO_vs_MP[results_Br.DMSO_vs_MP$adj.P.Val <= 0.05, ]
dim(results_Br.DMSO_vs_MP) #0 significant genes

results_Br.DMSO_vs_MPDMSO <- results_Br.DMSO_vs_MPDMSO[results_Br.DMSO_vs_MPDMSO$adj.P.Val <= 0.05, ]
dim(results_Br.DMSO_vs_MPDMSO) # 0 significant genes

results_Br.MP_vs_MPDMSO <- results_Br.MP_vs_MPDMSO[results_Br.MP_vs_MPDMSO$adj.P.Val <= 0.05, ]
dim(results_Br.MP_vs_MPDMSO) #0 significant genes




#### Brain F1 ####

F1.Br.design = model.matrix(~F1.Br_conds)
Lm_F1.Br = DGEList(counts=F1.Brain)
Lm_F1.Br = calcNormFactors(Lm_F1.Br)
keep <- filterByExpr(Lm_F1.Br, F1.Br.design)
Lm_F1.Br <- Lm_F1.Br[keep, , keep.lib.sizes=FALSE]

#voom transformation
v.F1Br <- voom(Lm_F1.Br, F1.Br.design, plot = TRUE) 
fit.F1Br <- lmFit(v.F1Br, F1.Br.design)
fit.F1Br <- eBayes(fit.F1Br)

#Extracting results
#DMSO vs High
results_Br.F1 <- topTable(fit.F1Br, coef = ncol(F1.Br.design), number = Inf, adjust.method = "BH")

#Significant p-values
results_Br.F1 <- results_Br.F1[results_Br.F1$adj.P.Val <= 0.05, ]
dim(results_Br.F1) #0 significant genes

#### Gonad High-Low-DMSO ####
HLD.Go.design = model.matrix(~0+HLD.Go_conds)
Lm_HLD.Go = DGEList(counts=H.L.D.Gonad)
Lm_HLD.Go = calcNormFactors(Lm_HLD.Go)
keep <- filterByExpr(Lm_HLD.Go, HLD.Go.design)
Lm_HLD.Go <- Lm_HLD.Go[keep, , keep.lib.sizes=FALSE]

#voom transformation
v.HLDGo <- voom(Lm_HLD.Go, HLD.Go.design, plot = TRUE) 
fit.HLDGo <- lmFit(v.HLDGo, HLD.Go.design)

#Defining Pairwise
HLD.Gocontrast <- makeContrasts(
  Go.DMSO_vs_High = HLD.Go_condsDMSO - HLD.Go_condsHigh,
  Go.DMSO_vs_Low = HLD.Go_condsDMSO - HLD.Go_condsLow,
  Go.High_vs_Low = HLD.Go_condsHigh - HLD.Go_condsLow,
  levels = HLD.Go.design
)

fit2HLD.Go <- contrasts.fit(fit.HLDGo, HLD.Gocontrast)
fit2HLD.Go <- eBayes(fit2HLD.Go)

#Extracting results
#DMSO vs High
results_Go.DMSO_vs_High <- topTable(fit2HLD.Go, coef = "Go.DMSO_vs_High", number = Inf, adjust.method = "BH")

#DMSO vs Low
results_Go.DMSO_vs_Low <- topTable(fit2HLD.Go, coef = "Go.DMSO_vs_Low", number = Inf, adjust.method = "BH")

#High vs Low
results_Go.High_vs_Low <- topTable(fit2HLD.Go, coef = "Go.High_vs_Low", number = Inf, adjust.method = "BH")


#Significant p-values
results_Go.DMSO_vs_High <- results_Go.DMSO_vs_High[results_Go.DMSO_vs_High$adj.P.Val <= 0.05, ]
dim(results_Go.DMSO_vs_High) #0 significant genes

results_Go.DMSO_vs_Low <- results_Go.DMSO_vs_Low[results_Go.DMSO_vs_Low$adj.P.Val <= 0.05, ]
dim(results_Go.DMSO_vs_Low) # 0 significant genes

results_Go.High_vs_Low <- results_Go.High_vs_Low[results_Go.High_vs_Low$adj.P.Val <= 0.05, ]
dim(results_Go.High_vs_Low) #0 significant genes

#### Gonad Microplastics ####
MP.Go.design = model.matrix(~0+MP.Go_conds)
Lm_MP.Go = DGEList(counts=MP.MPD.D.Gonad)
Lm_MP.Go = calcNormFactors(Lm_MP.Go)
keep <- filterByExpr(Lm_MP.Go, MP.Go.design)
Lm_MP.Go <- Lm_MP.Go[keep, , keep.lib.sizes=FALSE]
#voom transformation
v.MPGo <- voom(Lm_MP.Go, MP.Go.design, plot = TRUE) 
fit.MPGo <- lmFit(v.MPGo, MP.Go.design)

#Defining Pairwise
MP.Gocontrast <- makeContrasts(
  Go.DMSO_vs_MP = MP.Go_condsDMSO - MP.Go_condsMP,
  Go.DMSO_vs_MPDMSO = MP.Go_condsDMSO - MP.Go_condsMPDMSO,
  Go.MP_vs_MPDMSO = MP.Go_condsMP - MP.Go_condsMPDMSO,
  levels = MP.Go.design
)

fit2MP.Go <- contrasts.fit(fit.MPGo, MP.Gocontrast)
fit2MP.Go <- eBayes(fit2MP.Go)

#Extracting results
#DMSO vs MP
results_Go.DMSO_vs_MP <- topTable(fit2MP.Go, coef = "Go.DMSO_vs_MP", number = Inf, adjust.method = "BH")

#DMSO vs MPDMSO
results_Go.DMSO_vs_MPDMSO <- topTable(fit2MP.Go, coef = "Go.DMSO_vs_MPDMSO", number = Inf, adjust.method = "BH")

#MP vs MPDMSO
results_Go.MP_vs_MPDMSO <- topTable(fit2MP.Go, coef = "Go.MP_vs_MPDMSO", number = Inf, adjust.method = "BH")

#Significant p-values
results_Go.DMSO_vs_MP <- results_Go.DMSO_vs_MP[results_Go.DMSO_vs_MP$adj.P.Val <= 0.05, ]
dim(results_Go.DMSO_vs_MP) #0 significant genes

results_Go.DMSO_vs_MPDMSO <- results_Go.DMSO_vs_MPDMSO[results_Go.DMSO_vs_MPDMSO$adj.P.Val <= 0.05, ]
dim(results_Go.DMSO_vs_MPDMSO) # 0 significant genes

results_Go.MP_vs_MPDMSO <- results_Go.MP_vs_MPDMSO[results_Go.MP_vs_MPDMSO$adj.P.Val <= 0.05, ]
dim(results_Go.MP_vs_MPDMSO) #0 significant genes

#### Gonad F1 ####
F1.Go.design = model.matrix(~F1.Go_conds)
Lm_F1.Go = DGEList(counts=F1.Gonad)
Lm_F1.Go = calcNormFactors(Lm_F1.Go)
keep <- filterByExpr(Lm_F1.Go, F1.Go.design)
Lm_F1.Go <- Lm_F1.Go[keep, , keep.lib.sizes=FALSE]
#voom transformation
v.F1Go <- voom(Lm_F1.Go, F1.Go.design, plot = TRUE) 
fit.F1Go <- lmFit(v.F1Go, F1.Go.design)
fit.F1Go <- eBayes(fit.F1Go)

#Extracting results
#DMSO vs High
results_Go.F1 <- topTable(fit.F1Go, coef = ncol(F1.Go.design), number = Inf, adjust.method = "BH")

#Significant p-values
results_Go.F1 <- results_Go.F1[results_Go.F1$adj.P.Val <= 0.05, ]
dim(results_Go.F1) #0 significant genes



#### Liver High-Low-DMSO ####
HLD.Li.design = model.matrix(~0+HLD.Li_conds)
Lm_HLD.Li = DGEList(counts=H.L.D.Liver)
Lm_HLD.Li = calcNormFactors(Lm_HLD.Li)
keep <- filterByExpr(Lm_HLD.Li, HLD.Li.design)
Lm_HLD.Li <- Lm_HLD.Li[keep, , keep.lib.sizes=FALSE]
#voom transformation
v.HLDLi <- voom(Lm_HLD.Li, HLD.Li.design, plot = TRUE) 
fit.HLDLi <- lmFit(v.HLDLi, HLD.Li.design)

#Defining Pairwise
HLD.Licontrast <- makeContrasts(
  Li.DMSO_vs_High = HLD.Li_condsDMSO - HLD.Li_condsHigh,
  Li.DMSO_vs_Low = HLD.Li_condsDMSO - HLD.Li_condsLow,
  Li.High_vs_Low = HLD.Li_condsHigh - HLD.Li_condsLow,
  levels = HLD.Li.design
)

fit2HLD.Li <- contrasts.fit(fit.HLDLi, HLD.Licontrast)
fit2HLD.Li <- eBayes(fit2HLD.Li)

#Extracting results
#DMSO vs High
results_Li.DMSO_vs_High <- topTable(fit2HLD.Li, coef = "Li.DMSO_vs_High", number = Inf, adjust.method = "BH")

#DMSO vs Low
results_Li.DMSO_vs_Low <- topTable(fit2HLD.Li, coef = "Li.DMSO_vs_Low", number = Inf, adjust.method = "BH")

#High vs Low
results_Li.High_vs_Low <- topTable(fit2HLD.Li, coef = "Li.High_vs_Low", number = Inf, adjust.method = "BH")


#Significant p-values
results_Li.DMSO_vs_High <- results_Li.DMSO_vs_High[results_Li.DMSO_vs_High$adj.P.Val <= 0.05, ]
dim(results_Li.DMSO_vs_High) #0 significant genes

results_Li.DMSO_vs_Low <- results_Li.DMSO_vs_Low[results_Li.DMSO_vs_Low$adj.P.Val <= 0.05, ]
dim(results_Li.DMSO_vs_Low) # 0 significant genes

results_Li.High_vs_Low <- results_Li.High_vs_Low[results_Li.High_vs_Low$adj.P.Val <= 0.05, ]
dim(results_Li.High_vs_Low) #0 significant genes


#### Liver Microplastics ####
MP.Li.design = model.matrix(~0+MP.Li_conds)
Lm_MP.Li = DGEList(counts=MP.MPD.D.Liver)
Lm_MP.Li = calcNormFactors(Lm_MP.Li)
keep <- filterByExpr(Lm_MP.Li, MP.Li.design)
Lm_MP.Li <- Lm_MP.Li[keep, , keep.lib.sizes=FALSE]
#voom transformation
v.MPLi <- voom(Lm_MP.Li, MP.Li.design, plot = TRUE) 
fit.MPLi <- lmFit(v.MPLi, MP.Li.design)

#Defining Pairwise
MP.Licontrast <- makeContrasts(
  Li.DMSO_vs_MP = MP.Li_condsDMSO - MP.Li_condsMP,
  Li.DMSO_vs_MPDMSO = MP.Li_condsDMSO - MP.Li_condsMPDMSO,
  Li.MP_vs_MPDMSO = MP.Li_condsMP - MP.Li_condsMPDMSO,
  levels = MP.Li.design
)

fit2MP.Li <- contrasts.fit(fit.MPLi, MP.Licontrast)
fit2MP.Li <- eBayes(fit2MP.Li)

#Extracting results
#DMSO vs MP
results_Li.DMSO_vs_MP <- topTable(fit2MP.Li, coef = "Li.DMSO_vs_MP", number = Inf, adjust.method = "BH")

#DMSO vs MPDMSO
results_Li.DMSO_vs_MPDMSO <- topTable(fit2MP.Li, coef = "Li.DMSO_vs_MPDMSO", number = Inf, adjust.method = "BH")

#MP vs MPDMSO
results_Li.MP_vs_MPDMSO <- topTable(fit2MP.Li, coef = "Li.MP_vs_MPDMSO", number = Inf, adjust.method = "BH")

#Significant p-values
results_Li.DMSO_vs_MP <- results_Li.DMSO_vs_MP[results_Li.DMSO_vs_MP$adj.P.Val <= 0.05, ]
dim(results_Li.DMSO_vs_MP) #0 significant genes

results_Li.DMSO_vs_MPDMSO <- results_Li.DMSO_vs_MPDMSO[results_Li.DMSO_vs_MPDMSO$adj.P.Val <= 0.05, ]
dim(results_Li.DMSO_vs_MPDMSO) # 0 significant genes

results_Li.MP_vs_MPDMSO <- results_Li.MP_vs_MPDMSO[results_Li.MP_vs_MPDMSO$adj.P.Val <= 0.05, ]
dim(results_Li.MP_vs_MPDMSO) #0 significant genes


#### Liver F1 ####
F1.Li.design = model.matrix(~F1.Li_conds)
Lm_F1.Li = DGEList(counts=F1.Liver)
Lm_F1.Li = calcNormFactors(Lm_F1.Li)
keep <- filterByExpr(Lm_F1.Li, F1.Li.design)
Lm_F1.Li <- Lm_F1.Li[keep, , keep.lib.sizes=FALSE]
#voom transformation
v.F1Li <- voom(Lm_F1.Li, F1.Li.design, plot = TRUE) 
fit.F1Li <- lmFit(v.F1Li, F1.Li.design)
fit.F1Li <- eBayes(fit.F1Li)

#Extracting results
#DMSO vs High
results_Li.F1 <- topTable(fit.F1Li, coef = ncol(F1.Li.design), number = Inf, adjust.method = "BH")

#Significant p-values
results_Li.F1 <- results_Li.F1[results_Li.F1$adj.P.Val <= 0.05, ]
dim(results_Li.F1) #0 significant genes