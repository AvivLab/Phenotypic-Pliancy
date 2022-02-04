# Find all non-cancer cell types in Head & Neck Cancer single-cell gene expression data
# & run step 2 of RNA-seq preprocessing and analysis (document on mlambros@mail.einstein.yu.edu
# google drive titled "RNA-sequecing analysis notes")
# AND
# analyze single-cell RNA-seq datasets (head & neck cancer and ovarian cancer datasets) using
# Seurat package and our methods for analyzing phenotypic pliancy

# Set working directory for where all files (except head & neck cancer data, unless you put them all together) are kept:
workingDir = "/home/maryl/Desktop/Feb2018R01Grant/ResubmissionDocumentsAndAnalysis/OvarianScRnaSeqDataAndAnalysis/GSE118828_RAW/"
setwd(workingDir)

# set working directory where head and neck cancer data is kept:
headAndNeckWorkingDir = "/home/maryl/Desktop/Feb2018R01Grant/ResubmissionDocumentsAndAnalysis/SingleCellTranscriptomicAnalysisOfPrimaryAndMetastaticTumorEcosystemsInHeadAndNeckCancer_CorrectVersion"


### Can either run preprocessing using Seurat package (here) OR below starting at "run step 2..." which should end in the same results as the Seurat preprocessing:
## Seurat package preprocessing:
# Need gene expression matrix only for Seurat package, so re-format data here:
library("Seurat")
library("sctransform")
library("ggplot2")
library("SCnorm")
#library("slingshot") # for pseudotime trajectory & functional analysis
##################################################
# Functions used below:
convertLabelOvarian <- function(x) {
  if (grepl("Primary*",x)) {
    x<- "Primary"
  } else if (grepl("Metastatic*",x)) {
    x<- "Metastatic"
  }
  else if (grepl("Normal*",x)) {
    x<- "Normal"
  }
  x
}
convertLabel<- function(x) {
  if (grepl("normal*",x)) {
    x<- "Normal"
  } else if (grepl("lymph*",x)) {
    x<- "Lymph"
  } else if (grepl("metast*",x)) {
    x<- "Metastatic"
  } else if (grepl("primary*",x)) {
    x<- "Primary"
  }
  x
}
difExpBetweenTypes <- function(seuratObjData,geneSymbols,testToUse,ident1String,ident2String){
  differences <- FindMarkers(seuratObjData, ident.1 = ident1String, ident.2 = ident2String, test.use = testToUse)
  diffinds <- sapply(geneSymbols[1:length(geneSymbols)],function(x) {any(rownames(differences)==x)})
  difExpGenes <- geneSymbols[diffinds]
  return(difExpGenes)
}
difExpPcgGenes <- function(difExpGenes, geneList, listIdentString, ident1String, ident2String){
  # How many PcG mechanism genes (trithorax, polycomb and genes controlled by PcG),
  # trithorax alone, polycomb alone, or genes controlled by PcG alone (--> 4 cases are possible listIdent) are
  # differentially expressed:
  pcgIndsDifExp = sapply(difExpGenes,function(x) {any(geneList==x)})
  pcgGenesDifExp = difExpGenes[pcgIndsDifExp]
  print(paste("There are", length(pcgGenesDifExp), "out of", length(pcgnames), listIdentString, "genes differentially expressed btw", ident1String ,"&", ident2String))
  print(pcgGenesDifExp)
  return(pcgGenesDifExp)
}

difExpPcgMechAndTargetGenes <- function(difExpGenes, controlledGenes, ident1String, ident2String){
  pcgGenesDifExpNotUnderControl = difExpGenes[! difExpGenes %in% controlledGenes]
  pcgGenesDifExpUnderControl = difExpGenes[difExpGenes %in% controlledGenes]
  cat(paste("pcg mechanism genes differentially expressed btw", ident1String ,"&", ident2String, "=", length(pcgGenesDifExpNotUnderControl)))
  cat(" ", pcgGenesDifExpNotUnderControl)
  cat(paste("\npcg target genes differentially expressed btw", ident1String ,"&", ident2String, "=", length(pcgGenesDifExpUnderControl)))
  cat(" ", pcgGenesDifExpUnderControl)
}

#########################################################################################
## PREPROCESSING AND ANALYSES:
# Set what data doing analyses on:
ovarianCancerBool = FALSE # Set if looking at ovarian or head & neck cancer data
usePreprocessHeadAndNeckBool = FALSE # to use strict preprocessing like in Head & Neck paper or less strict like in Ovarian cancer paper
scaleBool = FALSE # If want to scale and center for each gene (so lowers variance of gene expression across all cells and also gives
                  # more weight to more noisy lowly expressed genes and decreases weight of higher expressed genes --> so will see in
                  # PlotCountDepth that density is lower for higher expressed genes when use scaling)
plotNormTestBool = FALSE # If to plot and save analysis of normalization methods accuracy
normTestBool = TRUE # If run code using our normalization method vs Seurat
allDataHeadAndNeckBool = TRUE # If to run with all patients' data OR all patients' normal cells and only primary and metastatic cells from patients with matching primary & metastatic samples
regressCellCycleGenesBool = TRUE
saveSeuratObjBool = TRUE
useSavedSeuratObjBool = TRUE
if (ovarianCancerBool==TRUE){
  seuratObjFilename = 'ovarianOurNormCellCycRegresOutNotStrictPreproc2_4_22.rds'#'ova#rianOurNormalizationCellCycleRegressedOutNotStrictPreprocessingSeuratObj.rds'
} else {
  seuratObjFilename = 'headAndNeckOurNormalizationCellCycleRegressedOutNotStrictPreprocessingSeuratObj.rds' # 'headNeckOurNormCellCycRegresOutNotStrictPreproc2_4_22.rds'
}
subclusteringBool = FALSE
savePlotsAsPdfBool = FALSE
restrictToDeGenesBool = TRUE


if (useSavedSeuratObjBool == FALSE) {
  if (ovarianCancerBool == FALSE) {
    #### For head and neck cancer data testing ####
    if (usePreprocessHeadAndNeckBool == TRUE) {
      # Already preprocessed like in head & neck cancer paper:
      setwd(headAndNeckWorkingDir)
      if (allDataHeadAndNeckBool == TRUE) {
        # If want to use all cell types and patients:
        input_file = "fileForAnalysis85AllCellTypes.gz"#args[1] # a tab-separated table where the first column holds the gene symbol followed by as many columns as data conditions to be analyzed
        labels_file = "labels85"#args[2] # labels for each column in data table (first column not important, the remaining columns will have the condition that the column holds (these will be the levels of the analysis)
      } else {
        # If want to use only primary and metastatic cells from patients with matching primary and metastatic tumors, plus normal & lymph cells from all patients
        input_file = "fileForAnalysis85AllNormalButOnlyMatchingPrimaryAndMetast.gz"#args[1] # a tab-separated table where the first column holds the gene symbol followed by as many columns as data conditions to be analyzed
        labels_file = "labels85OnlyMatchingPrimaryAndMetast" # labels for each column in data table (first column not important, the remaining columns will have the condition that the column holds (these will be the levels of the analysis)
      }
      data=read.delim(input_file, header=F, stringsAsFactors=F, sep="\t")
      labels=as.character(read.delim(labels_file, header=F, stringsAsFactors=F, sep="\t"))
      numRows<-dim(data)[1]
      myinds = labels != "normalMyo" & labels!="normal-Fibro" & labels!="normalMast" & labels!= "normalB"
      mydata = data[1:numRows,myinds]
      allDataAllCellTypesOvarian <- mydata
    } else {
      # Preprocess like below only:
      setwd(headAndNeckWorkingDir)
      if (allDataHeadAndNeckBool == TRUE) {
        # If want to use all cell types and patients:
        allDataAllCellTypes = read.delim("HNSCC_all_data_all_cell_types_TestWithHeaders.gz", header=F, stringsAsFactors=F, sep="\t") # with either no sample ID's ("HNSCC_all_data_all_cell_types_ML.gz") or sample ID's in first row ("HNSCC_all_data_all_cell_types_TestWithHeaders.gz")
      } else {
        # If want to use only primary and metastatic cells from patients with matching primary and metastatic tumors, plus normal & lymph cells from all patients
        allDataAllCellTypes = read.delim("HNSCC_all_data_only_matching_cell_types_with_all_normals_WithHeaders.gz", header=F, stringsAsFactors=F, sep="\t") # with sampleIds
      }
      allDataAllCellTypes[,1] <- gsub("'", "", allDataAllCellTypes[,1])
      allDataAllCellTypes = allDataAllCellTypes[-1,]
      allDataAllCellTypesOvarian <- allDataAllCellTypes
    }
  } else {
    #### For Ovarian Cancer preprocessing and analysis:
    setwd(workingDir)
    sampleTitle <- c("PN1-P",	"BN1-P",	"HG3-M1",	"NA1-P",	"NM1",	"LG1-P",	"HG2F-P",	"HG2F-M",	"HG3-P",	"HG3-M2",	"HG4-P1",	"HG4-P2",	"HG4-P3",	"HG4-P4",	"LG2-M",	"LG2-P",	"HG1-P",	"HG1-M")
    geoAccessions <- c("GSM3348303",	"GSM3348304",	"GSM3348305",	"GSM3348306",	"GSM3348307",	"GSM3348308",	"GSM3348309",	"GSM3348310",	"GSM3348311",	"GSM3348312",	"GSM3348313",	"GSM3348314",	"GSM3348315",	"GSM3348316",	"GSM3348317",	"GSM3348318",	"GSM3348319",	"GSM3348320")
    sampleSite <- c("Primary",	"Primary",	"Metastatic",	"Primary",	"Normal",	"Primary",	"Primary",	"Metastatic",	"Primary",	"Metastatic",	"Primary",	"Primary",	"Primary",	"Primary",	"Metastatic",	"Primary",	"Primary",	"Metastatic")
    allCellTypes = unique(sampleSite) # cell type --> either primary, metastatic or normal
    allDataAllCellTypesOvarian = read.delim("ovarian_all_data_all_cell_types_with_headers", header=F, stringsAsFactors=F, sep="\t") # with sample ID's in first row
    allDataAllCellTypesOvarian[,1] <- gsub("'", "", allDataAllCellTypesOvarian[,1])
    allDataAllCellTypesOvarian = allDataAllCellTypesOvarian[-1,] # remove row with sample ID's
  }
  # Create numeric matrix with only gene expression values:
  geneExpMatrix <-sapply(allDataAllCellTypesOvarian[2:dim(allDataAllCellTypesOvarian)[1],2:dim(allDataAllCellTypesOvarian)[2]],as.numeric)
  rownames(geneExpMatrix) <- allDataAllCellTypesOvarian$V1[-1]
  colnames(geneExpMatrix) <- allDataAllCellTypesOvarian[1, 2:dim(allDataAllCellTypesOvarian)[2]]
  # Create Seurat Obj using the data above:
  if (usePreprocessHeadAndNeckBool == TRUE) {
    seuratObjData = CreateSeuratObject(counts = geneExpMatrix, min.cells = 0, min.features = 0)
  } else {
    # Initial preprocessing happens when CreateSeuratObject creates Seurat object: min.cells = 3
    # means keep gene is has expression in at least 3 cells, min.features = 200 means keep cell
    # if has at least 200 genes with expression:
    seuratObjData = CreateSeuratObject(counts = geneExpMatrix, min.cells = 3, min.features = 200)
  }
  if (normTestBool == TRUE) {
    geneExpMatrixFromSeurat = as.matrix(seuratObjData[["RNA"]]@data)
  }
  # No mitocondrial gene expression in both head & neck cancer dataset and ovarian dataset
  # (but this is code to test if there is mitocondrial expression):
  #seuratObjData[["percent.mt"]] <- PercentageFeatureSet(seuratObjData, pattern = "^MT-")
  #VlnPlot(seuratObjData, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  #########################################################################################
  # Normalize data, find most variable gene features to use for downstream analysis like PCA, and scale the data:
  all.genes <- rownames(seuratObjData)
  ## Using NormalizeData and ScaleData functions:
  seuratObjDataNormalizationMethod <- seuratObjData
  seuratObjDataNormalizationMethod <- NormalizeData(seuratObjDataNormalizationMethod)
  if (scaleBool == TRUE) {
    # If want to scale and center for each gene (so lowers variance of gene expression across all cells and also gives
    # more weight to more noisy lowly expressed genes and decreases weight of higher expressed genes --> so will see in
    # PlotCountDepth that density is lower for higher expressed genes when use scaling)
    seuratObjDataNormalizationMethod <- ScaleData(seuratObjDataNormalizationMethod, features = all.genes, do.scale = TRUE, do.center = TRUE)
  } else {
    seuratObjDataNormalizationMethod <- ScaleData(seuratObjDataNormalizationMethod, features = all.genes, do.scale = FALSE, do.center = FALSE)
  }
  ## Using SCTransform function to normalize and scale data in one function:
  #seuratObjData <- SCTransform(seuratObjData, verbose = FALSE, variable.features.n = 3000, do.scale = scaleBool, do.center = TRUE)
  # Find most variable genes if need to use less genes for downstream analysis (don't need to do this so we use all the genes):
  #seuratObjData = FindVariableFeatures(seuratObjData, nfeatures = 5000)
  #VariableFeaturePlot(seuratObjData)
  ########################################################################################
  # Test normalization results:
  if (normTestBool == TRUE) {
    ##### SCTransform normalized and Scaled/Centered results are stored here:
    #seuratNormalizedAndCenteredData <- as.matrix(seuratObjData[["SCT"]]@scale.data) # @data hold log-normalized of corrected UMI counts; @scale.data contains the residuals (normalized values)
    ##### Seurat normalization with NormalizeData and ScaleData functions:
    seuratNormalizationMethodData <- as.matrix(seuratObjDataNormalizationMethod[["RNA"]]@scale.data)
    ##### Normalization like in Head & neck paper (log transform and center data only)
    normalizeAndCenterGeneExpMatrix <- scale(log(geneExpMatrixFromSeurat+1),center=TRUE,scale=TRUE) # center and scale by cell NOT gene here
    seuratObjData <- CreateSeuratObject(counts = normalizeAndCenterGeneExpMatrix, min.cells = 0, min.features = 0)
    all.genes <- rownames(seuratObjData)
    if (scaleBool == TRUE) {
      # If want to scale and center for each gene (so lowers variance of gene expression across all cells and also gives
      # more weight to more noisy lowly expressed genes and decreases weight of higher expressed genes --> so will see in
      # PlotCountDepth that density is lower for higher expressed genes when use scaling)
      seuratObjData <- ScaleData(seuratObjData, features = all.genes, do.scale = TRUE, do.center = TRUE)
    } else {
      seuratObjData <- ScaleData(seuratObjData, features = all.genes, do.scale = FALSE, do.center = FALSE)
    }
  } else {
    seuratObjData <- seuratObjDataNormalizationMethod
    seuratNormalizationMethodData <- as.matrix(seuratObjDataNormalizationMethod[["RNA"]]@scale.data)
  }
  #########################################################################################
  # Set cell type identities to use for labeling cells in analyses below:
  y <- seuratObjData@assays$RNA@counts@Dimnames[[2]]
  if (ovarianCancerBool == TRUE) {
    # For Ovarian cancer Data:
    labelsToUse<-as.character(sapply(y,convertLabelOvarian))
  } else {
    # For Head & Neck cancer data:
    labelsToUse<-as.character(sapply(y,convertLabel))
  }
  # Add labels as CellType in Seurat Obj:
  seuratObjData$CellType <- labelsToUse
  Idents(seuratObjData) <- "CellType"
  # set normalization matrices column names if wanting to test different normalization methods:
  if (normTestBool == TRUE) {
    colnames(normalizeAndCenterGeneExpMatrix) <- labelsToUse
    #colnames(seuratNormalizedAndCenteredData) <- labelsToUse
    colnames(seuratNormalizationMethodData) <- labelsToUse
  } else {
    colnames(seuratNormalizationMethodData) <- labelsToUse
  }
  #########################################################################################
  # Plot Normalization results for seurat normalization method and head & neck cancer paper's normalization method:
  if (plotNormTestBool) {
    setwd(workingDir)
    if (ovarianCancerBool == TRUE) {
      titleToAdd = "Ovarian"
      matchingTitle = ""
    } else {
        if (usePreprocessHeadAndNeckBool == TRUE) {
          titleToAdd = "PreprocessedStrictHeadAndNeck"
        } else {
          titleToAdd = "PreprocessedNotStrictHeadAndNeck"
        }
        if (allDataHeadAndNeckBool == TRUE) {
          matchingTitle = "AllCells"
        } else {
          matchingTitle = "OnlyMatchingCells"
        }
    }
    if (scaleBool == TRUE) {
      scalingTitle = "AndScaling"
    } else {
      scalingTitle = ""
    }
    pdf(paste("countDepth_EvaluationFor_OurNormalizationMethodWithCentering",scalingTitle,"_",titleToAdd,matchingTitle,"Test.pdf", sep = ""), height=5, width=7)
    plotresults <- plotCountDepth(normalizeAndCenterGeneExpMatrix, Conditions = labelsToUse) # test if need to do single-cell normalization or not (do values line up with 1 or curves not lining up on x axis? If not lining up, then need to use scnorm)
    # Since gene expression increases proportionally with sequencing depth, we expect to find the
    # estimated count-depth relationships near 1 for all genes. This is typically true for bulk RNAseq datasets. However, it does not hold in most single-cell RNA-seq datasets. In this example
    # data, the relationship is quite variable across genes.
    dev.off()
    #pdf(paste("countDepth_EvaluationFor_SeuratSctTransformMethodWithCentering",scalingTitle,"_",titleToAdd,".pdf", sep = ""), height=5, width=7)
    #plotresults <- plotCountDepth(seuratNormalizedAndCenteredData, Conditions = labelsToUse, FilterExpression = -0.5) # test if need to do single-cell normalization or not (do values line up with 1 or curves not lining up on x axis? If not lining up, then need to use scnorm)
    # Since gene expression increases proportionally with sequencing depth, we expect to find the
    # estimated count-depth relationships near 1 for all genes. This is typically true for bulk RNAseq datasets. However, it does not hold in most single-cell RNA-seq datasets. In this example
    # data, the relationship is quite variable across genes.
    #dev.off()
    pdf(paste("countDepth_EvaluationFor_SeuratNormalizeDataFuncMethodWithCentering",scalingTitle,"_",titleToAdd,"Test.pdf", sep = ""), height=5, width=7)
    plotresults <- plotCountDepth(seuratNormalizationMethodData, Conditions = labelsToUse) # test if need to do single-cell normalization or not (do values line up with 1 or curves not lining up on x axis? If not lining up, then need to use scnorm)
    # Since gene expression increases proportionally with sequencing depth, we expect to find the
    # estimated count-depth relationships near 1 for all genes. This is typically true for bulk RNAseq datasets. However, it does not hold in most single-cell RNA-seq datasets. In this example
    # data, the relationship is quite variable across genes.
    dev.off()
  }
  # Plot number of RNA counts vs number of genes:
  dev.new()
  FeatureScatter(seuratObjData, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  # Based on plot above can subset to get rid of outliers if want (head & neck cancer has
  # nFeature_RNA between 3000 and 6000, nCount_RNA is between 500,000 to 1,000,000) (ovarian
  # cancer has nFeature_RNA between ___ and ___, nCount_RNA is between ___ to ___)
  #seuratObjData <- subset(seuratObjData, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
  #########################################################################################
  # Regress out cell cycle genes and total # of UMIs for Ovarian dataset (UMIs not used in Head&Neck data)
  # A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
  # segregate this list into markers of G2/M phase and markers of S phase
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  #seuratObjData <- RunPCA(seuratObjData, features = c(s.genes, g2m.genes), verbose = FALSE)
  #DimPlot(seuratObjData)
  if (regressCellCycleGenesBool == TRUE) {
    seuratObjData <- CellCycleScoring(seuratObjData, s.features = s.genes, g2m.features = g2m.genes)
    seuratObjData <- ScaleData(seuratObjData, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(seuratObjData), do.center = FALSE, do.scale = FALSE)
  }
  # visualize with PCA what data looks like now with respect to cell cycle genes:
  seuratObjData <- RunPCA(seuratObjData, features = c(s.genes, g2m.genes), verbose = FALSE)
  dev.new()
  DimPlot(seuratObjData)
  # Save Seurat object:
  if (saveSeuratObjBool == TRUE) {
    saveRDS(seuratObjData, file = seuratObjFilename)
  }
}
if (useSavedSeuratObjBool == TRUE) {
  setwd(workingDir)
# to open saved seurat object:
  seuratObjData <- readRDS(file = seuratObjFilename)
  # Set cell type identities to use for labeling cells in analyses below:
  y <- seuratObjData@assays$RNA@counts@Dimnames[[2]]
  if (ovarianCancerBool == TRUE) {
    # For Ovarian cancer Data:
    labelsToUse<-as.character(sapply(y,convertLabelOvarian))
  } else {
    # For Head & Neck cancer data:
    labelsToUse<-as.character(sapply(y,convertLabel))
  }
  # Add labels as CellType in Seurat Obj:
  seuratObjData$CellType <- labelsToUse
  Idents(seuratObjData) <- "CellType"
}
cat(paste(seuratObjData@assays$RNA@counts@Dim, c("genes", "cells"), "after preprocessing"))
# print off number of cells that are cancer and non-cancer for both primary and metastatic site:
if (ovarianCancerBool == TRUE) {
  print(paste("Ovarian: # cancer cells in primary site =", length(labelsToUse[labelsToUse=="Primary"])))
  print(paste("Ovarian: # non-cancer cells in primary site =", length(labelsToUse[labelsToUse=="Normal"]))) 
  print(paste("Ovarian: # cancer cells in metastatic site =", length(labelsToUse[labelsToUse=="Metastatic"])))
} else {
  print(paste("Head&Neck: # cancer cells in primary site =", length(labelsToUse[labelsToUse=="Primary"])))
  print(paste("Head&Neck: # non-cancer cells in primary site =", length(labelsToUse[labelsToUse=="Normal"]))) 
  print(paste("Head&Neck: # cancer cells in metastatic site =", length(labelsToUse[labelsToUse=="Metastatic"])))
  print(paste("Head&Neck: # non-cancer cells in metastatic site =", length(labelsToUse[labelsToUse=="Lymph"])))
}
#########################################################################################
# Test over all modality of data set after normalization:
library("modes")
geneExpMatrixToUseForModes = as.matrix(seuratObjData[["RNA"]]@scale.data)
colnames(geneExpMatrixToUseForModes) <- labelsToUse
numGenes = dim(geneExpMatrixToUseForModes)[1]
modesOfGenesPrimary = rep(0, numGenes)
modesOfGenesMetast = rep(0, numGenes)
modesOfGenesNormal = rep(0, numGenes)
for (i in 1:numGenes) {
  modesOfGenesPrimary[i] <- modes(geneExpMatrixToUseForModes[i,labelsToUse=="Primary"])[1,1] # find mode ("Value" output of modes func) of each gene across all different cell types
  modesOfGenesMetast[i] <- modes(geneExpMatrixToUseForModes[i,labelsToUse=="Metastatic"])[1,1]
  modesOfGenesNormal[i] <- modes(geneExpMatrixToUseForModes[i,labelsToUse=="Normal" | labelsToUse=="Lymph"])[1,1]
}
cat(paste("Mean modality of all genes is", mean(c(mean(modesOfGenesNormal),mean(modesOfGenesPrimary),mean(modesOfGenesMetast))), "with std =", sd(c(mean(modesOfGenesNormal),mean(modesOfGenesPrimary),mean(modesOfGenesMetast)))))
hist(c(modesOfGenesPrimary,modesOfGenesMetast,modesOfGenesNormal))
# Test some by visualizing with RidgePlot
uniqueModesOfGenes = unique(c(modesOfGenesPrimary,modesOfGenesMetast,modesOfGenesNormal))
cat(uniqueModesOfGenes)
# Frequency of modes greater than 0 (not unimodal):
cat(paste("# of genes with multi-modes for primary:",length(modesOfGenesPrimary[modesOfGenesPrimary > 0]),"; avg multimodal primary:",mean(modesOfGenesPrimary[modesOfGenesPrimary > 0])))
cat(paste("# of genes with multi-modes for metastatic:",length(modesOfGenesMetast[modesOfGenesMetast > 0]),"; avg multimodal metastatic:",mean(modesOfGenesMetast[modesOfGenesMetast > 0])))
cat(paste("# of genes with multi-modes for normal:",length(modesOfGenesNormal[modesOfGenesNormal > 0]),"; avg multimodal normal:",mean(modesOfGenesNormal[modesOfGenesNormal > 0])))
## Investigate visually the modality if want to:
# c(1:numGenes)[modesOfGenesNormal == 10] # --> use this to help pick which gene position to use in ridgeplot to visualize what having these different modes Values looks like
# dev.new()
# RidgePlot(seuratObjData, features = all.genes[c(2986, 3603)])#all.genes[c(1:numGenes)[modesOfGenesNormal == 10]])
# RidgePlot(seuratObjData, features = c("BTG1","RPL10","GADD45B","HSPA1A","CD81","CALM1"))#rownames(geneExpMatrixToUseForModes)[c(1:numGenes)[modesOfGenes == 7]])
## Test what genes usually mutated in cancer look like:
# dev.new()
# RidgePlot(seuratObjData, features = c("TP53", "PIK3CA", "NOTCH1", "TP63", "CDKN2A")) # --> these genes are not present in Ovarian data
#########################################################################################
# Plot Expression changes of PcG genes:
setwd(workingDir)
pcg_file="completeListPcgGenes.csv"
pcglst<-read.csv(pcg_file, header=FALSE)
geneSymbols = seuratObjData@assays$RNA@counts@Dimnames[[1]]
pcginds<-sapply(geneSymbols[1:length(geneSymbols)],function(x) {any(pcglst==x)})
pcgnames<-geneSymbols[pcginds]
# Pull out complete list of genes controlled by Polycomb or Trithorax complexes so that
# can investigate change in DE and also remove from pliancy z-score analysis
genesControlledByPcg = as.matrix(read.csv("controlledByPcgGenesList.csv", header = FALSE))
controlledByPcg <- sapply(geneSymbols,function(x) {any(genesControlledByPcg==x)})
genesControlledByPcgPresent = geneSymbols[controlledByPcg] # genes controlled by PcG mechanism actually present after preprocessing
# Pull out complete list of Polycomb genes for DE analysis:
polycombGenes = as.matrix(read.csv("polycombGenesList.csv", header = FALSE))
polycombGenesPos <- sapply(geneSymbols,function(x) {any(polycombGenes==x)})
polycombGenesPresent = geneSymbols[polycombGenesPos] # polycomb actually present after preprocessing
# Pull out complete list of trithorax genes for DE analysis:
trithoraxGenes = as.matrix(read.csv("trithoraxGenesList.csv", header = FALSE))
trithoraxGenesPos <- sapply(geneSymbols,function(x) {any(trithoraxGenes==x)})
trithoraxGenesPresent = geneSymbols[trithoraxGenesPos] # trithorax genes actually present after preprocessing
# Print results:
cat(paste("There are", length(pcgnames), "out of", length(pcglst),"PcG genes after preprocessing"))
cat(paste("There are", length(genesControlledByPcgPresent), "out of", length(genesControlledByPcg),"genes controlled by PcG mechanism after preprocessing"))
cat(paste("There are", length(polycombGenesPresent), "out of", length(polycombGenes),"polycomb genes after preprocessing"))
cat(paste("There are", length(trithoraxGenesPresent), "out of", length(trithoraxGenes),"trithorax genes after preprocessing"))
#dev.new()
#VlnPlot(seuratObjData, features = pcgnames[1:9])
#########################################################################################
# Differential Expression Analysis btw Primary & Metastatic
library("MAST")
library("DESeq2")
testUsing = "wilcox"
genesDifExpPrimaryAndMetast = difExpBetweenTypes(seuratObjData,geneSymbols,testUsing,"Primary","Metastatic")
pcgGenesDifExpPrimaryAndMetast = difExpPcgGenes(genesDifExpPrimaryAndMetast, pcglst, "PcG mechanism", "Primary","Metastatic")
cat(paste(length(genesDifExpPrimaryAndMetast),"genes differentially expressed between primary & metastatic"))
difExpPcgMechAndTargetGenes(pcgGenesDifExpPrimaryAndMetast, genesControlledByPcgPresent, "Primary", "Metastatic")
#dev.new()
#VlnPlot(seuratObjData, features = pcgGenesDifExpPrimaryAndMetast)
#########################################################################################
# Differential Expression Analysis btw Primary & Normal
difExpGenesPrimaryAndNormal = difExpBetweenTypes(seuratObjData,geneSymbols,testUsing,"Primary","Normal")
pcgGenesDifExpPrimaryAndNormal = difExpPcgGenes(difExpGenesPrimaryAndNormal,pcglst,"PcG mechanism","Primary","Normal")
cat(paste(length(difExpGenesPrimaryAndNormal),"genes differentially expressed between primary & normal"))
difExpPcgMechAndTargetGenes(pcgGenesDifExpPrimaryAndNormal, genesControlledByPcgPresent, "Primary", "Normal")
#dev.new()
#VlnPlot(seuratObjData, features = pcgGenesDifExpPrimaryAndNormal)
# #########################################################################################
# # Differential Expression Analysis btw Metastatic & Normal (either Lymph if head & neck OR normal if ovarian)
if (ovarianCancerBool == TRUE) {
  difExpGenesMetastAndNormal = difExpBetweenTypes(seuratObjData,geneSymbols,testUsing,"Normal","Metastatic")
  pcgGenesDifExpMetastAndNormal = difExpPcgGenes(difExpGenesMetastAndNormal,pcglst,"PcG mechanism","Normal","Metastatic")
  print(paste(length(difExpGenesMetastAndNormal),"genes differentially expressed between Normal & Metastatic"))
  difExpPcgMechAndTargetGenes(pcgGenesDifExpMetastAndNormal, genesControlledByPcgPresent, "Normal", "Metastatic")
  #dev.new()
  #VlnPlot(seuratObjData, features = pcgGenesDifExpMetastAndNormal)
} else {
  difExpGenesMetastAndNormal = difExpBetweenTypes(seuratObjData,geneSymbols,testUsing,"Normal","Metastatic")
  pcgGenesDifExpMetastAndNormal = difExpPcgGenes(difExpGenesMetastAndNormal,pcglst,"PcG mechanism","Normal","Metastatic")
  print(paste(length(difExpGenesMetastAndNormal),"genes differentially expressed between Normal & Metastatic"))
  difExpPcgMechAndTargetGenes(pcgGenesDifExpMetastAndNormal, genesControlledByPcgPresent, "Normal", "Metastatic")

  difExpGenesMetastAndLymph = difExpBetweenTypes(seuratObjData,geneSymbols,testUsing,"Lymph","Metastatic")
  pcgGenesDifExpMetastAndLymph = difExpPcgGenes(difExpGenesMetastAndLymph,pcglst,"PcG mechanism","Lymph","Metastatic")
  print(paste(length(difExpGenesMetastAndLymph),"genes differentially expressed between Lymph & Metastatic"))
  difExpPcgMechAndTargetGenes(pcgGenesDifExpMetastAndLymph, genesControlledByPcgPresent, "Lymph", "Metastatic")

  difExpGenesLymphAndNormal = difExpBetweenTypes(seuratObjData,geneSymbols,testUsing,"Normal","Lymph")
  pcgGenesDifExpLymphAndNormal = difExpPcgGenes(difExpGenesLymphAndNormal,pcglst,"PcG mechanism","Normal","Lymph")
  print(paste(length(difExpGenesLymphAndNormal),"genes differentially expressed between Normal & Lymph"))
  difExpPcgMechAndTargetGenes(pcgGenesDifExpLymphAndNormal, genesControlledByPcgPresent, "Normal", "Lymph")
}

#########################################################################################
# Run PCA analysis on genes differentially expressed btw primary and metastatic cells:
# remove pcg genes that are NOT targets of pcg:
pcgGenesDifExpPrimaryMetastNotUnderControl = pcgGenesDifExpPrimaryAndMetast[! pcgGenesDifExpPrimaryAndMetast %in% genesControlledByPcgPresent]
geneDifExpPrimaryMetastatNotIncludingPcg = genesDifExpPrimaryAndMetast[! genesDifExpPrimaryAndMetast %in% pcgGenesDifExpPrimaryMetastNotUnderControl]
## If want to scale and center subset of geneDifExpPrimaryMetastatNotIncludingPcg genes before running PCA:
# seuratObjData <- ScaleData(seuratObjData, features = geneDifExpPrimaryMetastatNotIncludingPcg, do.scale = TRUE, do.center = TRUE)
# seuratObjData <- RunPCA(seuratObjData, features = geneDifExpPrimaryMetastatNotIncludingPcg, verbose = FALSE)
# # Plot PCA:
# dev.new()
# DimPlot(seuratObjData, reduction = "pca")
# #title(main = paste("PCA for",titleToAdd,"when use",matchingTitle,scalingTitle,"the data for genes DE btw primary & metastatic cells"))
# # Plot information about PCA:
# dev.new()
# VizDimLoadings(seuratObjData, dims = 1:2, reduction = "pca")
# dev.new()
# ElbowPlot(seuratObjData)
# dev.new()
# FeaturePlot(seuratObjData, reduction = "pca", features = c("CDH13", "PCDHB15"))
# dev.new()
# RidgePlot(seuratObjData, features = c("CDH13", "PCDHB15", "RAP1B", "CDH23"))
# dev.new()
# VlnPlot(seuratObjData, features = c("CDH13", "PCDHB15"))
#########################################################################################
# Remove multimodal genes from analysis:
# multimodalGenes = all.genes[modesOfGenes %in% uniqueModesOfGenes[uniqueModesOfGenes > 0]]
# nonMultimodalGenes = geneDifExpPrimaryMetastatNotIncludingPcg[! geneDifExpPrimaryMetastatNotIncludingPcg %in% multimodalGenes]
# seuratObjData <- RunPCA(seuratObjData, features = nonMultimodalGenes, verbose = FALSE)
# # Plot PCA:
# dev.new()
# DimPlot(seuratObjData, reduction = "pca")
#########################################################################################
# # Run PCA analysis on all genes:
# seuratObjData <- RunPCA(seuratObjData, features = all.genes, verbose = FALSE)
# # Plot PCA:
# dev.new()
# DimPlot(seuratObjData, reduction = "pca")
# # Look at how specific gene's expression changes in each cell on the PCA plot:
# dev.new()
# FeaturePlot(seuratObjData, features = pcgGenesDifExpPrimaryMetastNotUnderControl)
#########################################################################################
# # Run UMAP and plot:
# seuratObjData <- RunUMAP(seuratObjData, dims = 1:3)
# dev.new()
# DimPlot(seuratObjData, reduction = "umap")
# dev.new()
# FeaturePlot(seuratObjData, reduction = "umap", features = c("CDH13", "PCDHB15"))
###################################################################################################
# Pull out normalized gene expression matrix data from Seurat Object to run rest of grant analysis:
normalizedData = seuratObjData@assays$RNA@scale.data
# Set the column names back to either metastatic, primary, lymph, or normal --> so changing the unique
# names (ike metast.1, metast.2, etc.) that seurat package gave them back to our names:
colnames(normalizedData) <- labelsToUse
# Pull out gene expression matrix for only those genes that were differentially expressed between
# primary and metastatic that excludes the PcG mechanism genes but not those controlled by PcG mechanism:
differences = normalizedData[geneDifExpPrimaryMetastatNotIncludingPcg,]
if (ovarianCancerBool == TRUE) {
  levels<-c('Normal','Primary','Metastatic')
  pts<- c(1,1,1)
  colorsToUse = c("red","green","blue")
  titleToUseToSave = "ovarian"
} else {
  levels<-c('Lymph','Normal','Primary','Metastatic')
  pts<- c(16,1,1,1)
  colorsToUse = c("black","red","green","blue")
  titleToUseToSave = "headAndNeck"
}
f<- factor(colnames(differences),levels)
# Plot PCA to test that same as Seurat Package PCA --> is somewhat different rotation on PC1 but
# otherwise the same. This version is actually better to use for grant!
pca<-prcomp(differences,center=TRUE,scale=TRUE)
if (savePlotsAsPdfBool == TRUE){
  pdf(paste(titleToUseToSave,"PCA.pdf",sep=""), height=5, width=7)
  plot(pca$rotation[,1],pca$rotation[,2],col=colorsToUse[f], pch=pts[f],xlab="Principal Component 1",ylab="Principal Component 2")
  legend("bottomleft",NULL,levels(f),col=colorsToUse,pch=pts)
  dev.off()
} else {
  dev.new()
  plot(pca$rotation[,1],pca$rotation[,2],col=colorsToUse[f], pch=pts[f],xlab="Principal Component 1",ylab="Principal Component 2")
  legend("bottomleft",NULL,levels(f),col=colorsToUse,pch=pts)
}
#########################################################################################
# Run mapper analysis:
# top10pca<-pca$rotation[,1:5]
# library(TDAmapper)
# if (ovarianCancerBool==TRUE){
#   mapper1<-mapper1D(distance_matrix = dist(top10pca),
#                     filter_values=pca$rotation[,1],
#                     num_intervals=25,
#                     percent_overlap = 50,
#                     num_bins_when_clustering = 40)
# } else {
#   mapper1<-mapper1D(distance_matrix = dist(top10pca),
#                     filter_values=pca$rotation[,1],
#                     num_intervals=30,
#                     percent_overlap = 50,
#                     num_bins_when_clustering = 40)
# }
# library(igraph)
# g1 <- graph.adjacency(mapper1$adjacency, mode="undirected")
# for (i in 1:mapper1$num_vertices) {
#   if (ovarianCancerBool==TRUE){
#     curpts<-mapper1$points_in_vertex[[i]]
#     numm<-sum(labelsToUse[curpts]=='Metastatic')
#     nump<-sum(labelsToUse[curpts]=='Primary')
#     numn<-sum(labelsToUse[curpts]=='Normal')
#     clrs<-c("blue","green","red")
#     j<-which.max(c(numm,nump,numn))
#     V(g1)$color[i]<-clrs[j]
#   } else {
#     curpts<-mapper1$points_in_vertex[[i]]
#     numm<-sum(labelsToUse[curpts]=='Metastatic')
#     nump<-sum(labelsToUse[curpts]=='Primary')
#     numn<-sum(labelsToUse[curpts]=='Normal')
#     numl<-sum(labelsToUse[curpts]=='Lymph')
#     clrs<-c("blue","green","red","black")
#     j<-which.max(c(numm,nump,numn,numl))
#     V(g1)$color[i]<-clrs[j]
#   }
# }
# if (savePlotsAsPdfBool == TRUE){
#   pdf(paste(titleToUseToSave,"Mapper.pdf",sep=""), height=5, width=7)
#   plot(g1)
#   legend("topright",c("Normal","Primary","Metastatic"),fill=c("red","green","blue"),cex=0.5)
#   dev.off()
# } else {
#   dev.new()
#   plot(g1)
#   legend("topright",c("Normal","Primary","Metastatic"),fill=c("red","green","blue"),cex=0.5)
# }

# # Mapper 2D function (using PCA)
# mapper2 <- mapper2D(distance_matrix=dist(top10pca),
#                     filter_values=list(pca$rotation[,1],pca$rotation[,2]),
#                     num_intervals=c(7,7),
#                     percent_overlap=40,
#                     num_bins_when_clustering=15)
# library(igraph)
# g2 <- graph.adjacency(mapper2$adjacency, mode="undirected")
# for (i in 1:mapper2$num_vertices) {
#   if (ovarianCancerBool==TRUE){
#     curpts<-mapper2$points_in_vertex[[i]]
#     numm<-sum(labelsToUse[curpts]=='Metastatic')
#     nump<-sum(labelsToUse[curpts]=='Primary')
#     numn<-sum(labelsToUse[curpts]=='Normal')
#     clrs<-c("blue","green","red")
#     j<-which.max(c(numm,nump,numn))
#     V(g2)$color[i]<-clrs[j]
#   } else {
#     curpts<-mapper2$points_in_vertex[[i]]
#     numm<-sum(labelsToUse[curpts]=='Metastatic')
#     nump<-sum(labelsToUse[curpts]=='Primary')
#     numn<-sum(labelsToUse[curpts]=='Normal')
#     numl<-sum(labelsToUse[curpts]=='Lymph')
#     clrs<-c("blue","green","red","black")
#     j<-which.max(c(numm,nump,numn,numl))
#     V(g2)$color[i]<-clrs[j]
#   }
# }
# dev.new()
# plot(g2)
# legend("topright",c("Normal","Primary","Metastatic"),fill=c("red","green","blue"),cex=0.8)
#########################################################################################
# Cluster Primary & Metastatic into different subclusters based on PCA results,
# then redo DE analysis:
if (subclusteringBool == TRUE) {
  if (ovarianCancerBool == TRUE) {
    levelsSubClustering <- c('Normal','Primary1','Primary2','Metastatic')
    #levelsSubClustering <- c('Normal','Primary','Metastatic1','Metastatic2')
    colors <- c('red','cyan','green','blue')
    ptsVec <- c(16,1,1,1)
  } else {
    levelsSubClustering <- c('Lymph','Normal','Primary1','Primary2','Metastatic')
    #levelsSubClustering <- c('Lymph','Normal','Primary','Metastatic1','Metastatic2')
    colors <- c('black','red','cyan','green','blue')
    ptsVec <- c(16,1,1,1,1)
  }
  # Set cell type identities to use for labeling cells in analyses below:
  y <- seuratObjData@assays$RNA@counts@Dimnames[[2]]
  if (ovarianCancerBool == TRUE) {
    # For Ovarian cancer Data:
    labelsToUse<-as.character(sapply(y,convertLabelOvarian))
  } else {
    # For Head & Neck cancer data:
    labelsToUse<-as.character(sapply(y,convertLabel))
  }
  # Add labels as CellType in Seurat Obj:
  seuratObjData$CellType <- labelsToUse
  Idents(seuratObjData) <- "CellType"
  mylabelsSubClustering <- labelsToUse
  differences = normalizedData[geneDifExpPrimaryMetastatNotIncludingPcg,]
  pca<-prcomp(differences,center=TRUE,scale=TRUE)
  pcaLoadings = seuratObjData@reductions$pca@feature.loadings
  mylabelsSubClustering[labelsToUse=='Primary' & pca$rotation[,2]>=0] <-'Primary1'
  mylabelsSubClustering[labelsToUse=='Primary' & pca$rotation[,2]<0] <-'Primary2'
  #mylabelsSubClustering[labelsToUse=='Metastatic' & pca$rotation[,1]>=(-0.02)] <-'Metastatic1'
  #mylabelsSubClustering[labelsToUse=='Metastatic' & pca$rotation[,1]<(-0.02)] <-'Metastatic2'
  factorSubClustering <- factor(mylabelsSubClustering,levelsSubClustering)
  designSubclustering <- model.matrix(~0+factorSubClustering)
  colnames(designSubclustering)=levelsSubClustering
  if (savePlotsAsPdfBool == TRUE) {
    pdf(paste(titleToUseToSave,"subclusteringPCA.pdf",sep=""), height=5, width=7)
    plot(pca$rotation[,1],pca$rotation[,2],col=colors[factorSubClustering], main="PCA with Primary Spilt", pch = ptsVec[factorSubClustering])
    legend("topleft",NULL,levels(factorSubClustering),col=rep(colors),pch = ptsVec)
    dev.off()
  } else {
    dev.new()
    plot(pca$rotation[,1],pca$rotation[,2],col=colors[factorSubClustering], main="PCA with Primary Spilt", pch = ptsVec[factorSubClustering])
    legend("topleft",NULL,levels(factorSubClustering),col=rep(colors),pch = ptsVec)
  }
  ###################################################################################################
  ### UMAP with subclustering of primary cells: ###
  # Add subclustering labels in Seurat Obj:
  seuratDataSubcluster <- seuratObjData
  seuratDataSubcluster$CellType <- mylabelsSubClustering
  Idents(seuratDataSubcluster) <- "CellType"
  # Run UMAP and plot:
  seuratDataSubcluster <- RunUMAP(seuratDataSubcluster, dims = 1:3)
  if (ovarianCancerBool==TRUE){
    #pdf(paste(titleToUseToSave,"subclusteringUMAP.pdf",sep=""), height=5, width=7)
    dev.new()
    DimPlot(seuratDataSubcluster, reduction = "umap", cols = c("green","cyan", "blue","red"), pt.size = 0.8)
    #dev.off()
  } else {
    #pdf(paste(titleToUseToSave,"subclusteringUMAP.pdf",sep=""), height=5, width=7)
    dev.new()
    DimPlot(seuratDataSubcluster, reduction = "umap", cols = c("black", "red","blue", "cyan","green"), pt.size = 0.8)
    #dev.off()
  }
  
  
  ### Differential expression analysis after subclustering: ###
  seuratDataSubcluster <- seuratObjData
  seuratDataSubcluster$CellType <- mylabelsSubClustering
  Idents(seuratDataSubcluster) <- "CellType"
  
  testUsing = "wilcox"
  genesDifExpPrimary2AndMetastSubclustering = difExpBetweenTypes(seuratDataSubcluster,geneSymbols,testUsing,"Primary2","Metastatic")
  pcgGenesDifExpPrimary2AndMetastSubclustering = difExpPcgGenes(genesDifExpPrimary2AndMetastSubclustering, pcglst, "PcG mechanism", "Primary2","Metastatic")
  cat(paste(length(genesDifExpPrimary2AndMetastSubclustering),"genes differentially expressed between primary2 subcluster & metastatic"))
  
  genesDifExpPrimary1AndMetastSubclustering = difExpBetweenTypes(seuratDataSubcluster,geneSymbols,testUsing,"Primary1","Metastatic")
  pcgGenesDifExpPrimary1AndMetastSubclustering = difExpPcgGenes(genesDifExpPrimary1AndMetastSubclustering, pcglst, "PcG mechanism", "Primary1","Metastatic")
  cat(paste(length(genesDifExpPrimary1AndMetastSubclustering),"genes differentially expressed between primary1 subcluster & metastatic"))
  
  primary2Metast <- setdiff(pcgGenesDifExpPrimary1AndMetastSubclustering,pcgGenesDifExpPrimary2AndMetastSubclustering)
  dev.new()
  RidgePlot(seuratDataSubcluster, features = primary2Metast)
  dev.new()
  VlnPlot(seuratDataSubcluster, features = primary2Metast, slot = "counts", log = TRUE)
  dev.new()
  FeaturePlot(seuratObjData, features = primary2Metast)
  
  difExpGenesPrimary2AndNormalSubclustering = difExpBetweenTypes(seuratDataSubcluster,geneSymbols,testUsing,"Primary2","Normal")
  pcgGenesDifExpPrimary2AndNormalSubclustering = difExpPcgGenes(difExpGenesPrimary2AndNormalSubclustering,pcglst,"PcG mechanism","Primary2","Normal")
  cat(paste(length(difExpGenesPrimary2AndNormalSubclustering),"genes differentially expressed between primary2 subcluster & normal"))
  difExpPcgMechAndTargetGenes(pcgGenesDifExpPrimary2AndNormalSubclustering, genesControlledByPcgPresent, "Primary2", "Normal")
  
  difExpGenesPrimary1AndNormalSubclustering = difExpBetweenTypes(seuratDataSubcluster,geneSymbols,testUsing,"Primary1","Normal")
  pcgGenesDifExpPrimary1AndNormalSubclustering = difExpPcgGenes(difExpGenesPrimary1AndNormalSubclustering,pcglst,"PcG mechanism","Primary1","Normal")
  cat(paste(length(difExpGenesPrimary1AndNormalSubclustering),"genes differentially expressed between primary1 subcluster & normal"))
  difExpPcgMechAndTargetGenes(pcgGenesDifExpPrimary1AndNormalSubclustering, genesControlledByPcgPresent, "Primary1", "Normal")
  
  primaryNormal <- setdiff(pcgGenesDifExpPrimary1AndNormalSubclustering,pcgGenesDifExpPrimary2AndNormalSubclustering)
  dev.new()
  RidgePlot(seuratDataSubcluster, features = primaryNormal)
  dev.new()
  VlnPlot(seuratDataSubcluster, features = primaryNormal[1:3], slot = "counts", log = TRUE)
  dev.new()
  FeaturePlot(seuratObjData, features = primaryNormal)
}

###################################################################################################
# Aim 2 H2 quantify distance from metastatic to
# normal lymph vs primary and normal primary
#####################################################
### Quantify distance by using the log fold change from
# the differential expression analysis between the
# different pairs and find the mean & standard deviation
# of these log fold change values to compare for the different pairs:
foldChangeAnalysisFunc <- function(seuratObjectData, identity1String, identity2String){
  differentialExpressionMetastaticNormal <- FindMarkers(seuratObjectData, ident.1 = identity1String, ident.2 = identity2String, test.use = "wilcox")
  return(differentialExpressionMetastaticNormal)
}
differentialExpressedGenesAnalysisFunc <- function(differentialExpressionMetastaticNormal, identity1String, identity2String) {
  diffindsMetastaticNormal <- sapply(geneSymbols[1:length(geneSymbols)],function(x) {any(rownames(differentialExpressionMetastaticNormal)==x)})
  difExpGenesMetastaticNormal <- geneSymbols[diffindsMetastaticNormal]
  pcgdiffsMetastaticNormal<- difExpGenesMetastaticNormal %in% polycombGenesPresent
  trithoraxdiffsMetastaticNormal <- difExpGenesMetastaticNormal %in% trithoraxGenesPresent
  controldiffsMetastaticNormal <- difExpGenesMetastaticNormal %in% genesControlledByPcgPresent
  pcgMechanismdiffs <- difExpGenesMetastaticNormal %in% c(polycombGenesPresent,trithoraxGenesPresent)
  cat(length(difExpGenesMetastaticNormal[pcgdiffsMetastaticNormal]), "Polycomb genes dif. expressed in", identity1String, "&", identity2String, difExpGenesMetastaticNormal[pcgdiffsMetastaticNormal],
      "\nSum of log fold change:",sum(differentialExpressionMetastaticNormal$avg_logFC[pcgdiffsMetastaticNormal]),
      "\nStd of log fold change:",sd(differentialExpressionMetastaticNormal$avg_logFC[pcgdiffsMetastaticNormal]),
      "\n",length(difExpGenesMetastaticNormal[trithoraxdiffsMetastaticNormal]), "Trithorax genes dif. expressed in", identity1String, "&", identity2String, difExpGenesMetastaticNormal[trithoraxdiffsMetastaticNormal],
      "\nSum of log fold change:",sum(differentialExpressionMetastaticNormal$avg_logFC[trithoraxdiffsMetastaticNormal]),
      "\nStd of log fold change:",sd(differentialExpressionMetastaticNormal$avg_logFC[trithoraxdiffsMetastaticNormal]),
      "\n",length(difExpGenesMetastaticNormal[controldiffsMetastaticNormal]), "PcG controlled genes dif. expressed in", identity1String, "&", identity2String, difExpGenesMetastaticNormal[controldiffsMetastaticNormal],
      "\nSum of log fold change:",sum(differentialExpressionMetastaticNormal$avg_logFC[controldiffsMetastaticNormal]),
      "\nStd of log fold change:",sd(differentialExpressionMetastaticNormal$avg_logFC[controldiffsMetastaticNormal]),
      "\n",length(difExpGenesMetastaticNormal[pcgMechanismdiffs]), "Polycomb & trithorax genes dif. expressed in", identity1String, "&", identity2String, difExpGenesMetastaticNormal[pcgMechanismdiffs],
      "\nSum of log fold change:",sum(differentialExpressionMetastaticNormal$avg_logFC[pcgMechanismdiffs]),
      "\nStd of log fold change:",sd(differentialExpressionMetastaticNormal$avg_logFC[pcgMechanismdiffs]))
  return(difExpGenesMetastaticNormal)
}
# Primary Vs Normal Results:
differentialExpressionPrimaryNormal <- foldChangeAnalysisFunc(seuratObjData, "Primary", "Normal") # contains differential expression analysis information like log fold change for those genes differentially expressed
#write.csv(differentialExpressionPrimaryNormal, paste(titleToUseToSave,"PrimaryVsNormalDiffExpResultsDf.csv", sep=""))
#differentialExpressionPrimaryNormal <- read.csv(paste(titleToUseToSave,"PrimaryVsNormalDiffExpResultsDf.csv", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
difExpGenesPrimaryNormal <- differentialExpressedGenesAnalysisFunc(differentialExpressionPrimaryNormal, "Primary", "Normal")
pcgdiffsPrimaryNormal<- difExpGenesPrimaryNormal[difExpGenesPrimaryNormal %in% c(polycombGenesPresent,trithoraxGenesPresent)]

# Metastatic vs Normal:
differentialExpressionMetastaticNormal <- foldChangeAnalysisFunc(seuratObjData, "Metastatic", "Normal")
#write.csv(differentialExpressionMetastaticNormal, paste(titleToUseToSave,"MetastaticVsNormalDiffExpResultsDf.csv", sep=""))
#differentialExpressionMetastaticNormal <- read.csv(paste(titleToUseToSave,"MetastaticVsNormalDiffExpResultsDf.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE, row.names = 1)
difExpGenesMetastaticNormal <- differentialExpressedGenesAnalysisFunc(differentialExpressionMetastaticNormal, "Metastatic", "Normal")
pcgdiffsMetastaticNormal<- difExpGenesMetastaticNormal[difExpGenesMetastaticNormal %in% c(polycombGenesPresent,trithoraxGenesPresent)]

if (ovarianCancerBool == FALSE){
  # Metastatic vs Lymph:
  differentialExpressionMetastaticLymph <- foldChangeAnalysisFunc(seuratObjData, "Metastatic", "Lymph")
  #write.csv(differentialExpressionMetastaticLymph, paste(titleToUseToSave,"MetastaticVsLymphDiffExpResultsDf.csv", sep=""))
  #differentialExpressionMetastaticLymph <- read.csv(paste(titleToUseToSave,"MetastaticVsLymphDiffExpResultsDf.csv", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
  difExpGenesMetastaticLymph <- differentialExpressedGenesAnalysisFunc(differentialExpressionMetastaticLymph, "Metastatic", "Lymph")
  pcgdiffsMetastaticLymph<- difExpGenesMetastaticLymph[difExpGenesMetastaticLymph %in% c(polycombGenesPresent,trithoraxGenesPresent)]
  
  MetastaticLymphLogFoldChangeDf <- as.matrix(differentialExpressionMetastaticLymph$avg_logFC)
  rownames(MetastaticLymphLogFoldChangeDf) <- difExpGenesMetastaticLymph
  pcgMechGenesMetastaticLymphLogFoldChangeDf <- MetastaticLymphLogFoldChangeDf[difExpGenesMetastaticLymph %in% c(polycombGenesPresent,trithoraxGenesPresent),]
}

# Results above when subset data into primary1 and primary2:
# Primary2 Vs Normal Results:
differentialExpressionPrimary2Normal <- foldChangeAnalysisFunc(seuratDataSubcluster, "Primary2", "Normal")
#write.csv(differentialExpressionPrimary2Normal, paste(titleToUseToSave,"Primary2VsNormalDiffExpResultsDf.csv", sep=""))
#differentialExpressionPrimary2Normal <- read.csv(paste(titleToUseToSave,"Primary2VsNormalDiffExpResultsDf.csv", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
difExpGenesPrimary2Normal <- differentialExpressedGenesAnalysisFunc(differentialExpressionPrimary2Normal, "Primary2", "Normal")
pcgdiffsPrimary2Normal<- difExpGenesPrimary2Normal[difExpGenesPrimary2Normal %in% c(polycombGenesPresent,trithoraxGenesPresent)]
# Primary1 Vs Normal Results:
differentialExpressionPrimary1Normal <- foldChangeAnalysisFunc(seuratDataSubcluster, "Primary1", "Normal")
#write.csv(differentialExpressionPrimary1Normal, paste(titleToUseToSave,"Primary1VsNormalDiffExpResultsDf.csv", sep=""))
#differentialExpressionPrimary1Normal <- read.csv(paste(titleToUseToSave,"Primary1VsNormalDiffExpResultsDf.csv", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
difExpGenesPrimary1Normal <- differentialExpressedGenesAnalysisFunc(differentialExpressionPrimary1Normal, "Primary1", "Normal")
pcgdiffsPrimary1Normal<- difExpGenesPrimary1Normal[difExpGenesPrimary1Normal %in% c(polycombGenesPresent,trithoraxGenesPresent)]


### Find significance of differences for PcG genes ###
allPcGMechanismGenes <- c(polycombGenesPresent,trithoraxGenesPresent)
numPcgMechanismOnlyGenes <- length(allPcGMechanismGenes)

MetastaticNormalLogFoldChangeDf <- as.matrix(differentialExpressionMetastaticNormal$avg_logFC)
rownames(MetastaticNormalLogFoldChangeDf) <- difExpGenesMetastaticNormal
pcgMechGenesMetastaticNormalLogFoldChangeDf <- MetastaticNormalLogFoldChangeDf[difExpGenesMetastaticNormal %in% c(polycombGenesPresent,trithoraxGenesPresent),]

PrimaryNormalLogFoldChangeDf <- as.matrix(differentialExpressionPrimaryNormal$avg_logFC)
rownames(PrimaryNormalLogFoldChangeDf) <- difExpGenesPrimaryNormal
pcgMechGenesPrimaryNormalLogFoldChangeDf <- PrimaryNormalLogFoldChangeDf[difExpGenesPrimaryNormal %in% c(polycombGenesPresent,trithoraxGenesPresent),]

dfForSignificance = data.frame(matrix(0, 2, numPcgMechanismOnlyGenes))
colnames(dfForSignificance) <- allPcGMechanismGenes
for (i in 1:length(pcgdiffsPrimaryNormal)) {
  dfForSignificance[1,pcgdiffsPrimaryNormal[i]] <- pcgMechGenesPrimaryNormalLogFoldChangeDf[pcgdiffsPrimaryNormal[i]]
}
for (i in 1:length(pcgdiffsMetastaticNormal)) {
  dfForSignificance[2,pcgdiffsMetastaticNormal[i]] <- pcgMechGenesMetastaticNormalLogFoldChangeDf[pcgdiffsMetastaticNormal[i]]
}
# Build random log fold change vectors to compare:
differencesOfRandomSample <- rep(0,1000)
for (k in 1:1000) {
  for (j in 1:numPcgMechanismOnlyGenes) {
    dfForSignificance[,j]<-sample(dfForSignificance[,j])
  }
  sumOfRandomSample <- rowSums(dfForSignificance)
  differencesOfRandomSample[k] <- sumOfRandomSample[1] - sumOfRandomSample[2]
}
dev.new()
hist(differencesOfRandomSample)
actualDistance <- (sum(pcgMechGenesPrimaryNormalLogFoldChangeDf)-sum(pcgMechGenesMetastaticNormalLogFoldChangeDf))
abline(v = actualDistance, col="red")

library("fitdistrplus")
fit <- fitdist(differencesOfRandomSample[!is.na(differencesOfRandomSample)], "norm")
pValHist <- pnorm(q = actualDistance, mean= fit$estimate[['mean']], sd= fit$estimate[['sd']])
cat(paste("P-value for Dif. in Log Fold Changes btw Primary & Metastatic:", 1-pValHist))



### Quantifying distance from metastatic to normal
# lymph vs primary to normal primary for just PcG
# genes after pre-processing:
# only use PcG genes that control others (so not including genes that are under control of PcG mechanism)
# When only 1 PcG mechanism gene differentially expressed btw primary and metastatic, then will not be able to calculate distance...
pcgGenesDifExpPrimaryMetastUnderControl = pcgGenesDifExpPrimaryAndMetast[pcgGenesDifExpPrimaryAndMetast %in% genesControlledByPcgPresent]
pcgGenesDifExpPrimaryMetastNotUnderControl = geneSymbols[! geneSymbols %in% genesControlledByPcgPresent]
pcgdiffs<-normalizedData[geneSymbols[geneSymbols %in% c(polycombGenesPresent,trithoraxGenesPresent)],]
if (ovarianCancerBool==TRUE) {
  pcgdatm<-pcgdiffs[,labelsToUse=='Metastatic']
  pcgdatp<-pcgdiffs[,labelsToUse=='Primary']
  pcgdatn<-pcgdiffs[,labelsToUse=='Normal']
  pcgmeansm<-rowMeans(pcgdatm)
  pcgmeansp<-rowMeans(pcgdatp)
  pcgmeansn<-rowMeans(pcgdatn)
  pcgdiff1<-pcgmeansp-pcgmeansn
  pcgdiff2<-pcgmeansm-pcgmeansn
  cat("Ovarian: Distance btw primary & normal for PcG genes:",dist(rbind(pcgmeansp,pcgmeansn)),
      "\nDistance btw metastatic & normal for PcG genes:",dist(rbind(pcgmeansm,pcgmeansn)),"\n")
  cat("Ovarian: Mean of differences btw primary & normal for PcG genes:", mean(pcgdiff1),
      "\nMean of differences btw metastatic & normal for PcG genes:",mean(pcgdiff2),"\n")
} else {
  pcgdatm<-pcgdiffs[,labelsToUse=='Metastatic']
  pcgdatp<-pcgdiffs[,labelsToUse=='Primary']
  pcgdatn<-pcgdiffs[,labelsToUse=='Normal']
  pcgdatl<-pcgdiffs[,labelsToUse=='Lymph']
  pcgmeansm<-rowMeans(pcgdatm)
  pcgmeansp<-rowMeans(pcgdatp)
  pcgmeansn<-rowMeans(pcgdatn)
  pcgmeansl<-rowMeans(pcgdatl)
  pcgdiff1<-pcgmeansp-pcgmeansn
  pcgdiff2<-pcgmeansm-pcgmeansl
  cat("Head&Neck: Distance btw primary & normal for PcG genes:",dist(rbind(pcgmeansp,pcgmeansn)),
      "\nDistance btw metastatic & normal for PcG genes:",dist(rbind(pcgmeansm,pcgmeansl)),"\n")
  cat("Head&Neck: Mean of differences btw primary & normal for PcG genes:", mean(pcgdiff1), "Std:",sd(pcgdiff1),
      "\nMean of differences btw metastatic & normal for PcG genes:",mean(pcgdiff2), "Std:",sd(pcgdiff2))
}
## Gene by gene differences for PcG genes so for H1:
if (ovarianCancerBool==TRUE) {
  diff1<-pcgmeansp-pcgmeansn
  diff2<-pcgmeansm-pcgmeansn
  colorsGeneDiffs <- rep("black", length(diff2))
  colorsGeneDiffs[abs(diff2) <= abs(diff1)] <- "red"
  if (savePlotsAsPdfBool==TRUE){
    pdf(paste(titleToUseToSave,"GeneByGeneDifferences.pdf",sep=""), height=5, width=7)
    plot(diff1,diff2,main = "Gene by Gene Differences", xlab = "Distance btw primary & normal", ylab = "Distance btw metastatic and normal", col=colorsGeneDiffs)
    abline(0,1)
    abline(0,-1)
    abline(h=0)
    abline(v=0)
    dev.off()
  } else {
    dev.new()
    plot(diff1,diff2,main = "Gene by Gene Differences", xlab = "Distance btw primary & normal", ylab = "Distance btw metastatic and normal", col=colorsGeneDiffs)
    abline(0,1)
    abline(0,-1)
    abline(h=0)
    abline(v=0)
  }
  c(sqrt(sum(diff1^2)),sqrt(sum(diff2^2)),mean(abs(diff2)<abs(diff1)))
  print(paste("Ovarian: # of genes in red =",length(colorsGeneDiffs[colorsGeneDiffs=="red"]), "out of", length(colorsGeneDiffs)))
} else {
  diff1<-pcgmeansp-pcgmeansn
  diff2<-pcgmeansm-pcgmeansl
  colorsGeneDiffs <- rep("black", length(diff2))
  colorsGeneDiffs[abs(diff2) <= abs(diff1)] <- "red"
  if (savePlotsAsPdfBool==TRUE){
    pdf(paste(titleToUseToSave,"GeneByGeneDifferences.pdf",sep=""), height=5, width=7)
    plot(diff1,diff2,main = "Gene by Gene Differences", xlab = "Distance btw primary & normal", ylab = "Distance btw metastatic and normal", col=colorsGeneDiffs)
    abline(0,1)
    abline(0,-1)
    abline(h=0)
    abline(v=0)
    dev.off()
  } else {
    dev.new()
    plot(diff1,diff2,main = "Gene by Gene Differences", xlab = "Distance btw primary & normal", ylab = "Distance btw metastatic and normal", col=colorsGeneDiffs)
    abline(0,1)
    abline(0,-1)
    abline(h=0)
    abline(v=0)
  }
  c(sqrt(sum(diff1^2)),sqrt(sum(diff2^2)),mean(abs(diff2)<abs(diff1)))
  print(paste("Head&Neck: # of genes in black =",length(colorsGeneDiffs[colorsGeneDiffs=="black"]), "out of", length(colorsGeneDiffs)))
}


### Quantifying distance from metastatic to normal
# lymph vs primary to normal primary for genes
# differentially expressed btw metastatic and
# primary that are NOT PcG mechanism genes or multimodal genes:
restrictToDeGenesBool = FALSE
if (subclusteringBool == TRUE) {
  if (restrictToDeGenesBool == TRUE) {
    differences = normalizedData[genesDifExpPrimaryAndMetast,mylabelsSubClustering!="Primary1"]
  } else {
    differences = normalizedData[,mylabelsSubClustering!="Primary1"]
  }
  labelsToUse = labelsToUse[mylabelsSubClustering!="Primary1"]
} else {
   if (restrictToDeGenesBool == TRUE) {
     differences = normalizedData[genesDifExpPrimaryAndMetast,]
     if (ovarianCancerBool == TRUE) {
       # For Ovarian cancer Data:
       labelsToUse<-as.character(sapply(y,convertLabelOvarian))
     } else {
       # For Head & Neck cancer data:
       labelsToUse<-as.character(sapply(y,convertLabel))
     }
   } else {
    differences = normalizedData
    if (ovarianCancerBool == TRUE) {
      # For Ovarian cancer Data:
      labelsToUse<-as.character(sapply(y,convertLabelOvarian))
    } else {
      # For Head & Neck cancer data:
      labelsToUse<-as.character(sapply(y,convertLabel))
    }
   }
}
if (ovarianCancerBool==TRUE) {
  datm<-differences[,labelsToUse=='Metastatic']
  datp<-differences[,labelsToUse=='Primary']
  datn<-differences[,labelsToUse=='Normal']
  meansm<-rowMeans(datm)
  meansp<-rowMeans(datp)
  meansn<-rowMeans(datn)
  # Aim 2 H2 quantify distance results:
  if (restrictToDeGenesBool == TRUE) {
    cat("Ovarian: Distance btw primary & normal for DE genes btw primary & metastatic:",dist(rbind(meansp,meansn)),
        "\nDistance btw metastatic & normal for DE genes btw primary & metastatic:",dist(rbind(meansm,meansn)),"\n")
  } else {
    cat("Ovarian: Distance btw primary & normal for all genes:",dist(rbind(meansp,meansn)),
        "\nDistance btw metastatic & normal for all genes:",dist(rbind(meansm,meansn)),"\n")
  }
} else {
  datm<-differences[,labelsToUse=='Metastatic']
  datp<-differences[,labelsToUse=='Primary']
  datn<-differences[,labelsToUse=='Normal']
  datl<-differences[,labelsToUse=='Lymph']
  meansm<-rowMeans(datm)
  meansp<-rowMeans(datp)
  meansn<-rowMeans(datn)
  meansl<-rowMeans(datl)
  # Aim 2 H2 quantify distance results:
  if (restrictToDeGenesBool == TRUE) {
    cat("Head&Neck: Distance btw primary & normal for DE genes btw primary & metastatic:",dist(rbind(meansp,meansn)),
        "\nDistance btw metastatic & lymph normal for DE genes btw primary & metastatic:",dist(rbind(meansm,meansl)),"\n")
  } else {
    cat("Head&Neck: Distance btw primary & normal for all genes:",dist(rbind(meansp,meansn)),
        "\nDistance btw metastatic & lymph normal for all genes:",dist(rbind(meansm,meansl)),"\n")
  }
}


### Aim 2 H2 quantify gene by gene differences:
if (ovarianCancerBool==TRUE) {
  diff1<-meansp-meansn
  diff2<-meansm-meansn
  colorsGeneDiffs <- rep("black", length(diff2))
  colorsGeneDiffs[abs(diff2) <= abs(diff1)] <- "red"
  #pdf(paste(titleToUseToSave,"GeneByGeneDifferences.pdf",sep=""), height=5, width=7)
  dev.new()
  plot(diff1,diff2,main = "Gene by Gene Differences", xlab = "Distance btw primary & normal", ylab = "Distance btw metastatic and normal", col=colorsGeneDiffs)
  #plot(fulldiff1,fulldiff2)
  abline(0,1)
  abline(0,-1)
  abline(h=0)
  abline(v=0)
  #abline(lm(diff2~diff1))
  #dev.off()
  c(sqrt(sum(diff1^2)),sqrt(sum(diff2^2)),mean(abs(diff2)<abs(diff1)))
  print(paste("Ovarian: # of genes in red =",length(colorsGeneDiffs[colorsGeneDiffs=="red"]), "out of", length(colorsGeneDiffs)))
} else {
  diff1<-meansp-meansn
  diff2<-meansm-meansl
  colorsGeneDiffs <- rep("black", length(diff2))
  colorsGeneDiffs[abs(diff2) <= abs(diff1)] <- "red"
  #pdf(paste(titleToUseToSave,"GeneByGeneDifferences.pdf",sep=""), height=5, width=7)
  dev.new()
  plot(diff1,diff2,main = "Gene by Gene Differences", xlab = "Distance btw primary & normal", ylab = "Distance btw metastatic and lymph normal", col = colorsGeneDiffs)
  #plot(fulldiff1,fulldiff2)
  abline(0,1)
  abline(0,-1)
  abline(h=0)
  abline(v=0)
  #abline(lm(diff2~diff1))
  #dev.off()
  c(sqrt(sum(diff1^2)),sqrt(sum(diff2^2)),mean(abs(diff2)<abs(diff1)))
  print(paste("Head&Neck: # of genes in red =",length(colorsGeneDiffs[colorsGeneDiffs=="red"]), "out of", length(colorsGeneDiffs)))
}

### Differential expression for primary1 and primary2 groups separately:

############################################################################
### Aim 2 H3: Pliancy Z-score analysis
#######################################
restrictToDeGenesBool = TRUE
numPcgMechanismOnlyGenes = (length(pcgnames)-length(genesControlledByPcgPresent)) # number of pcg mechanism genes only after preprocessing, so not including genes controlled by PcG mechanism
pcgOnlyInds <- pcginds & !controlledByPcg
pcgOnlynames<-geneSymbols[pcgOnlyInds] # only PcG mechanism gene names and not genes controlled by PcG

# If want to run with only genes differentially expressed btw primary and metastatic:
diffinds = sapply(geneSymbols[1:length(geneSymbols)],function(x) {any(genesDifExpPrimaryAndMetast==x)})
diffnotpcginds<- diffinds & !pcgOnlyInds
pcgMechanismGenesForPliancyScore <- c(polycombGenesPresent,trithoraxGenesPresent)
numPcgMechanismOnlyGenes <- length(pcgMechanismGenesForPliancyScore)
if (restrictToDeGenesBool == TRUE){
  allGenesExceptPcgMechanism <- difExpGenesMetastaticNormal[! difExpGenesMetastaticNormal %in% pcgMechanismGenesForPliancyScore]
} else {
  allGenesExceptPcgMechanism <- geneSymbols[! geneSymbols %in% pcgMechanismGenesForPliancyScore]
}
mydata2<-normalizedData[allGenesExceptPcgMechanism,]
mydata2m<-mydata2[,labelsToUse=='Metastatic']
mydata2p<-mydata2[,labelsToUse=="Primary"]
mydata2n<-mydata2[,labelsToUse=="Normal" | labelsToUse=="Lymph"]

len<-length(allGenesExceptPcgMechanism)#len<-sum(diffnotpcginds)

### For when using the Pliancy z-score:
myscores<-mydata2m
for (i in 1:len) {
  myrowm<-as.numeric(mydata2m[i,])
  mymeanm<-mean(myrowm)
  mysdm<-sd(myrowm)
  # z-score of distribution across all cells in a given
  # category (like metastatic, as did for 2/4/19 grant):
  myzm<-(myrowm-mymeanm)/mysdm
  myrown<-as.numeric(mydata2n[i,])
  mymeann<-mean(myrown)
  mysign<-sign(mymeanm-mymeann)
  # cell's direction of movement toward normal:
  mysignedz<-myzm*mysign
  myscores[i,]<-mysignedz
}
# Overall pliancy z-score:
mymeanscores<-colMeans(myscores)


### For when using distance between metastatic and normal instead of pliancy z-score:
myscores<-rep(0,dim(mydata2m)[2])
meansn<-rowMeans(mydata2n)
for (i in 1:dim(mydata2m)[2]) { # loop through each metastatic cell to find the distance between each cell compared to mean of normal
  metastaticCellExpression <- as.numeric(mydata2m[,i])
  myscores[i] <- dist(rbind(metastaticCellExpression,meansn))
}
# Overall distances:
mymeanscores<-myscores


pcgdat<-normalizedData[pcgMechanismGenesForPliancyScore,labelsToUse=="Metastatic"]
pcgdat<-data.frame(cbind(t(pcgdat),mymeanscores))
colnames(pcgdat)<-c(pcgMechanismGenesForPliancyScore,"score")

# dataFrameKmeans <- pcgdat[,c("score",pcgGeneToLookAt)]
# kmeansPcG <- kmeans(dataFrameKmeans, centers = 2, nstart = 100)
# dev.new()
# plot(dataFrameKmeans,col=c("blue","green")[factorKmeans])
# Fit linear regression for overall pliancy
# score as a function of all (length(pcgnames)-3) PcG-like genes:
# mymodel<-lm(score~.,data=pcgdat[kmeansPcG$cluster==1,c("score",pcgGeneToLookAt)])
# to find R-squared value look at summary of mymodel
mymodel<-lm(score~.,data=pcgdat)
summary(mymodel)
# Pull out PcG mechanism genes that are significant:
pvalsOfCoefficients <- summary(mymodel)$coefficients[-1,4] # first element is the p-value for the intercept and rest are for the slopes of each PcG mechanism gene, so remove first element
significantCoefficientInds <- pvalsOfCoefficients <= 0.1
# PcG genes with positive coefficients, meaning
# decreased expression positively correlated
# with movement towards normal:
pcgMechanismSlopes <- mymodel$coefficients[-1] # remove first element that is the intercept estimate because we just want to look at estimates of the slope for each PcG mechanism gene
pcgMechanismSlopesSignificant <- pcgMechanismSlopes[significantCoefficientInds]# pull out only slope values that are significant
#pcgGenesPos<-pcgMechanismSlopesSignificant[pcgMechanismSlopesSignificant > 0]
pcgGenesPos<- mymodel$coefficients[mymodel$coefficients > 0]
cat("# PcG genes w/ positive coefficients:",length(pcgGenesPos),"out of",length(pcgMechanismSlopesSignificant),"with significant coefficients")


# for (i in 1:length(names(pcgMechanismSlopesSignificant))) {
#   pcgGeneToLookAt <- names(pcgMechanismSlopesSignificant)[i]
#   dev.new()
#   plot(pcgdat[,pcgGeneToLookAt],mymeanscores,xlab=pcgGeneToLookAt)
#   abline(mymodel$coefficients[1],mymodel$coefficients[pcgGeneToLookAt],col="red")
# }
# 
# for (i in 1:length(names(pcgMechanismSlopesSignificant))) {
#   pcgGeneToLookAt <- controlPos[i]
#   dev.new()
#   plot(normalizedData[pcgGeneToLookAt,labelsToUse=="Metastatic"],mymeanscores,xlab=pcgGeneToLookAt)
#   abline(mymodel$coefficients[1],mymodel$coefficients[pcgGeneToLookAt],col="red")
# }


### Bootstrapping for Aim 2 H3 results to find pval ###

# If using all genes except differentially expressed btw primary and metastatic and also pcg mechanism genes: 
# controlPos<-sample(which(!pcgOnlyInds & !diffinds), numPcgMechanismOnlyGenes)
# Else
controlPos<-sample(allGenesExceptPcgMechanism, numPcgMechanismOnlyGenes) 
controldat<-normalizedData[controlPos,labelsToUse=="Metastatic"]
controldat<-data.frame(cbind(t(controldat),mymeanscores))
# If using all genes except differentially expressed btw primary and metastatic and also pcg mechanism genes: 
# colnames(controldat)<-c(names(controlPos),"score")
# Else:
colnames(controldat)<-c(controlPos,"score") 
mymodel<-lm(score~.,data=controldat)

nposdist<-rep(0,1000)
for (i in 1:1000)
{
  
  # If using all genes except differentially expressed btw primary and metastatic and also pcg mechanism genes: 
  # controlPos<-sample(which(!pcgOnlyInds & !diffinds), numPcgMechanismOnlyGenes)
  # Else:
  controlPos<-sample(allGenesExceptPcgMechanism,numPcgMechanismOnlyGenes) 
  controldat<-normalizedData[controlPos,labelsToUse=="Metastatic"]
  controldat<-data.frame(cbind(t(controldat),mymeanscores))
  # If using all genes except differentially expressed btw primary and metastatic and also pcg mechanism genes: 
  # colnames(controldat)<-c(names(controlPos),"score")
  # Else:
  colnames(controldat)<-c(controlPos,"score") 
  controlmodel<-lm(score~.,data=controldat)
  #npositive<-sum(controlmodel$coefficients[2:(numPcgMechanismOnlyGenes+1)][summary(controlmodel)$coefficients[-1,4] <= 0.1]>0)
  npositive<-sum(controlmodel$coefficients[-1]>0)
  nposdist[i]<-npositive
}
if (ovarianCancerBool == TRUE) {
  ylimMax = 200
} else {
  ylimMax = 300
}
if (savePlotsAsPdfBool == TRUE){
  pdf(paste(titleToUseToSave,"Bootstrapping.pdf",sep=""), height=5, width=7)
  hist(nposdist, breaks=13,xlim=c(0,length(pcgnames)),ylim=c(0,ylimMax))
  abline(v=length(pcgGenesPos),col="red")
  dev.off()
} else {
  dev.new()
  hist(nposdist, breaks=13,xlim=c(0,length(pcgnames)),ylim=c(0,ylimMax))
  abline(v=length(pcgGenesPos),col="red")
}

cat(paste("Mean and std of positive distribution:", mean(nposdist[!is.na(nposdist)]), sd(nposdist[!is.na(nposdist)])))
library("fitdistrplus")
fit <- fitdist(nposdist[!is.na(nposdist)], "norm")
pValHist <- pnorm(q = length(pcgGenesPos), mean= fit$estimate[['mean']], sd= fit$estimate[['sd']])
cat(paste("P-value of # PcG genes w/ positive coefficients:", 1-pValHist))
#####################################################################################





########################################################################################################################
############ RUN PREPROCESSING JUST LIKE HEAD & NECK CANCER DATA DID IN THEIR PAPER AND IN GRANT FEB 2018 ###############
# Run step 2 of RNA-seq preprocessing and analysis (document on mlambros@mail.einstein.yu.edu
# google drive titled "RNA-sequecing analysis notes") in which we eliminate any column (single-cell)
# with more than 80% of its genes with zero expression.
# Filter out genes according to paper method also (gene must be expressed in at least 3 cells):
genesToGetRidOf = c()
for (i in 2:dim(allDataAllCellTypesOvarian)[1]) { # first row contains cell types so only start at second row
  genesExpInAllCells = allDataAllCellTypesOvarian[i,2:dim(allDataAllCellTypesOvarian)[2]] # first column is gene names so only start at second column that contains data
  if (sum(genesExpInAllCells!=0) < 3) {
    if (length(genesToGetRidOf)==0){
      genesToGetRidOf = i
    } else {
      genesToGetRidOf = c(genesToGetRidOf, i)
    }
  }
}
allDataAllCellTypesOvarian = allDataAllCellTypesOvarian[-(genesToGetRidOf),1:dim(allDataAllCellTypesOvarian)[2]]
numRows = dim(allDataAllCellTypesOvarian)[1]
numGenes = numRows - 2 # minus 1 when no sample ID's and minus 2 when sample ID's included
numColumns= dim(allDataAllCellTypesOvarian)[2]
numSingleCells = numColumns - 1
cellTypeNames = allDataAllCellTypesOvarian[2,2:numColumns]
difCellTypes = unique(matrix(cellTypeNames))
sampleIds = allDataAllCellTypesOvarian[1,2:numColumns]
geneExpressionAllCellTypes = allDataAllCellTypesOvarian[3:numRows,2:numColumns] # 2:numRows when no sample ID's and 3:numRows when sample ID's
geneExpressionMatrixAllCellTypes = matrix(as.numeric(unlist(geneExpressionAllCellTypes)), ncol=numSingleCells, byrow = FALSE)

# Eliminate any column (single-cell) with more than __% of its genes with zero expression
eliminateColsWithMoreThanPercentGenesWithZeroExpression <- function(allDataAllCellTypesOvarian, numSingleCells, geneExpressionMatrixAllCellTypes, cellTypeNames, percent, gzFileName, labelFileName){
  list80 = list()
  cellTypeNames80 = c("gene_symbol")
  list80[[1]] = allDataAllCellTypesOvarian[2:dim(allDataAllCellTypesOvarian)[1],1] # list of all genes (plus gene_symbol header name, but minus sample_id header name)
  count80 = 2
  for (i in 1:numSingleCells){
    if (quantile(geneExpressionMatrixAllCellTypes[,i], probs = c(percent))[[1]] > 0) {
      list80[[count80]] = c(cellTypeNames[i], geneExpressionMatrixAllCellTypes[,i])
      cellTypeNames80[count80] = cellTypeNames[i]
      count80 = count80 + 1
    }
  }
  geneExpression80MatrixAllCellTypes = matrix(as.character(unlist(list80)), ncol=(count80-1), byrow = FALSE)
  cat(paste("Number of cells keep after preprocessing is",dim(geneExpression80MatrixAllCellTypes)[2]))
  gz80 = gzfile(gzFileName, "w")
  write.table(geneExpression80MatrixAllCellTypes, gz80, sep="\t", quote=F, row.names = FALSE, col.names = FALSE)
  close(gz80)
  write.table(cellTypeNames80, labelFileName, row.names = FALSE, col.names = FALSE, sep="\t")
}
#### PREPROCESSING FROM PAPER:
# Eliminate any column (single-cell) exactly like in paper
# (kept cells with at least 200 genes expressed) so
# 100 - (200/dim(geneExpressionMatrixAllCellTypes)[1])*100
# so eliminate 98.94665% of its gene with zero expression
eliminateColsWithMoreThanPercentGenesWithZeroExpression(allDataAllCellTypesOvarian, numSingleCells, geneExpressionMatrixAllCellTypes, cellTypeNames, 0.9894665, "geneExpression98_946655MatrixAllCellTypesOvarian.gz", "labels98_946655")

## How we did Head and Neck Cancer data preprocessing
# Eliminate any column (single-cell) with more than 80%
# of its genes with zero expression
eliminateColsWithMoreThanPercentGenesWithZeroExpression(allDataAllCellTypesOvarian, numSingleCells, geneExpressionMatrixAllCellTypes, cellTypeNames, 0.80, "geneExpression80MatrixAllCellTypesOvarian.gz", "labels80")

# Eliminate any column (single-cell) with more
# than 85% of its genes with zero expression
eliminateColsWithMoreThanPercentGenesWithZeroExpression(allDataAllCellTypesOvarian, numSingleCells, geneExpressionMatrixAllCellTypes, cellTypeNames, 0.85, "geneExpression85MatrixAllCellTypesOvarian.gz", "labels85")

# Eliminate any column (single-cell) with more than 90%
# of its genes with zero expression
eliminateColsWithMoreThanPercentGenesWithZeroExpression(allDataAllCellTypesOvarian, numSingleCells, geneExpressionMatrixAllCellTypes, cellTypeNames, 0.90, "geneExpression90MatrixAllCellTypesOvarian.gz", "labels90")

# Eliminate any column (single-cell) with more than 98%
# of its genes with zero expression
eliminateColsWithMoreThanPercentGenesWithZeroExpression(allDataAllCellTypesOvarian, numSingleCells, geneExpressionMatrixAllCellTypes, cellTypeNames, 0.98, "geneExpression98MatrixAllCellTypesOvarian.gz", "labels98")
