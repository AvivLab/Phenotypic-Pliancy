# Find all non-cancer cell types in Head & Neck Cancer single-cell gene expression data
# & run step 2 of RNA-seq preprocessing and analysis (document on mlambros@mail.einstein.yu.edu
# google drive titled "RNA-sequecing analysis notes")
# AND
# analyze single-cell RNA-seq datasets (head & neck cancer and ovarian cancer datasets) using
# Seurat package and our methods for analyzing phenotypic pliancy

# Set working directory for where all files (except head & neck cancer data, unless you put them all together) are kept:
workingDir = "C:/Users/mmlam/Desktop/BergmanLabRotation/SingleCellSequencingData/IdentificationOfGradeAndOriginSpecificCellPopulationsInSerousEpithelialOvarianCancerBySingleCellRNAseq/GSE118828_RAW/"
setwd(workingDir)

# set working directory where head and neck cancer data is kept:
headAndNeckWorkingDir = "C:/Users/mmlam/Desktop/BergmanLabRotation/SingleCellSequencingData/SingleCellTranscriptomicAnalysisOfPrimaryAndMetastaticTumorEcosystemsInHeadAndNeckCancer_CorrectVersion/"


### Can either run preprocessing using Seurat package (here) OR below starting at "run step 2..." which should end in the same results as the Seurat preprocessing:
## Seurat package preprocessing:
# Need gene expression matrix only for Seurat package, so re-format data here:
library("Seurat")
library("sctransform")
library("ggplot2")
library("SCnorm")
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

#########################################################################################
## PREPROCESSING AND ANALYSES:
# Set what data doing analyses on:
ovarianCancerBool = FALSE # Set if looking at ovarian or head & neck cancer data
usePreprocessHeadAndNeckBool = FALSE # to use strict preprocessing like in Head & Neck paper or less strict like in Ovarian cancer paper
scaleBool = FALSE # If want to scale and center for each gene (so lowers variance of gene expression across all cells and also gives
                  # more weight to more noisy lowly expressed genes and decreases weight of higher expressed genes --> so will see in
                  # PlotCountDepth that density is lower for higher expressed genes when use scaling)
centerBool = FALSE # if choose to center or not; for grant analysis centerBool = TRUE; for heuristic target prediction it equals FALSE
plotNormTestBool = FALSE # If to plot and save analysis of normalization methods accuracy
normTestBool = TRUE # If run code using our normalization method vs Seurat --> TRUE if normalize by our method, and false if do Seurat method
allDataHeadAndNeckBool = TRUE # If to run with all patients' data OR all patients' normal cells and only primary and metastatic cells from patients with matching primary & metastatic samples
regressCellCycleGenesBool = FALSE
saveSeuratObjBool = TRUE
useSavedSeuratObjBool = FALSE
seuratObjFilename = 'headAndNeckLaptopWithoutCellCylceRegressOutOurNorm' 
#'ovarianOurNormalizationCellCycleRegressedOutNotStrictPreprocessingSeuratObjMarylsLaptop.rds'


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
  # Do our normalization method:
  if (normTestBool == TRUE) {
    ##### SCTransform normalized and Scaled/Centered results are stored here:
    #seuratNormalizedAndCenteredData <- as.matrix(seuratObjData[["SCT"]]@scale.data) # @data hold log-normalized of corrected UMI counts; @scale.data contains the residuals (normalized values)
    ##### Seurat normalization with NormalizeData and ScaleData functions:
    seuratNormalizationMethodData <- as.matrix(seuratObjDataNormalizationMethod[["RNA"]]@scale.data)
    ##### Normalization like in Head & neck paper (log transform and center data only)
    if (centerBool == TRUE) {
      normalizeAndCenterGeneExpMatrix <- scale(log(geneExpMatrixFromSeurat+1),center=TRUE,scale=TRUE) # center and scale by cell NOT gene here
    } else {
      normalizeAndCenterGeneExpMatrix <- scale(log(geneExpMatrixFromSeurat+1),center=FALSE,scale=FALSE)
    }
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
    saveRDS(seuratObjData, file = paste(seuratObjFilename,".rds",sep=""))
  }
}
if (useSavedSeuratObjBool == TRUE) {
# to open saved seurat object:
  seuratData <- readRDS(file = paste(seuratObjFilename,".rbs",sep=""))
}
#########################################################################################################################################################
# End data preparation
############################################################################################################################



#### Print off genes and cells after preprocessing:
cat(paste(seuratObjData@assays$RNA@counts@Dim, c("genes", "cells"), "after preprocessing"))

#### Save data as csv file:
if (normTestBool == TRUE) {
  # If our normalization method, then save that data:
  #write.csv(normalizeAndCenterGeneExpMatrix,paste(seuratObjFilename,".csv",sep=""), append=FALSE)
  write.csv(normalizeAndCenterGeneExpMatrix,paste("headAndNeckNonStrictProLogTransformOnlyNoCycleRegress",".csv",sep=""), append=FALSE)
} else {
  # if seurat normalization method, then save that data:
  write.csv(seuratNormalizationMethodData,paste(seuratObjFilename,".csv",sep=""), append=FALSE)
}


## Test amount of zeros for each category:
categoryString = "Primary"
y = rowSums(normalizeAndCenterGeneExpMatrix[1:21294,colnames(normalizeAndCenterGeneExpMatrix)==categoryString]==0)
length(y[y==0])/21294*100 # genes that have NO zero expression
length(y[y==length(colnames(normalizeAndCenterGeneExpMatrix)[colnames(normalizeAndCenterGeneExpMatrix)==categoryString])])/21294*100 # percent of genes that have ALL zero expression for all the cells in given category

# Looking at housekeeping gene expression in different categories:
x = c(1:21294)[str_detect(geneNames[[1]],"^B2M$")]
round(sum(normalizeAndCenterGeneExpMatrix[x,colnames(normalizeAndCenterGeneExpMatrix)=="Primary"]<2)/sum(colnames(normalizeAndCenterGeneExpMatrix)=="Primary")*100,2)
round(sum(normalizeAndCenterGeneExpMatrix[x,colnames(normalizeAndCenterGeneExpMatrix)=="Lymph"]<2)/sum(colnames(normalizeAndCenterGeneExpMatrix)=="Lymph")*100,2)
round(sum(normalizeAndCenterGeneExpMatrix[x,colnames(normalizeAndCenterGeneExpMatrix)=="Normal"]<2)/sum(colnames(normalizeAndCenterGeneExpMatrix)=="Normal")*100,2)
round(sum(normalizeAndCenterGeneExpMatrix[x,colnames(normalizeAndCenterGeneExpMatrix)=="Metastatic"]<2)/sum(colnames(normalizeAndCenterGeneExpMatrix)=="Metastatic")*100,2)



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
c(1:numGenes)[modesOfGenesNormal == 10] # --> use this to help pick which gene position to use in ridgeplot to visualize what having these different modes Values looks like
dev.new()
RidgePlot(seuratObjData, features = all.genes[c(2986, 3603)])#all.genes[c(1:numGenes)[modesOfGenesNormal == 10]])
#RidgePlot(seuratObjData, features = c("BTG1","RPL10","GADD45B","HSPA1A","CD81","CALM1"))#rownames(geneExpMatrixToUseForModes)[c(1:numGenes)[modesOfGenes == 7]])
# Test what genes usually mutated in cancer look like:
#dev.new()
#RidgePlot(seuratObjData, features = c("TP53", "PIK3CA", "NOTCH1", "TP63", "CDKN2A")) # --> these genes are not present in Ovarian data


#########################################################################################
# Polycomb Genes Present:
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
# Differential Expression Analysis:
#########################################################################################
# Differential Expression Analysis btw Primary & Metastatic
library("MAST") # To install: if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
# BiocManager::install("packagename")
library("DESeq2")
library("limma") #For a more efficient implementation of the Wilcoxon Rank Sum Test,
#(default method for FindMarkers) please install the limma package
#--------------------------------------------
 # install.packages('BiocManager')
#BiocManager::install('limma')
testUsing = "wilcox"
genesDifExpPrimaryAndMetast = difExpBetweenTypes(seuratObjData,geneSymbols,testUsing,"Primary","Metastatic")
pcgGenesDifExpPrimaryAndMetast = difExpPcgGenes(genesDifExpPrimaryAndMetast, pcglst, "PcG mechanism", "Primary","Metastatic")
cat(paste(length(genesDifExpPrimaryAndMetast),"genes differentially expressed between primary & metastatic"))
#dev.new()
#VlnPlot(seuratObjData, features = pcgGenesDifExpPrimaryAndMetast)
#########################################################################################
# Differential Expression Analysis btw Primary & Normal
difExpGenesPrimaryAndNormal = difExpBetweenTypes(seuratObjData,geneSymbols,testUsing,"Primary","Normal")
pcgGenesDifExpPrimaryAndNormal = difExpPcgGenes(difExpGenesPrimaryAndNormal,pcglst,"PcG mechanism","Primary","Normal")
cat(paste(length(difExpGenesPrimaryAndNormal),"genes differentially expressed between primary & normal"))
#dev.new()
#VlnPlot(seuratObjData, features = pcgGenesDifExpPrimaryAndNormal)
#########################################################################################
# Differential Expression Analysis btw Metastatic & Normal (either Lymph if head & neck OR normal if ovarian)
if (ovarianCancerBool == TRUE) {
  difExpGenesMetastAndNormal = difExpBetweenTypes(seuratObjData,geneSymbols,testUsing,"Normal","Metastatic")
  pcgGenesDifExpMetastAndNormal = difExpPcgGenes(difExpGenesMetastAndNormal,pcglst,"PcG mechanism","Normal","Metastatic")
  print(paste(length(difExpGenesMetastAndNormal),"genes differentially expressed between Normal & Metastatic"))
  #dev.new()
  #VlnPlot(seuratObjData, features = pcgGenesDifExpMetastAndNormal)
} else {
  difExpGenesMetastAndNormal = difExpBetweenTypes(seuratObjData,geneSymbols,testUsing,"Normal","Metastatic")
  pcgGenesDifExpMetastAndNormal = difExpPcgGenes(difExpGenesMetastAndNormal,pcglst,"PcG mechanism","Normal","Metastatic")
  print(paste(length(difExpGenesMetastAndNormal),"genes differentially expressed between Normal & Metastatic"))

  difExpGenesMetastAndLymph = difExpBetweenTypes(seuratObjData,geneSymbols,testUsing,"Lymph","Metastatic")
  pcgGenesDifExpMetastAndLymph = difExpPcgGenes(difExpGenesMetastAndLymph,pcglst,"PcG mechanism","Lymph","Metastatic")
  print(paste(length(difExpGenesMetastAndLymph),"genes differentially expressed between Lymph & Metastatic"))

  difExpGenesLymphAndNormal = difExpBetweenTypes(seuratObjData,geneSymbols,testUsing,"Normal","Lymph")
  pcgGenesDifExpLymphAndNormal = difExpPcgGenes(difExpGenesLymphAndNormal,pcglst,"PcG mechanism","Normal","Lymph")
  print(paste(length(difExpGenesLymphAndNormal),"genes differentially expressed between Normal & Lymph"))
}



#########################################################################################
# Differential expression using edgeR (limma) method:



#########################################################################################
# Run PCA analysis on genes differentially expressed btw primary and metastatic cells:
#########################################################################################
# remove pcg genes that are NOT targets of pcg:
pcgGenesDifExpPrimaryMetastNotUnderControl = pcgGenesDifExpPrimaryAndMetast[! pcgGenesDifExpPrimaryAndMetast %in% genesControlledByPcgPresent]
geneDifExpPrimaryMetastatNotIncludingPcg = genesDifExpPrimaryAndMetast[! genesDifExpPrimaryAndMetast %in% pcgGenesDifExpPrimaryMetastNotUnderControl]
# If want to scale and center subset of geneDifExpPrimaryMetastatNotIncludingPcg genes before running PCA:
#seuratObjData <- ScaleData(seuratObjData, features = geneDifExpPrimaryMetastatNotIncludingPcg, do.scale = TRUE, do.center = TRUE)
seuratObjData <- RunPCA(seuratObjData, features = geneDifExpPrimaryMetastatNotIncludingPcg, verbose = FALSE)
# Plot PCA:
dev.new()
DimPlot(seuratObjData, reduction = "pca")
#title(main = paste("PCA for",titleToAdd,"when use",matchingTitle,scalingTitle,"the data for genes DE btw primary & metastatic cells"))
# Plot information about PCA:
dev.new()
VizDimLoadings(seuratObjData, dims = 1:2, reduction = "pca")
dev.new()
ElbowPlot(seuratObjData)
#########################################################################################
# Remove multimodal genes from analysis:
multimodalGenes = all.genes[modesOfGenes %in% uniqueModesOfGenes[uniqueModesOfGenes > 0]]
nonMultimodalGenes = geneDifExpPrimaryMetastatNotIncludingPcg[! geneDifExpPrimaryMetastatNotIncludingPcg %in% multimodalGenes]
seuratObjData <- RunPCA(seuratObjData, features = nonMultimodalGenes, verbose = FALSE)
# Plot PCA:
dev.new()
DimPlot(seuratObjData, reduction = "pca")
#########################################################################################
# Run PCA analysis on all genes:
seuratObjData <- RunPCA(seuratObjData, features = all.genes, verbose = FALSE)
# Plot PCA:
dev.new()
DimPlot(seuratObjData, reduction = "pca")
# Look at how specific gene's expression changes in each cell on the PCA plot:
dev.new()
FeaturePlot(seuratObjData, features = pcgGenesDifExpPrimaryMetastNotUnderControl)
#########################################################################################
# Run UMAP and plot:
seuratObjData <- RunUMAP(seuratObjData, dims = 1:3)
dev.new()
DimPlot(seuratObjData, reduction = "umap")
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
} else {
  levels<-c('Lymph','Normal','Primary','Metastatic')
  pts<- c(16,1,1,1)
  colorsToUse = c("black","red","green","blue")
}
f<- factor(colnames(differences),levels)
# Plot PCA to test that same as Seurat Package PCA --> is somewhat different rotation on PC1 but
# otherwise the same. This version is actually better to use for grant!
pca<-prcomp(differences,center=TRUE,scale=TRUE)
dev.new()
plot(pca$rotation[,1],pca$rotation[,2],col=colorsToUse[f], pch=pts[f],xlab="Principal Component 1",ylab="Principal Component 2")
legend("bottomleft",NULL,levels(f),col=colorsToUse,pch=pts)

#########################################################################################
# Cluster Primary & Metastatic into different subclusters based on PCA results,
# then redo DE analysis:
if (subclusteringBool == TRUE) {
  if (ovarianCancerBool == TRUE) {
    levelsSubClustering <- c('Normal','Primary1','Primary2','Metastatic')
    colors <- c('red','cyan','green','blue')
    ptsVec <- c(1,1,1,1)
  } else {
    levelsSubClustering <- c('Lymph','Normal','Primary1','Primary2','Metastatic')
    colors <- c('black','red','cyan','green','blue')
    ptsVec <- c(16,1,1,1,1)
  }
  mylabelsSubClustering <- labelsToUse
  pca<-prcomp(differences,center=TRUE,scale=TRUE)
  pcaLoadings = seuratObjData@reductions$pca@feature.loadings
  mylabelsSubClustering[labelsToUse=='Primary' & pca$rotation[,2]>=0] <-'Primary1'
  mylabelsSubClustering[labelsToUse=='Primary' & pca$rotation[,2]<0] <-'Primary2'
  #mylabelsSubClustering[labelsToUse=='Metastatic' & pca$rotation[,1]>=(-0.02)] <-'Metastatic1'
  #mylabelsSubClustering[labelsToUse=='Metastatic' & pca$rotation[,1]<(-0.02)] <-'Metastatic2'
  factorSubClustering <- factor(mylabelsSubClustering,levelsSubClustering)
  designSubclustering <- model.matrix(~0+factorSubClustering)
  colnames(designSubclustering)=levelsSubClustering
  #colors <- c('black','red','cyan','green','blue','purple')
  dev.new()
  plot(pca$rotation[,1],pca$rotation[,2],col=colors[factorSubClustering], main="PCA with Primary Spilt", pch = ptsVec[factorSubClustering])
  legend("bottomleft",NULL,levels(factorSubClustering),col=rep(colors),pch = ptsVec)
  
  ###################################################################################################
  # Differential expression analysis after subclustering:
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
  
  difExpGenesPrimary1AndNormalSubclustering = difExpBetweenTypes(seuratDataSubcluster,geneSymbols,testUsing,"Primary1","Normal")
  pcgGenesDifExpPrimary1AndNormalSubclustering = difExpPcgGenes(difExpGenesPrimary1AndNormalSubclustering,pcglst,"PcG mechanism","Primary1","Normal")
  cat(paste(length(difExpGenesPrimary1AndNormalSubclustering),"genes differentially expressed between primary1 subcluster & normal"))
  
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
### Quantifying distance from metastatic to normal
# lymph vs primary to normal primary for just PcG
# genes after pre-processing (33 PcG-like genes):
# only use PcG genes that control others (so not including genes that are under control of PcG mechanism)
# When only 1 PcG mechanism gene differentially expressed btw primary and metastatic, then will not be able to calculate distance...
pcgdiffs<-normalizedData[pcgGenesDifExpPrimaryMetastNotUnderControl,]
if (ovarianCancerBool==TRUE) {
  pcgdatm<-pcgdiffs[,labelsToUse=='Metastatic']
  pcgdatp<-pcgdiffs[,labelsToUse=='Primary']
  pcgdatn<-pcgdiffs[,labelsToUse=='Normal']
  pcgmeansm<-rowMeans(pcgdatm)
  pcgmeansp<-rowMeans(pcgdatp)
  pcgmeansn<-rowMeans(pcgdatn)
  pcgdiff1<-pcgmeansp-pcgmeansn
  pcgdiff2<-pcgmeansm-pcgmeansn
  cat("Ovarian: Distance btw primary & normal for PcG genes:",dist(rbind(pcgmeansp,pcgmeansn)),"\nDistance btw metastatic & normal for PcG genes:",dist(rbind(pcgmeansm,pcgmeansn)),"\n")
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
  cat("Head&Neck: Distance btw primary & normal for PcG genes:",dist(rbind(pcgmeansp,pcgmeansn)),"\nDistance btw metastatic & lymph normal for PcG genes:",dist(rbind(pcgmeansm,pcgmeansl)),"\n")
}
### Quantifying distance from metastatic to normal
# lymph vs primary to normal primary for genes
# differentially expressed btw metastatic and
# primary that are NOT PcG mechanism genes or multimodal genes:
if (subclusteringBool == TRUE) {
  differences = differences[,mylabelsSubClustering!="Primary1"]
  labelsToUse = labelsToUse[mylabelsSubClustering!="Primary1"]
}
if (ovarianCancerBool==TRUE) {
  datm<-differences[,labelsToUse=='Metastatic']
  datp<-differences[,labelsToUse=='Primary']
  datn<-differences[,labelsToUse=='Normal']
  meansm<-rowMeans(datm)
  meansp<-rowMeans(datp)
  meansn<-rowMeans(datn)
  # Aim 2 H2 quantify distance results:
  cat("Ovarian: Distance btw primary & normal:",dist(rbind(meansp,meansn)),"\nDistance btw metastatic & normal:",dist(rbind(meansm,meansn)),"\n")
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
  cat("Head&Neck: Distance btw primary & normal for DE genes btw primary & metastatic:",dist(rbind(meansp,meansn)),"\nDistance btw metastatic & lymph normal for DE genes btw primary & metastatic:",dist(rbind(meansm,meansl)),"\n")
}

### Aim 2 H2 quantify gene by gene differences:
if (ovarianCancerBool==TRUE) {
  diff1<-meansp-meansn
  diff2<-meansm-meansn
  colorsGeneDiffs <- rep("black", length(diff2))
  colorsGeneDiffs[abs(diff2) < abs(diff1)] <- "red"
  dev.new()
  plot(diff1,diff2,main = "Gene by Gene Differences", xlab = "Distance btw primary & normal", ylab = "Distance btw metastatic and normal", col=colorsGeneDiffs)
  #plot(fulldiff1,fulldiff2)
  abline(0,1)
  abline(0,-1)
  abline(h=0)
  abline(v=0)
  #abline(lm(diff2~diff1))
  c(sqrt(sum(diff1^2)),sqrt(sum(diff2^2)),mean(abs(diff2)<abs(diff1)))
} else {
  diff1<-meansp-meansn
  diff2<-meansm-meansl
  colorsGeneDiffs <- rep("black", length(diff2))
  colorsGeneDiffs[abs(diff2) < abs(diff1)] <- "red"
  dev.new()
  plot(diff1,diff2,main = "Gene by Gene Differences", xlab = "Distance btw primary & normal", ylab = "Distance btw metastatic and lymph normal", col = colorsGeneDiffs)
  #plot(fulldiff1,fulldiff2)
  abline(0,1)
  abline(0,-1)
  abline(h=0)
  abline(v=0)
  #abline(lm(diff2~diff1))
  c(sqrt(sum(diff1^2)),sqrt(sum(diff2^2)),mean(abs(diff2)<abs(diff1)))
}

### Differential expression for primary1 and primary2 groups separately:

############################################################################
### Aim 2 H3: Pliancy Z-score analysis
#######################################
numPcgMechanismOnlyGenes = (length(pcgnames)-length(genesControlledByPcgPresent)) # number of pcg mechanism genes only after preprocessing, so not including genes controlled by PcG mechanism
pcgOnlyInds <- pcginds & !controlledByPcg
pcgOnlynames<-geneSymbols[pcgOnlyInds] # only PcG mechanism gene names and not genes controlled by PcG

diffinds = sapply(geneSymbols[1:length(geneSymbols)],function(x) {any(genesDifExpPrimaryAndMetast==x)})
diffnotpcginds<- diffinds & !pcgOnlyInds
mydata2<-normalizedData[geneDifExpPrimaryMetastatNotIncludingPcg,]
mydata2m<-mydata2[,labelsToUse=='Metastatic']
mydata2p<-mydata2[,labelsToUse=="Primary"]
mydata2n<-mydata2[,labelsToUse=="Normal" | labelsToUse=="Lymph"]

len<-sum(diffnotpcginds)
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
pcgdat<-normalizedData[pcgOnlyInds,labelsToUse=="Metastatic"]
pcgdat<-data.frame(cbind(t(pcgdat),mymeanscores))
colnames(pcgdat)<-c(pcgOnlynames,"score")
# Fit linear regression for overall pliancy
# score as a function of all (length(pcgnames)-3) PcG-like genes:
mymodel<-lm(score~.,data=pcgdat)
# PcG genes with positive coefficients, meaning
# decreased expression positively correlated
# with movement towards normal:
pcgGenesPos<-mymodel$coefficients[mymodel$coefficients > 0]
cat("# PcG genes w/ positive coefficients:",length(pcgGenesPos))

### Bootstrapping for Aim 2 H3 results to find pval ###
controlPos<-sample(which(!pcgOnlyInds & !diffinds), numPcgMechanismOnlyGenes)
controldat<-normalizedData[controlPos,labelsToUse=="Metastatic"]
controldat<-data.frame(cbind(t(controldat),mymeanscores))
colnames(controldat)<-c(names(controlPos),"score")
mymodel<-lm(score~.,data=controldat)

nposdist<-rep(0,1000)
for (i in 1:1000)
{
  controlPos<-sample(which(!pcgOnlyInds & !diffinds),numPcgMechanismOnlyGenes)
  controldat<-normalizedData[controlPos,labelsToUse=="Metastatic"]
  controldat<-data.frame(cbind(t(controldat),mymeanscores))
  colnames(controldat)<-c(names(controlPos),"score")
  controlmodel<-lm(score~.,data=controldat)
  npositive<-sum(controlmodel$coefficients[2:(numPcgMechanismOnlyGenes+1)]>0)
  nposdist[i]<-npositive
}
if (ovarianCancerBool == TRUE) {
  ylimMax = 20
} else {
  ylimMax = 200
}
dev.new()
hist(nposdist, breaks=13,xlim=c(0,length(pcgnames)),ylim=c(0,ylimMax))
abline(v=length(pcgGenesPos),col="red")

cat(paste("Mean and std of positive distribution:", mean(nposdist[!is.na(nposdist)]), sd(nposdist[!is.na(nposdist)])))
library("fitdistrplus")
fit <- fitdist(nposdist[!is.na(nposdist)], "norm")
pValHist <- dnorm(x= length(pcgGenesPos), mean= fit$estimate[['mean']], sd= fit$estimate[['sd']])
cat(paste("P-value of # PcG genes w/ positive coefficients:", pValHist))
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
