# Created by Maryl Lambros 9/24/19 to import Ovarian sc-RNA-seq data
# and combine all the raw data files so can preprocess and analyze
# like Head & Neck cancer dataset.
setwd("../Desktop/BergmanLabRotation/SingleCellSequencingData/IdentificationOfGradeAndOriginSpecificCellPopulationsInSerousEpithelialOvarianCancerBySingleCellRNAseq/GSE118828_RAW/")

dataset <- read.delim("GSM3348320_V00331151_Metastatic_S1.counts.umiCounts.aboveBackground.table.csv",header = FALSE, stringsAsFactors=F, sep=",")