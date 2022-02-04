# Created by Maryl Lambros 9/24/19 to import Ovarian sc-RNA-seq data
# and combine all the raw data files so can preprocess and analyze
# like Head & Neck cancer dataset.
# Set working directory:
workingDir = "" # USER TO INPUT
setwd(workingDir)


sampleTitle <- c("PN1-P",	"BN1-P",	"HG3-M1",	"NA1-P",	"NM1",	"LG1-P",	"HG2F-P",	"HG2F-M",	"HG3-P",	"HG3-M2",	"HG4-P1",	"HG4-P2",	"HG4-P3",	"HG4-P4",	"LG2-M",	"LG2-P",	"HG1-P",	"HG1-M")
geoAccessions <- c("GSM3348303",	"GSM3348304",	"GSM3348305",	"GSM3348306",	"GSM3348307",	"GSM3348308",	"GSM3348309",	"GSM3348310",	"GSM3348311",	"GSM3348312",	"GSM3348313",	"GSM3348314",	"GSM3348315",	"GSM3348316",	"GSM3348317",	"GSM3348318",	"GSM3348319",	"GSM3348320")
sampleSite <- c("Primary",	"Primary",	"Metastatic",	"Primary",	"Normal",	"Primary",	"Primary",	"Metastatic",	"Primary",	"Metastatic",	"Primary",	"Primary",	"Primary",	"Primary",	"Metastatic",	"Primary",	"Primary",	"Metastatic")
seqFiles <- list.files(path = ".")
datasetToTestThatGeneNamesInSameOrder <- read.delim("GSM3348304_565_Cystadenoma_S1.counts.umiCounts.aboveBackground.table.csv",header = FALSE, stringsAsFactors=F, sep=",")
for (i in 1:length(geoAccessions)) {
  fileOfDatasetToAddToDf <- grep(geoAccessions[i],seqFiles)
  # make sure only 1 file with that geo accession is present; if more than one, then throw an error
  if (length(fileOfDatasetToAddToDf) == 1) {
    dataset <- read.delim(seqFiles[fileOfDatasetToAddToDf],header = FALSE, stringsAsFactors=F, sep=",")
    print(paste("number of cells in geo accession number", geoAccessions[i], "=",dim(dataset)[1]))
    } else if (length(fileOfDatasetToAddToDf) < 1) {
      print(paste("Error: no file with geo accession number", geoAccessions[i]))
      break
    } else if (length(fileOfDatasetToAddToDf) > 1) {
      print(paste("Error: more than one file with geo accession number", geoAccessions[i]))
      break
  }
  # Test that gene names are in same order as other datasets:
  if (all(dataset[1,1:dim(dataset)[2]]==datasetToTestThatGeneNamesInSameOrder[1,1:dim(dataset)[2]]) == FALSE) { # equals TRUE when genes are in same order for these two files
    print(paste("Error: Gene names are NOT in the same order for geo accession number", geoAccessions[i], "and GSM3348304"))
    break
  }
  matrixNums <- sapply(dataset[2:dim(dataset)[1],2:dim(dataset)[2]], as.numeric)
  matrixNums <- t(matrixNums)
  # Add cell type label:
  sampleSiteLabels <- rep(sampleSite[i],dim(matrixNums)[2])
  matNames <- rbind(sampleSiteLabels,matrixNums)
  # Add sample Ids for these cells:
  sampleIds <- rep(sampleTitle[i],dim(matrixNums)[2])
  matIds <- rbind(sampleIds,matNames)
  # add Gene names only for first dataset to append to data frame that will
  # contain all the data for all the files so that gene names are only in the
  # first column of this data frame containing all the data:
  if (i == 1) {
    geneNames <- as.vector(as.character(c("sample_id","gene_symbol",dataset[1,2:dim(dataset)[2]])))
    matIds <- cbind(geneNames, matIds)
    allData <- matIds
  } else {
    allData <- cbind(allData,matIds)
  }
}
allDataAllCellTypesOvarian <- as.data.frame(allData)
rownames(allDataAllCellTypesOvarian) <- 1:dim(allDataAllCellTypesOvarian)[1]
colnames(allDataAllCellTypesOvarian) <- paste0("V", 1:dim(allDataAllCellTypesOvarian)[2])
gz80 = gzfile("ovarian_all_data_all_cell_types_with_headers", "w")
write.table(allDataAllCellTypesOvarian, gz80, sep="\t", quote=F, row.names = FALSE, col.names = FALSE)
close(gz80)