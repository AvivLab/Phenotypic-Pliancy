setwd("C:/Users/mmlam/Desktop/BergmanLabRotation/SingleCellSequencingData/SingleCellTranscriptomicAnalysisOfPrimaryAndMetastaticTumorEcosystemsInHeadAndNeckCancer_CorrectVersion")
input_file = "fileForAnalysis90AllCellTypes.gz"#args[1] # a tab-separated table where the first column holds the gene symbol followed by as many columns as data conditions to be analyzed
labels_file = "labels85"#args[2] # labels for each column in data table (first column not important, the remaining columns will have the condition that the column holds (these will be the levels of the analysis)
fdr = 0.01#args[4] # significance level  

pcg_file="completeListPcgGenes.csv"
pcglst<-read.csv(pcg_file, header=FALSE)
head(pcglst)

# to install if (!requireNamespace("BiocManager", quietly = TRUE))
#               install.packages("BiocManager")
# BiocManager::install("limma")
# Then load the libraries after installation:
library("limma")
library("SCnorm")
library("SingleCellExperiment")
library("scran")
library("scRecover")
data=read.delim(input_file, header=F, stringsAsFactors=F, sep="\t")
labels=as.character(read.delim(labels_file, header=F, stringsAsFactors=F, sep="\t"))
numRows<-dim(data)[1]
geneSymbols<- data[2:numRows,1]
head(geneSymbols)
pcginds<-sapply(geneSymbols,function(x) {any(pcglst==x)})
pcgPos<-(1:(numRows-1))
pcgPos<-pcgPos[pcginds]
pcgnames<-geneSymbols[pcgPos]


convertLabel<- function(x) {
  if (grepl("normal",x)) {
    x<- "normal"
  } else if (grepl("lymph",x)) {
    x<- "lymph"
  }
  x
}

eucdist<-function(v1,v2) {
  sqrt(sum((v1-v2)^2))
}

setdist<-function(dat1,dat2) {
  d1<-dim(dat1)[2]
  d2<-dim(dat2)[2]
  mindists<-c()
  for (i in 1:d1){
    mydists<-c()
    for (j in 1:d2){
      mydists<-c(mydists,eucdist(dat1[,i],dat2[,j]))
    }
    mindists<-c(mindists,min(mydists))
  }
  mean(mindists)
}

plotpca<- function (myinds,mytitle,geneinterestInds,geneinterestNames,data,plotNormTestBool)
{
  levels<-c('lymph','normal','primary','metast')
  mylabels<-labels[myinds]
  mylabels<-as.character(sapply(mylabels,convertLabel))
  mydata<-sapply(data[2:numRows,myinds],as.numeric)
  rownames(mydata)<-geneSymbols
  colnames(mydata)<-mylabels
  f<- factor(mylabels,levels)
  design <- model.matrix(~0+f)
  colnames(design)=levels
  if (plotNormTestBool) {
    pdf("check_Data_count-depth_evaluationAfterNormalizationWithScranNormalizeFunction.pdf", height=5, width=7)
    plotresults <- plotCountDepth(logNormalizedWithScran, Conditions = mylabels, FilterCellProportion = .1, NCores=1) # test if need to do single-cell normalization or not (do values line up with 1 or curves not lining up on x axis? If not lining up, then need to use scnorm)
    # Since gene expression increases proportionally with sequencing depth, we expect to find the
    # estimated count-depth relationships near 1 for all genes. This is typically true for bulk RNAseq datasets. However, it does not hold in most single-cell RNA-seq datasets. In this example
    # data, the relationship is quite variable across genes.
    dev.off()
    str(plotresults)
  }
  # Normalize using median:
  ExampleSimSCData <- as.data.matrix(mydata[]
  ExampleSimSCData.CPM <- t((t(ExampleSimSCData) / colSums(ExampleSimSCData)) *
                              mean(colSums(ExampleSimSCData)))
  # Normalize data using SCnorm:
  MedExpr <- apply(mydata, 1, function(c) median(c[c != 0]))
  plot(density(log(MedExpr), na.rm=T))
  abline(v=log(c(1,2,3,4,5)))
  par(mfrow=c(2,2))
  myDataSce <- SingleCellExperiment(list(counts=mydata))
  normalizeDataWithScran <- normalize(myDataSce)
  logNormalizedWithScran = logcounts(normalizeDataWithScran)
  dataNormMetast <- SCnorm(Data = myDataSce,
                           Conditions = mylabels,
                           PrintProgressPlots = TRUE)
  
  logData=voom(mydata, plot = FALSE)
  fit=lmFit(logData,design)
  logData<-data.frame(logData)
  # Differential expression btw metastatic and primary cancer cells
  contrastTumor=makeContrasts(ccDif=metast-primary, levels=design)
  fitTumor=contrasts.fit(fit, contrastTumor)
  fit2=eBayes(fitTumor)
  results <- decideTests(fit2,method="separate",lfc=1,p.value=fdr)
  diffinds<-results[,1]!=0
  diffPositions=which(diffinds) # get the row position in the original data set of genes with differential expression
  # Data matrix containing only those genes that differentailly
  # expressed btw metastatic and primary cancer cells:
  differences=logData[diffPositions, ] # from the original data table, extract the rows for the differentially-expressed genes
  cat("Dim. of [genes x cells] of metastatic vs primary differenital expression:", dim(differences))
  
  #####################################################
  # Aim 2 H2 quantify distance from metastatic to
  # normal lymph vs primary and normal primary
  #####################################################
  ### Quantifying distance from metastatic to normal
  # lymph vs primary to normal primary for just PcG
  # genes after pre-processing (33 PcG-like genes):
  pcgdiffinds<- geneinterestInds & diffinds # only 2 PcG genes differentially expressed btw metastatic & primary though, so this analysis isn't really relevant
  pcgdiffpos<-which(pcgdiffinds)[1:2]
  pcgdiffs<-logData[pcgdiffpos,]
  pcgdatm<-pcgdiffs[,mylabels=='metast']
  pcgdatp<-pcgdiffs[,mylabels=='primary']
  pcgdatn<-pcgdiffs[,mylabels=='normal']
  pcgdatl<-pcgdiffs[,mylabels=='lymph']
  pcgmeansm<-rowMeans(pcgdatm)
  pcgmeansp<-rowMeans(pcgdatp)
  pcgmeansn<-rowMeans(pcgdatn)
  pcgmeansl<-rowMeans(pcgdatl)
  pcgdiff1<-pcgmeansp-pcgmeansn
  pcgdiff2<-pcgmeansm-pcgmeansl
  cat("Distance btw primary & normal for PcG genes:",dist(rbind(pcgmeansp,pcgmeansn)),"\nDistance btw metastatic & lymph normal for PcG genes:",dist(rbind(pcgmeansm,pcgmeansl)),"\n")
  
  ### Quantifying distance from metastatic to normal
  # lymph vs primary to normal primary for genes
  # differentially expressed btw metastatic and
  # primary (575 genes):
  datm<-differences[,mylabels=='metast']
  datp<-differences[,mylabels=='primary']
  datn<-differences[,mylabels=='normal']
  datl<-differences[,mylabels=='lymph']
  meansm<-rowMeans(datm)
  meansp<-rowMeans(datp)
  meansn<-rowMeans(datn)
  meansl<-rowMeans(datl)
  # Aim 2 H2 quantify distance results:
  cat("Distance btw primary & normal:",dist(rbind(meansp,meansn)),"\nDistance btw metastatic & lymph normal:",dist(rbind(meansm,meansl)),"\n")
  
  ### Aim 2 H2 quantify gene by gene differences:
  diff1<-meansp-meansn
  diff2<-meansm-meansl
  dev.new()
  plot(diff1,diff2,main = "Gene by Gene Differences", xlab = "Distance btw primary & normal", ylab = "Distance btw metastatic and lymph normal")
  #plot(fulldiff1,fulldiff2)
  abline(0,1)
  abline(h=0)
  abline(v=0)
  abline(lm(diff2~diff1))
  c(sqrt(sum(diff1^2)),sqrt(sum(diff2^2)),mean(abs(diff2)<abs(diff1)))
  #c(sqrt(sum(fulldiff1^2)),sqrt(sum(fulldiff2^2)),mean(abs(fulldiff1)<abs(fulldiff2)))
  #####################################################
  
  
  #####################################################
  # Aim 2 H1 
  #####################################################
  # Constrast btw primary and normal at primary site:
  contrastPN=makeContrasts(ccDif=primary-normal, levels=design)
  fitPN=contrasts.fit(fit, contrastPN)
  fitPN2=eBayes(fitPN)
  results <- decideTests(fitPN2,method="separate",lfc=1,p.value=fdr)
  diffindsPN<-results[,1]!=0
  upDiffindsPN<-results[,1]>0
  downDiffindsPN<-results[,1]<0
  cat("Pcg genes differentially expressed btw primary and normal in primary:\n", geneSymbols[diffindsPN & pcginds],
      "\nNumber of Pcg genes differentially expressed:\n", length(geneSymbols[diffindsPN & pcginds]),
      "\nPcg Genes Upregulated:\n",geneSymbols[upDiffindsPN & pcginds],
      "\nPcg Genes Downregulated:\n",geneSymbols[downDiffindsPN & pcginds],
      "\n",
      "\nStemness genes differentially expressed btw primary and normal in primary:\n", geneSymbols[diffindsPN & steminds],
      "\nNumber of stemness genes differentially expressed:\n", length(geneSymbols[diffindsPN & steminds]),
      "\nStemness Genes Upregulated:\n",geneSymbols[upDiffindsPN & steminds],
      "\nStemness Genes Downregulated:\n",geneSymbols[downDiffindsPN & steminds],"\n")
  
  # Constrast btw metastatic and normal at lymph:
  contrastMN=makeContrasts(ccDif=metast-lymph, levels=design)
  fitMN=contrasts.fit(fit, contrastMN)
  fitMN2=eBayes(fitMN)
  results <- decideTests(fitMN2,method="separate",lfc=1,p.value=fdr)
  diffindsMN<-results[,1]!=0
  upDiffindsMN<-results[,1]>0
  downDiffindsMN<-results[,1]<0
  cat("Pcg genes differentially expressed btw metastatic and normal in lymph:\n", geneSymbols[diffindsMN & pcginds],
      "\nNumber of Pcg genes differentially expressed:\n", length(geneSymbols[diffindsMN & pcginds]),
      "\nPcg Genes Upregulated:\n",geneSymbols[upDiffindsMN & pcginds],
      "\nPcg Genes Downregulated:\n",geneSymbols[downDiffindsMN & pcginds],
      "\n",
      "\nStemness genes differentially expressed btw metastatic and normal in lymph:\n", geneSymbols[diffindsMN & steminds],
      "\nNumber of stemness genes differentially expressed:\n", length(geneSymbols[diffindsMN & steminds]),
      "\nStemness Genes Upregulated:\n",geneSymbols[upDiffindsMN & steminds],
      "\nStemness Genes Downregulated:\n",geneSymbols[downDiffindsMN & steminds],"\n")
  
  
  #### Spilt primary cells into 2 subclusters based on sign of PRC 2:
  levelsSubClustering <- c('lymph','normal','primary1','primary2','metast')
  mylabelsSubClustering <- mylabels
  pca<-prcomp(differences,center=TRUE,scale.=TRUE)
  mylabelsSubClustering[mylabels=='primary' & pca$rotation[,2]>=0] <-'primary1'
  mylabelsSubClustering[mylabels=='primary' & pca$rotation[,2]<0] <-'primary2'  
  factorSubClustering <- factor(mylabelsSubClustering,levelsSubClustering)
  designSubclustering <- model.matrix(~0+factorSubClustering)
  colnames(designSubclustering)=levelsSubClustering
  fitSubclustering=lmFit(logData,designSubclustering)
  
  # Constrast between Primary 1 (positive in PRC2 -> so closer to metastatic in PRC2) and normal in primary site:
  contrastP1N=makeContrasts(ccDif=primary1-normal, levels=designSubclustering)
  fitP1N1=contrasts.fit(fitSubclustering, contrastP1N)
  fitP1N2=eBayes(fitP1N1)
  results <- decideTests(fitP1N2,method="separate",lfc=1,p.value=fdr)
  diffindsP1N<-results[,1]!=0
  upDiffindsP1N<-results[,1]>0
  downDiffindsP1N<-results[,1]<0
  cat("Pcg genes differentially expressed btw primary 1 (positive PRC2) and normal:\n", geneSymbols[diffindsP1N & pcginds],
      "\nNumber of Pcg genes differentially expressed:\n", length(geneSymbols[diffindsP1N & pcginds]),
      "\nPcg Genes Upregulated:\n",geneSymbols[upDiffindsP1N & pcginds],
      "\nPcg Genes Downregulated:\n",geneSymbols[downDiffindsP1N & pcginds],
      "\n",
      "\nStemness genes differentially expressed btw primary 1 (positive PRC2) and normal:\n", geneSymbols[diffindsP1N & steminds],
      "\nNumber of stemness genes differentially expressed:\n", length(geneSymbols[diffindsP1N & steminds]),
      "\nStemness Genes Upregulated:\n",geneSymbols[upDiffindsP1N & steminds],
      "\nStemness Genes Downregulated:\n",geneSymbols[downDiffindsP1N & steminds],"\n")
  
  # Constrast between Primary 2 (negative in PRC2 -> so farther from metastatic in PRC2) and normal in primary site:
  contrastP2N=makeContrasts(ccDif=primary2-normal, levels=designSubclustering)
  fitP2N1=contrasts.fit(fitSubclustering, contrastP2N)
  fitP2N2=eBayes(fitP2N1)
  results <- decideTests(fitP2N2,method="separate",lfc=1,p.value=fdr)
  diffindsP2N<-results[,1]!=0
  upDiffindsP2N<-results[,1]>0
  downDiffindsP2N<-results[,1]<0
  cat("Pcg genes differentially expressed btw primary 2 (negative PRC2) and normal:\n", geneSymbols[diffindsP2N & pcginds],
      "\nNumber of Pcg genes differentially expressed:\n", length(geneSymbols[diffindsP2N & pcginds]),
      "\nPcg Genes Upregulated:\n",geneSymbols[upDiffindsP2N & pcginds],
      "\nPcg Genes Downregulated:\n",geneSymbols[downDiffindsP2N & pcginds],
      "\n",
      "\nStemness genes differentially expressed btw primary 2 (negative PRC2) and normal:\n", geneSymbols[diffindsP2N & steminds],
      "\nNumber of stemness genes differentially expressed:\n", length(geneSymbols[diffindsP2N & steminds]),
      "\nStemness Genes Upregulated:\n",geneSymbols[upDiffindsP2N & steminds],
      "\nStemness Genes Downregulated:\n",geneSymbols[downDiffindsP2N & steminds],"\n")
  #####################################################
  
  
  fulldatm<-logData[,mylabels=='metast']
  fulldatp<-logData[,mylabels=='primary']
  fulldatn<-logData[,mylabels=='normal']
  fulldatl<-logData[,mylabels=='lymph']
  #fullmeansm<-rowMeans(fulldatm)
  #fullmeansp<-rowMeans(fulldatp)
  #fullmeansn<-rowMeans(fulldatn)
  #fullmeansl<-rowMeans(fulldatl)
  
  #fulldiff1<-fullmeansp-fullmeansn
  #fulldiff2<-fullmeansm-fullmeansl
  
  #####################################################
  ### Aim 2 H2 visualization ###
  #####################################################
  ### Visualization with PCA ###
  pca<-prcomp(differences,center=TRUE,scale.=TRUE)
  pca2<-pca$rotation[,1:2]
  mykmeans<-kmeans(pca2,2)
  #png(paste(mytitle,".png",sep=""))
  pts<- c(16,1,1,1)
  dev.new()
  plot(pca$rotation[,1],pca$rotation[,2],col=f, pch=pts[f],xlab="Principal Component 1",ylab="Principal Component 2")
  legend("bottomleft",NULL,levels(f),col=1:length(f),pch=pts)
  
  # Split color primary cells via sign in PRC2 and plot:
  colors <- c('black','red','cyan','green','blue')
  #png(paste(mytitle,"SubClustering.png",sep=""))
  dev.new()
  plot(pca$rotation[,1],pca$rotation[,2],col=colors[factorSubClustering], main=mytitle, pch = c(16,1,1,1,1)[factorSubClustering])
  legend("bottomleft",NULL,levels(factorSubClustering),col=rep(colors),pch = c(16,1,1,1,1))
  ###
  
  
  ### Visualization with UMAP ###
  library("umap")
  umapResults <- umap(differences)
  pts<- c(16,1,1,1)
  dev.new()
  plot(umapResults,col=f, pch=pts[f],xlab="UMAP 1",ylab="UMAP 2")
  legend("bottomleft",NULL,levels(f),col=1:length(f),pch=pts)
  
  ### Aim 2 H2 visualization with Mapper ###
  top10pca<-pca$rotation[,1:5]
  library(TDAmapper)
  mapper1<-mapper(dist_object = dist(top10pca),
                  filter_values=pca$rotation[,1],
                  num_intervals=10,
                  percent_overlap = 50,
                  num_bins_when_clustering = 6)
  mapper2 <- mapper2D(distance_matrix=dist(top10pca),
                      filter_values=list(pca$rotation[,1],pca$rotation[,2]),
                      num_intervals=c(8,8),
                      percent_overlap=40,
                      num_bins_when_clustering=5)
  library(igraph)
  g1 <- graph.adjacency(mapper1$adjacency, mode="undirected")
  g2 <- graph.adjacency(mapper2$adjacency, mode="undirected")
  for (i in 1:mapper1$num_vertices) {
    curpts<-mapper1$points_in_vertex[[i]]
    numm<-sum(mylabels[curpts]=='metast')
    nump<-sum(mylabels[curpts]=='primary')
    numn<-sum(mylabels[curpts]=='normal')
    numl<-sum(mylabels[curpts]=='lymph')
    clrs<-c("blue","green","red","black")
    j<-which.max(c(numm,nump,numn,numl))
    V(g1)$color[i]<-clrs[j]
  }
  for (i in 1:mapper2$num_vertices) {
    curpts<-mapper2$points_in_vertex[[i]]
    numm<-sum(mylabels[curpts]=='metast')
    nump<-sum(mylabels[curpts]=='primary')
    numn<-sum(mylabels[curpts]=='normal')
    numl<-sum(mylabels[curpts]=='lymph')
    clrs<-c("blue","green","red","black")
    j<-which.max(c(numm,nump,numn,numl))
    V(g2)$color[i]<-clrs[j]
  }
  #g2 <- graph.adjacency(myMapper$adjacency, mode="undirected")
  dev.new()
  plot(g1)
  legend("right",c("Normal","Primary","Metastatic"),fill=c("red","green","blue"),cex=0.8)
  
  dev.new()
  plot(g2)
  #####################################################
  
  
  #####################################################
  # Aim 2 H3
  #####################################################
  controlledByPcg <- sapply(geneSymbols,function(x) {any(c("CDH13","CDH23","RAP1B")==x)})
  pcgOnlyInds <- geneinterestInds & !controlledByPcg
  pcgOnlyPos<-(1:(numRows-1))
  pcgOnlyPos<-pcgOnlyPos[pcgOnlyInds]
  pcgOnlynames<-geneSymbols[pcgOnlyPos]
  
  diffnotpcginds<-diffinds & !pcgOnlyInds
  mydata2<-logData[diffnotpcginds,]
  mydata2m<-mydata2[mylabels=="metast"]
  mydata2p<-mydata2[mylabels=="primary"]
  mydata2n<-mydata2[mylabels=="normal" | mylabels=="lymph"]
  
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
  mymeanscores<-colMeans(sapply(myscores,as.numeric))
  pcgdat<-logData[pcgOnlyInds,mylabels=="metast"]
  pcgdat<-data.frame(cbind(t(pcgdat),mymeanscores))
  colnames(pcgdat)<-c(pcgOnlynames,"score")
  # Fit linear regression for overall pliancy
  # score as a function of all 30 PcG-like genes:
  mymodel<-lm(score~.,data=pcgdat)
  # PcG genes with positive coefficients, meaning
  # decreased expression positively correlated
  # with movement towards normal:
  pcgGenesPos<-mymodel$coefficients[mymodel$coefficients > 0]
  cat("# PcG genes w/ positive coefficients:",length(pcgGenesPos))
  
  ### Bootstrapping for Aim 2 H3 results to find pval ###
  controlPos<-sample(which(!pcgOnlyInds & !diffinds),30)
  controldat<-logData[controlPos,mylabels=="metast"]
  controldat<-data.frame(cbind(t(controldat),mymeanscores))
  colnames(controldat)<-c(names(controlPos),"score")
  mymodel<-lm(score~.,data=controldat)
  
  nposdist<-rep(0,1000)
  for (i in 1:1000)
  {
    controlPos<-sample(which(!pcgOnlyInds & !diffinds),30)
    controldat<-logData[controlPos,mylabels=="metast"]
    controldat<-data.frame(cbind(t(controldat),mymeanscores))
    colnames(controldat)<-c(names(controlPos),"score")
    controlmodel<-lm(score~.,data=controldat)
    npositive<-sum(controlmodel$coefficients[2:31]>0)
    nposdist[i]<-npositive
  }
  dev.new()
  hist(nposdist, breaks=13,xlim=c(14,28),ylim=c(0,200))
  abline(v=27,col="red")
  ########################################################
  
}

##endo separately
indsendo<-labels=="metast" | labels=="primary" | labels=="normalEndo"| labels=="lymphEndo"
plotpca(indsendo,"Endothelial",pcginds,pcgnames)

##fibro separately
indsfibro<-labels=="metast" | labels=="primary" | labels=="normalFibro" | labels=="lymphFibro"
plotpca(indsfibro,"Fibroblast",pcginds,pcgnames)

##endo+fibro
indsendofibro<- indsendo | indsfibro
plotpca(indsendofibro,"EndoFibro",pcginds,pcgnames)

##all combined
indscombo<- labels !="gene_symbol" & labels != "normalMyo" & labels!="normal-Fibro" & labels!="normalMast" & labels!= "normalB"
plotpca(indscombo,"Combined",pcginds,pcgnames)
# to run all combined by running blocks of code for each Hypothesis (H_) in plotpca function (so not running whole function)
myinds = indscombo
mytitle="Combined"
geneinterestInds = pcginds
geneinterestNames = pcgnames