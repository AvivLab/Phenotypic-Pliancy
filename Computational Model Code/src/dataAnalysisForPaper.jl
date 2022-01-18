# Created by Maryl Lambros on November 2021
# This file takes in julia structure of type "Measure" to calculate averages
# and performs analyses across different starting populations & evolutionary trajectories
#---------------------------------------------------------------------------
# Load packagees and files needed: (Note: install by Pkg.add("Distributions") in julia)
using Statistics # for calculating mean/average
using StatsPlots # for plotting
using CSV
using DataFrames

include("dataAnalysisFunctions.jl")

configfile = "constantsUsedForDataGeneration.jl"
#configfile = "constantsForTesting.jl"
indir = joinpath("..","input")
constantsFile = joinpath(indir, configfile) # constants.jl used is in the input folder in julia folder
include(constantsFile)

#-----------------------------------------------------------------------------
# What we want to output from data analysis:
# 1) General Characteristics of Population plot
    # numPcgTargets # number of targets of any PRC
    # numTargetsRepressed # number of targets that are repressed by any PRC
    # totalFitness # fitness of individual under both environments
    # fitnessEnv1 # fitness in just env 1
    # fitnessEnv2 # fitness in just env 2
# 2) Network Architecture measure plots:
    ## Each is a 2 element Vector{Float64} where 1st element is mean/average for population and 2nd is the standard deviation:
    # connectivityIngoing # Average connectivity of ingoing connections of all genes
    # connectivityOutgoing # Average connectivity of outgoing connections of all genes
    # btwCentrality # Average Betweenness Centrality of all genes
# 3) Target gene characteristics:
    ## For PRC target genes:
    # connectivityIngoingTargets # Average connectivity of ingoing connections of all target PRC genes
    # connectivityOutgoingTargets # Average connectivity of outgoing connections of all target PRC genes
    # btwCentralityTargets # Average Betweenness Centrality of all target PRC genes
    ## For target genes repressed by a PRC:
    # connectivityIngoingTargetsRepressed # Average connectivity of ingoing connections of all target PRC genes that are repressed by a PRC
    # connectivityOutgoingTargetsRepressed # Average connectivity of outgoing connections of all target PRC genes that are repressed by a PRC
    # btwCentralityTargetsRepressed # Average Betweenness Centrality of all target PRC genes that are repressed by a PRC
# 4) Phenotypic pliancy
    # phenotypicPliancyScore1to2 # phenotypic pliancy when switch from env 1 to env 2, where columns represent different combinations of PRC to remove so 3 columns where [1], [2], and [1,2] are removed
                                    # row 1 is mean/average, row 2 is standard deviation of the individuals in the population
    # phenotypicPliancyScore2to1 # phenotypic pliancy when switch from env 2 to env 1, where columns represent different combinations of PRC to remove so 3 columns where [1], [2], and [1,2] are removed
                                    # row 1 is mean/average, row 2 is standard deviation of the individuals in the population
# 5) PCA and/or UMAP of final states when PRCs left intact right after evolution completed; vs when leave PRC intact but move environment and then when break PRCs and move environment (like sequencing data analysis)


#--------------------------------------------------------------------------------
dataOutputFolderName = joinpath("..", "dataOutputForPaper", "Summary1_10000")
#dataOutputFolderName = joinpath("..", "dataOutputForPaper")

#-------------------------------------------------------------------------------
# Average for ALL populations and ALL evolutionary trajectories:
popSeed = 0 # for this section keep beccause averaging over all populations
evolSeed = 0 # for this section keep zero because averaging all evol trajectories for all populations
# When 6 different results:
plotMeasurements("df_phenotypicPliancyScore2to1.csv", dataOutputFolderName, popSeed, evolSeed)
plotMeasurements("df_phenotypicPliancyScore1to2.csv", dataOutputFolderName, popSeed, evolSeed)
plotMeasurements("df_founderRepressedGenes.csv", dataOutputFolderName, popSeed, evolSeed)
# When 4 different results:
plotMeasurements("df_btwCentralityTargetsRepressed.csv", dataOutputFolderName, popSeed, evolSeed)
plotMeasurements("df_connectivityOutgoingTargetsRepressed.csv", dataOutputFolderName, popSeed, evolSeed)
plotMeasurements("df_connectivityIngoingTargetsRepressed.csv", dataOutputFolderName, popSeed, evolSeed)
plotMeasurements("df_effectiveConnectivityIngoing.csv", dataOutputFolderName, popSeed, evolSeed)
plotMeasurements("df_effectiveConnectivityOutgoing.csv", dataOutputFolderName, popSeed, evolSeed)
plotMeasurements("df_effectiveBtwCentrality.csv", dataOutputFolderName, popSeed, evolSeed)
plotMeasurements("df_numTargetsRepressed.csv", dataOutputFolderName, popSeed, evolSeed)
# When 3 different results:
plotMeasurements("df_unstableDisrupted2to1.csv", dataOutputFolderName, popSeed, evolSeed)
plotMeasurements("df_unstableDisrupted1to2.csv", dataOutputFolderName, popSeed, evolSeed)
# When 2 different results:
plotMeasurements("df_fitnessEnv2.csv", dataOutputFolderName, popSeed, evolSeed)
plotMeasurements("df_fitnessEnv1.csv", dataOutputFolderName, popSeed, evolSeed)
plotMeasurements("df_totalFitness.csv", dataOutputFolderName, popSeed, evolSeed)
plotMeasurements("df_numPcgTargets.csv", dataOutputFolderName, popSeed, evolSeed)
plotMeasurements("df_btwCentralityTargets.csv", dataOutputFolderName, popSeed, evolSeed)
plotMeasurements("df_connectivityOutgoingTargets.csv", dataOutputFolderName, popSeed, evolSeed)
plotMeasurements("df_connectivityIngoingTargets.csv", dataOutputFolderName, popSeed, evolSeed)

# Plot euclidean distance results:
plotMultipleDfsInOne(["df_difIntact2_BrokenEnv1ToEnv2.csv",
    "df_difIntact2_IntactEnv1ToEnv2.csv", "df_difIntact1_BrokenEnv1ToEnv2.csv",
     "df_difIntact1_IntactEnv1ToEnv2.csv"], dataOutputFolderName, "PRC1and2")

# Plot together for avg for all genes, target genes, and repressed genes:
plotMultipleNetworkArchsInOne(["df_effectiveConnectivityIngoing.csv", "df_connectivityIngoingTargetsRepressed.csv", "df_connectivityIngoingTargets.csv"],
    dataOutputFolderName, "Incoming Connectivity")
plotMultipleNetworkArchsInOne(["df_effectiveConnectivityOutgoing.csv", "df_connectivityOutgoingTargetsRepressed.csv", "df_connectivityOutgoingTargets.csv"],
    dataOutputFolderName, "Outgoing Connectivity")
plotMultipleNetworkArchsInOne(["df_effectiveBtwCentrality.csv", "df_btwCentralityTargetsRepressed.csv", "df_btwCentralityTargets.csv"],
    dataOutputFolderName, "btwness Centrality")
# END ALL POPS AND EVOL TRAJS
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Average for SINGLE population but ALL evolutionary trajectories:
# User updates popSeed for single population want to plot
popSeed = 1 # User can change to input which population want to plot
evolSeed = 0 # keep zero because averaging all evol trajectories for single population
# When 6 different results:
plotMeasurements("df_phenotypicPliancyScore2to1.csv", dataOutputFolderName, popSeed, evolSeed)
plotMeasurements("df_phenotypicPliancyScore1to2.csv", dataOutputFolderName, popSeed, evolSeed)
# When 4 different results:
plotMeasurements("df_btwCentralityTargetsRepressed.csv", dataOutputFolderName, popSeed, evolSeed)
plotMeasurements("df_connectivityOutgoingTargetsRepressed.csv", dataOutputFolderName, popSeed, evolSeed)
plotMeasurements("df_connectivityIngoingTargetsRepressed.csv", dataOutputFolderName, popSeed, evolSeed)
plotMeasurements("df_effectiveConnectivityIngoing.csv", dataOutputFolderName, popSeed, evolSeed)
plotMeasurements("df_effectiveConnectivityOutgoing.csv", dataOutputFolderName, popSeed, evolSeed)
plotMeasurements("df_effectiveBtwCentrality.csv", dataOutputFolderName, popSeed, evolSeed)
plotMeasurements("df_numTargetsRepressed.csv", dataOutputFolderName, popSeed, evolSeed)
plotMeasurements("df_founderRepressedGenes.csv", dataOutputFolderName, popSeed, evolSeed)
# When 3 different results:
plotMeasurements("df_unstableDisrupted2to1.csv", dataOutputFolderName, popSeed, evolSeed)
plotMeasurements("df_unstableDisrupted1to2.csv", dataOutputFolderName, popSeed, evolSeed)
# When 2 different results:
plotMeasurements("df_fitnessEnv2.csv", dataOutputFolderName, popSeed, evolSeed)
plotMeasurements("df_fitnessEnv1.csv", dataOutputFolderName, popSeed, evolSeed)
plotMeasurements("df_totalFitness.csv", dataOutputFolderName, popSeed, evolSeed)
plotMeasurements("df_numPcgTargets.csv", dataOutputFolderName, popSeed, evolSeed)
plotMeasurements("df_btwCentralityTargets.csv", dataOutputFolderName, popSeed, evolSeed)
plotMeasurements("df_connectivityOutgoingTargets.csv", dataOutputFolderName, popSeed, evolSeed)
plotMeasurements("df_connectivityIngoingTargets.csv", dataOutputFolderName, popSeed, evolSeed)
#-------------------------------------------------------------------------------



#-------------------------------------------------------------------------------
# Average for SINGLE population and SINGLE evolutionary trajectory:
popSeed = 10 # User can change to input which population want to plot
evolSeed = 1 # User can change to input which evolutionary trajectory want to plot
# When 6 different results:
plotMeasurements("df_phenotypicPliancyScore2to1.csv", dataOutputFolderName, popSeed, evolSeed)
plotMeasurements("df_phenotypicPliancyScore1to2.csv", dataOutputFolderName, popSeed, evolSeed)
# When 4 different results:
plotMeasurements("df_btwCentralityTargetsRepressed.csv", dataOutputFolderName, popSeed, evolSeed)
plotMeasurements("df_connectivityOutgoingTargetsRepressed.csv", dataOutputFolderName, popSeed, evolSeed)
plotMeasurements("df_connectivityIngoingTargetsRepressed.csv", dataOutputFolderName, popSeed, evolSeed)
plotMeasurements("df_effectiveConnectivityIngoing.csv", dataOutputFolderName, popSeed, evolSeed)
plotMeasurements("df_effectiveConnectivityOutgoing.csv", dataOutputFolderName, popSeed, evolSeed)
plotMeasurements("df_effectiveBtwCentrality.csv", dataOutputFolderName, popSeed, evolSeed)
plotMeasurements("df_numTargetsRepressed.csv", dataOutputFolderName, popSeed, evolSeed)
plotMeasurements("df_founderRepressedGenes.csv", dataOutputFolderName, popSeed, evolSeed)
# When 3 different results:
plotMeasurements("df_unstableDisrupted2to1.csv", dataOutputFolderName, popSeed, evolSeed)
plotMeasurements("df_unstableDisrupted1to2.csv", dataOutputFolderName, popSeed, evolSeed)
# When 2 different results:
plotMeasurements("df_fitnessEnv2.csv", dataOutputFolderName, popSeed, evolSeed)
plotMeasurements("df_fitnessEnv1.csv", dataOutputFolderName, popSeed, evolSeed)
plotMeasurements("df_totalFitness.csv", dataOutputFolderName, popSeed, evolSeed)
plotMeasurements("df_numPcgTargets.csv", dataOutputFolderName, popSeed, evolSeed)
plotMeasurements("df_btwCentralityTargets.csv", dataOutputFolderName, popSeed, evolSeed)
plotMeasurements("df_connectivityOutgoingTargets.csv", dataOutputFolderName, popSeed, evolSeed)
plotMeasurements("df_connectivityIngoingTargets.csv", dataOutputFolderName, popSeed, evolSeed)
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Connectivity and betweenness centrality for each population plotted together:
# "df_btwcentrality.csv", "df_connectivityIngoing.csv", and
# "df_connectivityOutgoing.csv" contain one value for each population
# since each population's remains the same throughout evolution
plotPopLevelNetworkArchMeasures("df_btwcentrality.csv", dataOutputFolderName)
plotPopLevelNetworkArchMeasures("df_connectivityIngoing.csv", dataOutputFolderName)
plotPopLevelNetworkArchMeasures("df_connectivityOutgoing.csv", dataOutputFolderName)


#----------------------------------------------------------------------------
# 3) After evolution:
# Compare accross different starting populations and evolutionary trajectories for the 1,000th generation (final generation of evolution/after evolution):
plotFinalGenMeasurements("df_phenotypicPliancyScore2to1.csv", dataOutputFolderName, 1000)
plotFinalGenMeasurements("df_phenotypicPliancyScore1to2.csv", dataOutputFolderName, 1000)

plotFinalGenMeasurements("df_phenotypicPliancyScore2to1.csv", dataOutputFolderName, 250)
plotFinalGenMeasurements("df_phenotypicPliancyScore1to2.csv", dataOutputFolderName, 250)
