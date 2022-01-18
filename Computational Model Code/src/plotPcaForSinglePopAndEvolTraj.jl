# Load packagees and files needed:
# install by Pkg.add("Distributions") in julia
using JLD
using HDF5
using Statistics # for calculating mean/average
using MultivariateStats # for PCA plotting
using UMAP: umap
using StatsPlots
using DataFrames

include("types.jl")
include("individuals.jl")
include("population.jl")
include("pliancyFunctions.jl")

include("pcaPlotFunctions.jl")

#configfile = "constantsUsedForDataGeneration.jl"
configfile = "constantsForTesting.jl"
indir = joinpath("..","input")
constantsFile = joinpath(indir, configfile) # constants.jl used is in the input folder in julia folder
include(constantsFile)

#-------------------------------------------------------------------------------
# PCA results:
# constants: c = 0.3, selstr = 1.0, and sigstr = 6.0
plotPcaForSinglePopAndEvolTraj(6, 2, 1000, 2, 1, "..\\dataOutputForPaper\\dataOutputTesting")
# constants: c = 0.1, selstr = 0.5, and sigstr = 1.0
plotPcaForSinglePopAndEvolTraj(6, 1, 1000, 2, 1, "..\\dataOutputForPaper")
plotPcaForSinglePopAndEvolTraj(100, 1, 1000, 2, 1, "..\\dataOutputForPaper")
plotPcaForSinglePopAndEvolTraj(500, 1, 1000, 2, 1, "..\\dataOutputForPaper")
plotPcaForSinglePopAndEvolTraj(950, 1, 1000, 2, 1, "..\\dataOutputForPaper")

# 10,000 Generations, population size = 1,000
plotPcaForSinglePopAndEvolTraj(6, 1, 10000, 2, 1, "..\\dataOutputForPaper\\dataOutputTenThousandGens")

# 10,000 Generations, population size = 500 (same as plotPcaForSinglePopAndEvolTraj(6, 1, 1000, 2, 1, "..\\dataOutputForPaper"))
plotPcaForSinglePopAndEvolTraj(6, 1, 1000, 2, 1, "..\\dataOutputForPaper\\dataOutputTenThousandGens\\PopulationSize500")

#-------------------------------------------------------------------------------
# UMAP results:
# UMAP in R default settings are the same as Julia default settings for UMAP
# constants: c = 0.3, selstr = 1.0, and sigstr = 6.0
plotUmapForSinglePopAndEvolTraj(6, 1, 1000, 2, 1, "..\\dataOutputForPaper\\dataOutputTesting")

# constants: c = 0.1, selstr = 0.5, and sigstr = 1.0
plotUmapForSinglePopAndEvolTraj(6, 1, 1000, 2, 1, "..\\dataOutputForPaper")
plotUmapForSinglePopAndEvolTraj(100, 2, 1000, 2, 1, "..\\dataOutputForPaper")
plotUmapForSinglePopAndEvolTraj(500, 2, 1000, 2, 1, "..\\dataOutputForPaper")
plotUmapForSinglePopAndEvolTraj(950, 2, 1000, 2, 1, "..\\dataOutputForPaper")


#-------------------------------------------------------------------------------
# Supplemental Figure of PCA when switch back to original environment after
# switching to other environment 1st (to test sustainability of switching):
suppSwitchBackPcaForSinglePopAndEvolTraj(6, 1, 1000, 2, 1, "..\\dataOutputForPaper")

# Supplemental Figure of PCA to test gene expression phenotype when stay in
# same environment after evolution so "switch" to environment already currently in
suppStayInEnvPcaForSinglePopAndEvolTraj(6, 1, 1000, 2, 1, "..\\dataOutputForPaper")
