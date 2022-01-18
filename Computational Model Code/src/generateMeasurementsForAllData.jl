#---------------------------------------------------------------------------
# Load packagees and files needed:
using Distributions # install by Pkg.add("Distributions") in julia
using JLD #is there a JLD2 now for saving the population? Looks like JLD2 is no longer being maintained as of July 2020...JLD is still being maintained and was last updated on March 30, 2021 (as of 11/2/21)
using HDF5
using Statistics # for calculating mean/average
using LightGraphs
using SimpleWeightedGraphs # for generating graph object from network and calculating betweenness centrality
using GraphPlot
using Graphs
using Combinatorics
using Distances # for Euclidean distances between gene expression vectors

include("types.jl")
include("individuals.jl")
include("population.jl")
include("measure.jl")
include("measurements.jl")
include("pliancyFunctions.jl")

#---------------------------------------------------------------------------
# Read in .jld files for each run.
# Organization of the saved data is as follows:
    # Top level folder is named after the Starting Population random seed ("Pop_Seed_RANDSEEDINPUT"), where RANDSEEDNUM is the single input to the code to generate the data
    # Second level folders are named after the random seed setting the evolutionary trajectory (RANDSEEDVEC[j]) and also still the Starting Pop rand seed (RANDSEEDINPUT) for full name of: "Pop_Seed_RANDSEEDINPUTEvolSeed_RANDSEEDVEC[j]"
    # Third level folders are named after the Generation in evolution that taking the measurement (SAVEDATA[k]) and still starting pop rand seed and evolutionry trajectory rand seed for full name of: "Pop_Seed_RANDSEEDINPUTEvolSeed_RANDSEEDVEC[j]_Geni"

    # In these third level folders are the saved .jld files for the population
    # at that generation, which contains all the population attributes to use
    # for the measures so want to loop through all of these to take measurements
    # for 1) at different generations during evolution (so compare same generation
    # for different starting populations and trajectories), 2) after evolution
    # (so compare final 1,000th generation for different starting populations and trajectories)

    # Code below gets measurements for all the .jld files, then
    # file dataAnalysisForPaper.jl analyzes and compares the measurements

# Loop through all the top level folders for the different starting populations:
dataOutputFolderName = joinpath("..", "dataOutputForPaper")
allFolders = readdir(dataOutputFolderName)
allFolders = allFolders[occursin.("Pop_Seed",allFolders)]

#------------------------------------------------------------------------------
# Loop through all the data output folders (so different starting populations generated)
for i = 1:length(allFolders) #*******consult with Carlos how to make so can have input arguments for which folders to open so can generate measurements for multiple different folders at same time (one folder per node)
    startingPopFolderName = joinpath(dataOutputFolderName, allFolders[i])
    allSecondLevelFolders = readdir(startingPopFolderName)
    evolFolders = allSecondLevelFolders[occursin.("EvolSeed",allSecondLevelFolders)]
    # Load in constants file for generating Measurements:
    include(joinpath(startingPopFolderName,"constantsForTesting.jl"))

    # Loop through all the different evolutionary trajectories for the given starting population folder under
    for j = 1:length(evolFolders)
        evolFolderName = joinpath(startingPopFolderName,evolFolders[j])
        allThirdLevelFolders = readdir(evolFolderName)
        genFolders = allThirdLevelFolders[occursin.("Gen",allThirdLevelFolders) .& occursin.("Pop",allThirdLevelFolders)]

        # Loop through all the generations where save results:
        for k = 1:length(genFolders)
            # Read the .jld population file to generate Measure julia structure based on measurements collected in meaure.jl file
            genFolderName = joinpath(evolFolderName, genFolders[k])
            populationFileName = string(genFolders[k],".jld")
            popRead = jldopen(joinpath(genFolderName,populationFileName), "r") do file
                read(file, "Population")
            end
            # Generate Measure object/take measurements:
            popMeasurements = genPopMeasurements()
            # Collect measurements and assign results to each individual in popRead:
            map(x -> collectMeasures(popRead.individuals[x], popMeasurements.individuals[x], popRead.founder), 1:N)
            # Then calcuate the averagees and standard deviations for of the given population and store in the Measurements types
            collectPopMeasures(popMeasurements)
            # Store in the measure structure or wait to be pulled out in the
            # dataAnalysisForPaper.jl file to analyize for various comparisons/Save Measurements
            # (Population of Individuals measures) as a .jld file:
            measurementsFileName = string("Measures_",genFolders[k],".jld")
            jldopen(joinpath(genFolderName,measurementsFileName), "w") do file
	            write(file, "Measurements", popMeasurements)  # alternatively, say "@write file A"
	            close(file)
	        end
        end
    end
end
