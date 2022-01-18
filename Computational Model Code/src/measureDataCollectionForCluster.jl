# Written by Maryl Lambros November 2021
# Used for collecting and storing all measurement data for Polycomb model project

#---------------------------------------------------------------------------
# Load packagees and files needed:
#install by Pkg.add("PackageName") in julia
using JLD
using HDF5
using DataFrames
using CSV
using Statistics # for taking mean

include("types.jl")
include("measureCollectionFunctions.jl")

if isempty(ARGS)
    configfile = "constantsForTesting.jl"
    evolFolderName = "..\\dataOutputForPaper\\Pop_Seed_100\\Pop_Seed_100EvolSeed_2"
else
    #configfile = "constantsUsedForDataGeneration.jl"
    configfile = "constantsForTesting.jl"
    evolFolderName = ARGS[1] # example: "..\\dataOutputForPaper\\Pop_Seed_1\\Pop_Seed_1EvolSeed_1"
end
indir = joinpath("..","input")
constantsFile = joinpath(indir, configfile) # constants.jl used is in the input folder in julia folder
include(constantsFile)

#-------------------------------------------------------------------------------
# Initialize DataFrame for each measurement:
# PopSeed represents the 1000 different starting populations that are the
# columns of the dataframe where each column is comprised of 1 - 10 EvolSeed[1-28 Generation Vector[]]

# Generation Vector containing 28 entries of a given measurement so each entry of this vector may be #, vector, or matrix

# EvolSeed vector is composed of 10 Generation Vectors (for the 10 different evolutionary trajectories ran for each starting population)
# EvolSeeds will be the rows of the dataFrame
# Loop through all the top level folders for the different starting populations and evol seeds:
popSeed = parse(Int64,match(r"(?<=Pop_Seed_)(\d+)", evolFolderName).match)
evolSeed = parse(Int64,match(r"(?<=EvolSeed_)(\d+)", evolFolderName).match)
numGens = length(SAVEDATA)
allThirdLevelFolders = readdir(evolFolderName)
genFolders = allThirdLevelFolders[occursin.("Gen",allThirdLevelFolders) .& occursin.("Pop",allThirdLevelFolders)]
gensSaveData = map(x -> SAVEDATA[x], 1:length(SAVEDATA))
if numGens != length(genFolders)
    println("Error: Not correct number of generations to match evolSeed Folder and SAVEDATA vector for evolSeed=", evolSeed, " and popSeed=", popSeed)
end
#----------------------------------------------------------------------------
# Initialize data frames for each measurement:
df_connectivityIngoing = DataFrame(Pop=popSeed, In=0.0)
df_connectivityOutgoing = DataFrame(Pop=popSeed, Out=0.0)
df_btwcentrality = DataFrame(Pop=popSeed, Btw=0.0)

# Each generation results in 2 element vector:
df_connectivityIngoingTargets = DataFrame(Pop=popSeed, Evol=evolSeed, Gen=gensSaveData, Avg=fill(0.0, numGens), Std=fill(0.0, numGens))
df_connectivityOutgoingTargets = deepcopy(df_connectivityIngoingTargets)
df_btwCentralityTargets = deepcopy(df_connectivityIngoingTargets)
df_numPcgTargets = deepcopy(df_connectivityIngoingTargets)
df_totalFitness = deepcopy(df_connectivityIngoingTargets)
df_fitnessEnv1 = deepcopy(df_connectivityIngoingTargets)
df_fitnessEnv2 = deepcopy(df_connectivityIngoingTargets)

# Each generation results in a 3 element vector:
df_unstableDisrupted1to2 = DataFrame(Pop=popSeed, Evol=evolSeed, Gen=gensSaveData, PRC1=fill(0.0, numGens), PRC2=fill(0.0, numGens), PRC1and2=fill(0.0, numGens))
df_unstableDisrupted2to1 = deepcopy(df_unstableDisrupted1to2)

# Each generation results in 2 x 2 matrix:
df_effectiveConnectivityIngoing = DataFrame(Pop=popSeed, Evol=evolSeed, Gen=gensSaveData, AvgEnv1=fill(0.0, numGens), StdEnv1=fill(0.0, numGens), AvgEnv2=fill(0.0, numGens), StdEnv2=fill(0.0, numGens))
df_effectiveConnectivityOutgoing = deepcopy(df_effectiveConnectivityIngoing)
df_effectiveBtwCentrality = deepcopy(df_effectiveConnectivityIngoing)
df_connectivityIngoingTargetsRepressed = deepcopy(df_effectiveConnectivityIngoing)
df_connectivityOutgoingTargetsRepressed = deepcopy(df_effectiveConnectivityIngoing)
df_btwCentralityTargetsRepressed = deepcopy(df_effectiveConnectivityIngoing)
df_numTargetsRepressed = deepcopy(df_effectiveConnectivityIngoing)

# Each generation results in 2 x 3 matrix:
df_phenotypicPliancyScore1to2 = DataFrame(Pop=popSeed, Evol=evolSeed, Gen=gensSaveData, AvgPRC1=fill(0.0, numGens), StdPRC1=fill(0.0, numGens), AvgPRC2=fill(0.0, numGens), StdPRC2=fill(0.0, numGens), AvgPRC1and2=fill(0.0, numGens), StdPRC1and2=fill(0.0, numGens))
df_phenotypicPliancyScore2to1 = deepcopy(df_phenotypicPliancyScore1to2)
df_founderRepressedGenes = DataFrame(Pop=popSeed, Evol=evolSeed, Gen=gensSaveData, AvgEnv1=fill(0.0, numGens), StdEnv1=fill(0.0, numGens), AvgEnv2Dev=fill(0.0, numGens), StdEnv2Dev=fill(0.0, numGens), AvgEnv2Opt=fill(0.0, numGens), StdEnv2Opt=fill(0.0, numGens))
# Split up geneDifferencesEnv1ToEnv2 and geneDifferencesEnv2ToEnv1 into dataframes of each of the 11 euclidean distances with avg and std for each for given population:
## geneDifferencesEnv1ToEnv2
df_difIntact1_BrokenEnv1ToEnv2 = DataFrame(Pop=popSeed, Evol=evolSeed, Gen=gensSaveData, AvgPRC1=fill(0.0, numGens), StdPRC1=fill(0.0, numGens), AvgPRC2=fill(0.0, numGens), StdPRC2=fill(0.0, numGens), AvgPRC1and2=fill(0.0, numGens), StdPRC1and2=fill(0.0, numGens))
df_difIntact2_BrokenEnv1ToEnv2 = deepcopy(df_difIntact1_BrokenEnv1ToEnv2)
df_difIntact1_IntactEnv1toEnv2 = deepcopy(df_difIntact1_BrokenEnv1ToEnv2)
df_difIntact2_IntactEnv1toEnv2 = deepcopy(df_difIntact1_BrokenEnv1ToEnv2)
df_difBrokenEnv1ToEnv2_IntactEnv1ToEnv2 = deepcopy(df_difIntact1_BrokenEnv1ToEnv2)
df_difIntactEnv1_IntactEnv2 = deepcopy(df_difIntact1_BrokenEnv1ToEnv2)
## geneDifferencesEnv2ToEnv1
df_difIntact1_BrokenEnv2ToEnv1 = deepcopy(df_difIntact1_BrokenEnv1ToEnv2)
df_difIntact2_BrokenEnv2ToEnv1 = deepcopy(df_difIntact1_BrokenEnv1ToEnv2)
df_difIntact1_IntactEnv2toEnv1 = deepcopy(df_difIntact1_BrokenEnv1ToEnv2)
df_difIntact2_IntactEnv2toEnv1 = deepcopy(df_difIntact1_BrokenEnv1ToEnv2)
df_difBrokenEnv2ToEnv1_IntactEnv2ToEnv1 = deepcopy(df_difIntact1_BrokenEnv1ToEnv2)

try
    #---------------------------------------------------------------------------
    # Loop through all the generations where save results:
    for k = 1:length(genFolders)
        # Read the .jld population file to generate Measure julia structure based on measurements collected in meaure.jl file
        genFolderName = joinpath(evolFolderName, genFolders[k])
        genIndex = findall(x->x==parse(Int64,match(r"(?<=Gen)(\d+)", genFolderName).match),SAVEDATA)[1]
        genFiles = readdir(genFolderName)
        measuresFile = joinpath(genFolderName,genFiles[occursin.("Measures",genFiles)][1])
        measurementDataFromFile = jldopen(measuresFile, "r") do file
            read(file, "Measurements")
        end

        # Add starting population network characteristics results to data frames for each starting population (doesn't change for evolSeed and generations)
        if (evolSeed == 1 & k == 1)
            df_connectivityIngoing[:,:In] = [round(measurementDataFromFile.connectivityIngoing, digits=3)]
            df_connectivityOutgoing[:,:Out] = [round(measurementDataFromFile.connectivityOutgoing, digits=3)]
            df_btwcentrality[:,:Btw] = [round(measurementDataFromFile.btwCentrality, digits=3)]
        end

        # each result is 2 element vector:
        setDfValueValueVector(df_connectivityIngoingTargets, measurementDataFromFile.connectivityIngoingTargets, genIndex, gensSaveData)
        setDfValueValueVector(df_connectivityOutgoingTargets, measurementDataFromFile.connectivityOutgoingTargets, genIndex, gensSaveData)
        setDfValueValueVector(df_btwCentralityTargets, measurementDataFromFile.btwCentralityTargets, genIndex, gensSaveData)
        setDfValueValueVector(df_numPcgTargets, measurementDataFromFile.numPcgTargets, genIndex, gensSaveData)
        setDfValueValueVector(df_totalFitness, measurementDataFromFile.totalFitness, genIndex, gensSaveData)
        setDfValueValueVector(df_fitnessEnv1, measurementDataFromFile.fitnessEnv1, genIndex, gensSaveData)
        setDfValueValueVector(df_fitnessEnv2, measurementDataFromFile.fitnessEnv2, genIndex, gensSaveData)

        # each result is 3 element vector:
        setDfValueValueVector(df_unstableDisrupted1to2, measurementDataFromFile.unstableDisrupted1to2, genIndex, gensSaveData)
        setDfValueValueVector(df_unstableDisrupted2to1, measurementDataFromFile.unstableDisrupted2to1, genIndex, gensSaveData)

        # each result is 2x2 matrix:
        setDfValueValueMatrix(df_effectiveConnectivityIngoing, measurementDataFromFile.effectiveConnectivityIngoing[:,1], measurementDataFromFile.effectiveConnectivityIngoing[:,2], genIndex, gensSaveData)
        setDfValueValueMatrix(df_effectiveConnectivityOutgoing, measurementDataFromFile.effectiveConnectivityOutgoing[:,1], measurementDataFromFile.effectiveConnectivityOutgoing[:,2], genIndex, gensSaveData)
        setDfValueValueMatrix(df_effectiveBtwCentrality, measurementDataFromFile.effectiveBtwCentrality[:,1], measurementDataFromFile.effectiveBtwCentrality[:,2], genIndex, gensSaveData)
        setDfValueValueMatrix(df_connectivityIngoingTargetsRepressed, measurementDataFromFile.connectivityIngoingTargetsRepressed[:,1], measurementDataFromFile.connectivityIngoingTargetsRepressed[:,2], genIndex, gensSaveData)
        setDfValueValueMatrix(df_connectivityOutgoingTargetsRepressed, measurementDataFromFile.connectivityOutgoingTargetsRepressed[:,1], measurementDataFromFile.connectivityOutgoingTargetsRepressed[:,2], genIndex, gensSaveData)
        setDfValueValueMatrix(df_btwCentralityTargetsRepressed, measurementDataFromFile.btwCentralityTargetsRepressed[:,1], measurementDataFromFile.btwCentralityTargetsRepressed[:,2], genIndex, gensSaveData)
        setDfValueValueMatrix(df_numTargetsRepressed, measurementDataFromFile.numTargetsRepressed[:,1], measurementDataFromFile.numTargetsRepressed[:,2], genIndex, gensSaveData)

        # each result is 2x3 matrix:
        setDfValueValueMatrix(df_phenotypicPliancyScore1to2, measurementDataFromFile.phenotypicPliancyScore1to2[1,:], measurementDataFromFile.phenotypicPliancyScore1to2[2,:], genIndex, gensSaveData) # row 1 is mean, row 2 is std for PRC1, PRC2, and PRC1&2
        setDfValueValueMatrix(df_phenotypicPliancyScore2to1, measurementDataFromFile.phenotypicPliancyScore2to1[1,:], measurementDataFromFile.phenotypicPliancyScore2to1[2,:], genIndex, gensSaveData) # row 1 is mean, row 2 is std for PRC1, PRC2, and PRC1&2
        setDfValueValueMatrix(df_founderRepressedGenes,  measurementDataFromFile.founderRepressedGenes[1,:], measurementDataFromFile.founderRepressedGenes[2,:], genIndex, gensSaveData)

        # geneDifferencesEnv1To2
        avgsOfGeneDistancesEnv1ToEnv2 = distancePopAvgs(measurementDataFromFile, "geneDifferencesEnv1ToEnv2")
        stdsOfGeneDistancesEnv1ToEnv2 = distancePopStd(measurementDataFromFile, "geneDifferencesEnv1ToEnv2")
        setDfValueValueMatrix(df_difIntact1_BrokenEnv1ToEnv2, avgsOfGeneDistancesEnv1ToEnv2[1,:], stdsOfGeneDistancesEnv1ToEnv2[1,:], genIndex, gensSaveData)
        setDfValueValueMatrix(df_difIntact2_BrokenEnv1ToEnv2, avgsOfGeneDistancesEnv1ToEnv2[2,:], stdsOfGeneDistancesEnv1ToEnv2[2,:], genIndex, gensSaveData)
        setDfValueValueMatrix(df_difIntact1_IntactEnv1toEnv2, avgsOfGeneDistancesEnv1ToEnv2[3,:], stdsOfGeneDistancesEnv1ToEnv2[3,:], genIndex, gensSaveData)
        setDfValueValueMatrix(df_difIntact2_IntactEnv1toEnv2, avgsOfGeneDistancesEnv1ToEnv2[4,:], stdsOfGeneDistancesEnv1ToEnv2[4,:], genIndex, gensSaveData)
        setDfValueValueMatrix(df_difBrokenEnv1ToEnv2_IntactEnv1ToEnv2, avgsOfGeneDistancesEnv1ToEnv2[5,:], stdsOfGeneDistancesEnv1ToEnv2[5,:], genIndex, gensSaveData)
        setDfValueValueMatrix(df_difIntactEnv1_IntactEnv2, avgsOfGeneDistancesEnv1ToEnv2[6,:], stdsOfGeneDistancesEnv1ToEnv2[6,:], genIndex, gensSaveData)
        # geneDifferencesEnv2ToEnv1
        avgsOfGeneDistancesEnv2ToEnv1 = distancePopAvgs(measurementDataFromFile, "geneDifferencesEnv2ToEnv1")
        stdsOfGeneDistancesEnv2ToEnv1 = distancePopStd(measurementDataFromFile, "geneDifferencesEnv2ToEnv1")
        setDfValueValueMatrix(df_difIntact2_BrokenEnv2ToEnv1, avgsOfGeneDistancesEnv2ToEnv1[1,:], stdsOfGeneDistancesEnv2ToEnv1[1,:], genIndex, gensSaveData)
        setDfValueValueMatrix(df_difIntact1_BrokenEnv2ToEnv1, avgsOfGeneDistancesEnv2ToEnv1[2,:], stdsOfGeneDistancesEnv2ToEnv1[2,:], genIndex, gensSaveData)
        setDfValueValueMatrix(df_difIntact2_IntactEnv2toEnv1, avgsOfGeneDistancesEnv2ToEnv1[3,:], stdsOfGeneDistancesEnv2ToEnv1[3,:], genIndex, gensSaveData)
        setDfValueValueMatrix(df_difIntact1_IntactEnv2toEnv1, avgsOfGeneDistancesEnv2ToEnv1[4,:], stdsOfGeneDistancesEnv2ToEnv1[4,:], genIndex, gensSaveData)
        setDfValueValueMatrix(df_difBrokenEnv2ToEnv1_IntactEnv2ToEnv1, avgsOfGeneDistancesEnv2ToEnv1[5,:], stdsOfGeneDistancesEnv2ToEnv1[5,:], genIndex, gensSaveData)
    end

    # Save data to respective CSV files:
    if (evolSeed == 1) #only save for each population seed (doesn't change for different evolution seeds and generations)
        foldersStruct = splitpath(evolFolderName)
        popSeedFolderName = joinpath(foldersStruct[1],foldersStruct[2],foldersStruct[3])
        CSV.write(joinpath(popSeedFolderName,"df_connectivityIngoing.csv"), df_connectivityIngoing)
        CSV.write(joinpath(popSeedFolderName,"df_connectivityOutgoing.csv"), df_connectivityOutgoing)
        CSV.write(joinpath(popSeedFolderName,"df_btwcentrality.csv"), df_btwcentrality)
    end

    CSV.write(joinpath(evolFolderName,"df_connectivityIngoingTargets.csv"), df_connectivityIngoingTargets)
    CSV.write(joinpath(evolFolderName,"df_connectivityOutgoingTargets.csv"), df_connectivityOutgoingTargets)
    CSV.write(joinpath(evolFolderName,"df_btwCentralityTargets.csv"), df_btwCentralityTargets)
    CSV.write(joinpath(evolFolderName,"df_numPcgTargets.csv"), df_numPcgTargets)
    CSV.write(joinpath(evolFolderName,"df_totalFitness.csv"), df_totalFitness)
    CSV.write(joinpath(evolFolderName,"df_fitnessEnv1.csv"), df_fitnessEnv1)
    CSV.write(joinpath(evolFolderName,"df_fitnessEnv2.csv"), df_fitnessEnv2)

    CSV.write(joinpath(evolFolderName,"df_unstableDisrupted1to2.csv"), df_unstableDisrupted1to2)
    CSV.write(joinpath(evolFolderName,"df_unstableDisrupted2to1.csv"), df_unstableDisrupted2to1)

    CSV.write(joinpath(evolFolderName,"df_effectiveConnectivityIngoing.csv"), df_effectiveConnectivityIngoing)
    CSV.write(joinpath(evolFolderName,"df_effectiveConnectivityOutgoing.csv"), df_effectiveConnectivityOutgoing)
    CSV.write(joinpath(evolFolderName,"df_effectiveBtwCentrality.csv"), df_effectiveBtwCentrality)
    CSV.write(joinpath(evolFolderName,"df_connectivityIngoingTargetsRepressed.csv"), df_connectivityIngoingTargetsRepressed)
    CSV.write(joinpath(evolFolderName,"df_connectivityOutgoingTargetsRepressed.csv"), df_connectivityOutgoingTargetsRepressed)
    CSV.write(joinpath(evolFolderName,"df_btwCentralityTargetsRepressed.csv"), df_btwCentralityTargetsRepressed)
    CSV.write(joinpath(evolFolderName,"df_numTargetsRepressed.csv"), df_numTargetsRepressed)

    CSV.write(joinpath(evolFolderName,"df_phenotypicPliancyScore1to2.csv"), df_phenotypicPliancyScore1to2)
    CSV.write(joinpath(evolFolderName,"df_phenotypicPliancyScore2to1.csv"), df_phenotypicPliancyScore2to1)
    CSV.write(joinpath(evolFolderName,"df_founderRepressedGenes.csv"), df_founderRepressedGenes)

    CSV.write(joinpath(evolFolderName, "df_difIntact1_BrokenEnv1ToEnv2.csv"), df_difIntact1_BrokenEnv1ToEnv2)
    CSV.write(joinpath(evolFolderName, "df_difIntact2_BrokenEnv1ToEnv2.csv"), df_difIntact2_BrokenEnv1ToEnv2)
    CSV.write(joinpath(evolFolderName, "df_difIntact1_IntactEnv1toEnv2.csv"), df_difIntact1_IntactEnv1toEnv2)
    CSV.write(joinpath(evolFolderName, "df_difIntact2_IntactEnv1toEnv2.csv"), df_difIntact2_IntactEnv1toEnv2)
    CSV.write(joinpath(evolFolderName, "df_difBrokenEnv1ToEnv2_IntactEnv1ToEnv2.csv"), df_difBrokenEnv1ToEnv2_IntactEnv1ToEnv2)
    CSV.write(joinpath(evolFolderName, "df_difIntactEnv1_IntactEnv2.csv"), df_difIntactEnv1_IntactEnv2)
    CSV.write(joinpath(evolFolderName, "df_difIntact2_BrokenEnv2ToEnv1.csv"), df_difIntact2_BrokenEnv2ToEnv1)
    CSV.write(joinpath(evolFolderName, "df_difIntact1_BrokenEnv2ToEnv1.csv"), df_difIntact1_BrokenEnv2ToEnv1)
    CSV.write(joinpath(evolFolderName, "df_difIntact2_IntactEnv2toEnv1.csv"), df_difIntact2_IntactEnv2toEnv1)
    CSV.write(joinpath(evolFolderName, "df_difIntact1_IntactEnv2toEnv1.csv"), df_difIntact1_IntactEnv2toEnv1)
    CSV.write(joinpath(evolFolderName, "df_difBrokenEnv2ToEnv1_IntactEnv2ToEnv1.csv"), df_difBrokenEnv2ToEnv1_IntactEnv2ToEnv1)
catch
    println("Error occured in ", evolFolderName)
    display(stacktrace(catch_backtrace()))
end
