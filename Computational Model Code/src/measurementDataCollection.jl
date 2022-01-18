# Written by Maryl Lambros December 2021
# Used for collecting and storing all measurement data for Polycomb model project when in Atom locally for testing

dataOutputFolderName = joinpath("..", "dataOutputForPaper")
allFolders = readdir(dataOutputFolderName)
allFolders = allFolders[occursin.("Pop_Seed",allFolders)]
PopSeeds = parse.(Int64,map(x->match(r"(?<=Pop_Seed_)(\d+)", allFolders[x]).match,1:length(allFolders)))

#---------------------------------------------------------------------------
# Loop through all the top level folders for the different starting populations:
for i = 1:length(allFolders)
    startingPopFolderName = joinpath(dataOutputFolderName, allFolders[i])
    popSeed = PopSeeds[i]
    allSecondLevelFolders = readdir(startingPopFolderName)
    evolFolders = allSecondLevelFolders[occursin.("EvolSeed",allSecondLevelFolders)]
    # Loop through all the different evolutionary trajectories for the given starting population folder under
    for j = 1:length(evolFolders)
        evolFolderName = joinpath(startingPopFolderName,evolFolders[j])
        empty!(ARGS)
        push!(ARGS,evolFolderName)
        include("measureDataCollectionForCluster.jl")
    end
end
