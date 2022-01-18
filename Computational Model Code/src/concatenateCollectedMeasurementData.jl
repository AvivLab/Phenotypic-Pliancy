# Written by Maryl Lambros December 2021
# For concatenating all the collected measurement csv files into single large file for each measurement to store in ../dataOutputForPaper folder
#-------------------------------------------------------------------------------
# Functions:
function concatFirstLine(evolFolderName::String, dataOutputFolderName::String, measurementCsvFile::String)
    if isfile(joinpath(dataOutputFolderName, measurementCsvFile))
        # if the concatenated file already exists and trying to add very first line, then remove it to start building from scratch
        rm(joinpath(dataOutputFolderName, measurementCsvFile))
    end
    open(joinpath(evolFolderName, measurementCsvFile)) do input
        fileToWrite =  open(joinpath(dataOutputFolderName, measurementCsvFile),"a")
        write(fileToWrite, read(input))
        close(fileToWrite)
    end
end

function concatAllRestLines(evolFolderName::String, dataOutputFolderName::String, measurementCsvFile::String)
    open(joinpath(evolFolderName, measurementCsvFile)) do input
        readuntil(input, '\n')
        fileToWrite =  open(joinpath(dataOutputFolderName, measurementCsvFile),"a")
        write(fileToWrite, read(input))
        close(fileToWrite)
    end
end

#-------------------------------------------------------------------------------
# Code:
dataOutputFolderName = joinpath("..", "dataOutputForPaper")
allFolders = readdir(dataOutputFolderName)
allFolders = allFolders[occursin.("Pop_Seed",allFolders)]
PopSeeds = parse.(Int64,map(x->match(r"(?<=Pop_Seed_)(\d+)", allFolders[x]).match,1:length(allFolders)))

for i = 1:length(allFolders)
    startingPopFolderName = joinpath(dataOutputFolderName, allFolders[i])
    popSeed = PopSeeds[i]
    println("In folder: ", startingPopFolderName)
    allSecondLevelFolders = readdir(startingPopFolderName)
    evolFolders = allSecondLevelFolders[occursin.("EvolSeed",allSecondLevelFolders)]
    csvFilesToConcatenateInPopSeed = allSecondLevelFolders[occursin.("df_",allSecondLevelFolders)]
    if (length(csvFilesToConcatenateInPopSeed) > 0)
        if (i == 1)
            map(x -> concatFirstLine(startingPopFolderName, dataOutputFolderName, x), csvFilesToConcatenateInPopSeed)
        else
            map(x -> concatAllRestLines(startingPopFolderName, dataOutputFolderName, x), csvFilesToConcatenateInPopSeed)
        end
        # Loop through all the different evolutionary trajectories for the given starting population folder under
        for j = 1:length(evolFolders)
            evolFolderName = joinpath(startingPopFolderName,evolFolders[j])
            foldersInEvolDir = readdir(evolFolderName)
            csvFilesToConcatenate = foldersInEvolDir[occursin.("df_",foldersInEvolDir)]
            if (i == 1 && j == 1)
                # if the first pop seed and evol running, then keep header
                map(x -> concatFirstLine(evolFolderName, dataOutputFolderName, x), csvFilesToConcatenate)
            else
                # remove header and concatenate just the row of measurement data
                map(x -> concatAllRestLines(evolFolderName, dataOutputFolderName, x), csvFilesToConcatenate)
            end
        end
    end
end
