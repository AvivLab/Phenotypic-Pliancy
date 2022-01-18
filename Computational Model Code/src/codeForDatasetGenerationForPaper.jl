# Code to generate the full dataset used in ___(title)______ polycomb paper
# Written by Maryl Lambros from 2019 - 2021

# Load Packages: ********Check to make sure which packages we actually need for this code************
using Distributions # install by Pkg.add("Distributions") in julia
using JLD #is there a JLD2 now for saving the population? Looks like JLD2 is no longer being maintained as of July 2020...JLD is still being maintained and was last updated on March 30, 2021 (as of 11/2/21)
using HDF5
using Combinatorics
using Distributed
using Random
using DelimitedFiles # for writedlm function
#------------------------------------------------------------------------------
# Setup
include("utilities.jl")
include("types.jl")
include("individuals.jl")
include("population.jl")

#------------------------------------------------------------------------------------------------------
# Code to generate and save full simulated dataset:
# Generate folder to save datasets in:
if isempty(ARGS) # If running script directly from Atom
	configfile = "constantsForTesting.jl"
	indir = joinpath("..","input")
	constantsFile = joinpath(indir, configfile) # constants.jl used is in the input folder in julia folder
	include(constantsFile)
	RANDOMSEEDINPUT = RANDSEEDNUM
else # If running script on cluster and inputing randomSeed Argument through that
	configfile = "constantsUsedForDataGeneration.jl"
	indir = joinpath("..","input")
	constantsFile = joinpath(indir, configfile) # constants.jl used is in the input folder in julia folder
	include(constantsFile)
	RANDOMSEEDINPUT = parse(Int64,ARGS[1]) # Argument that represents the random seed
end

# to use for the generating the starting population and also naming the folder
topFolderName = string("Pop_Seed_",RANDOMSEEDINPUT)
dataPath = joinpath("..", "dataOutputForPaper")
#pathTesting = joinpath(dataPath,"dataOutputTesting")
pathTesting = joinpath(dataPath, "dataOutputTenThousandGens\\PopulationSize500")
dataDir = joinpath(pathTesting, topFolderName)
#dataDir = joinpath(dataPath, topFolderName)
# save constants.jl file used for run of this code
configpath = joinpath(dataDir, configfile)
if !isdir(dataDir)
    mkdir(dataDir)
end
# save constants file in this folder
cp(constantsFile, configpath)
# Return directory/folder name:
# print("\n",dataDir)
# Save the random seed vectors
file1 = open(joinpath(dataDir, "RandomSeedVector.csv"), "a");
writedlm(file1, RANDSEEDVEC, ',')
close(file1);

# if start with the same starting population for every independent trial but have different evolution (set in line Random.seed!(RANDSEEDVEC[j+1]);)
if SAMESTARTINGPOPFORINDTRIALSFLAG == true
    Random.seed!(RANDOMSEEDINPUT);
    pop = genpop();
    #println("\nRANDSEEDVEC = ", RANDSEEDVEC[1])
end


#------------------------------------------------------------------------------
# Generate dataset:
# Start looping through all independent trials to undergo evolution and track results for analyses:
# *****To run for different number of PRCs during evolution, import in the same RANDSEEDVEC vector and then evolve with 1, 2, or 3 PRCs. Need to add this to the code
for j = 1:INDTRIALS
	#print("\nevolution for pop. number: ", j, " gens = ", GENS)
	# Set random seed for evolution:
	Random.seed!(RANDSEEDVEC[j]);

	## Generate population and save results:
	# Generate new population with new starting w matrix, etc and different random seed for evolution, then create folder to save results under:
	if SAMESTARTINGPOPFORINDTRIALSFLAG == false
	    pop1 = genpop();
	    # create folder for this trial/separate population:
		dataDirTrial = dataDir

	# Copy same population so has same starting point, but different random seed so undergoes different evolutionary trajectory, then create folders to save results:
	else
	    pop1 = deepcopy(pop);
	    # create folder for this starting population and then save folders under this with the same starting populations but different evolutions (random seeds)
	    # Save under new folder for this trial when have the same starting population but different evolution under the starting population file
	    sameStartPopRandSeedForEvol = string(topFolderName,"EvolSeed_",RANDSEEDVEC[j])
	    dataDirTrial = joinpath(dataDir, sameStartPopRandSeedForEvol)
	    if !isdir(dataDirTrial)
	        mkdir(dataDirTrial)
	    end
	end

	## Population undergoes evolution:
	k = 1; # counter for saving population data
	for i = 0:GENS
	    ## Save population data at the specified generations in the constants file:
	    if i == SAVEDATA[k]
	        ## Save population data during evolution as julia object file (.jld)
	        # Filename for population saved ("Trial#_Gen#.csv")
	        populationFolderName = joinpath(dataDirTrial, string(sameStartPopRandSeedForEvol, "_Gen", i))
			if !isdir(populationFolderName)
	            mkdir(populationFolderName)
	        end
	        populationFilePath = joinpath(populationFolderName, string(sameStartPopRandSeedForEvol, "_Gen", i, ".jld"))
			# Save Population data for given generation in SAVEDATA:
	        jldopen(populationFilePath, "w") do file
	            write(file, "Population", pop1)  # alternatively, say "@write file A"
	            close(file)
	        end

			# Increase counter for saving population data:
	        k = k + 1;
	    end

		# Population undergoes one generation of evolution:
	    update(pop1);
	end
end
