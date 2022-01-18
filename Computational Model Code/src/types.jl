abstract type Individuals end
#
mutable struct Individual{T} <: Individuals
    network::Matrix{T} # W matrix with GenesXGenes representing genes interactions
    initstate::Vector{Float64} # Sw1 vector: initial state vector
    EnvIndInter::Matrix{T} # W x E with GenesxEnivornments representing environment gene interactions
    EnvState::Matrix{T} # Se1 vector; environment state vectors where each column is different environment state vector, like 2nd column would be Se2 vector; environment state vector/input that is 70% different than EnvState
    develstate::Matrix{T} # output of individual state vector of Genes length, after maxConv iterations with entries that are expression level of each gene
    optstate::Matrix{T} # optimum state - stable founder's development state vector (found using generate founder function at beginning of simulation)
    stable::BitArray{} # True if Individual is stable determined if current initial state is equal to the updated initial state vector (so last iteration equal to previous iteration state vector)
    fitness::Float64 # total fitness of individual under all different environmental states
    fitnessUnderEachEnv::Vector{Float64} # fitness's of an individual under each different environmental condition separately
    pathlength::Vector{Int64} # number iterations needed to reach stability in E1
    polycombvec::Matrix{T} # theta vector if only one PRC of polycomb susceptible genes: 0 if gene *not* susceptible to polycomb and 1 if gene is susceptible to polycomb; matrix when have multiple PRCs where column is each PRC and each row is that PRC's effect on a given gene
    polycombstate::Matrix{T} # after TIMECRIT what polycomb susceptible genes are expressed or supressed: 0 if susceptible gene is suppressed (meaning polycomb response element active) in given env, and 1 if susceptible gene is expressed (meaning polycomb response element is inactive) in given env. So matrix where rows are genes are columns and envs
end

mutable struct Population{I <: Individuals}
    individuals::Vector{I}
    founder::I
end

abstract type Measures end
#
mutable struct Measure{T} <: Measures
    # Measurements to take of generated full dataset:
    ## Phenotypic pliancy measures
    geneDifferencesEnv1ToEnv2::Array{Union{T, Missing}} # Euclidean distances for: 1) IntactEnv1 compared to BrokenEnv1ToEnv2; 2) IntactEnv2 compared to BrokenEnv1ToEnv2; 3) IntactEnv1 compared to IntactEnv1toEnv2; 4) IntactEnv2 compared to IntactEnv1toEnv2; 5) BrokenEnv1ToEnv2 compared to IntactEnv1ToEnv2; and 6) IntactEnv1 compared to IntactEnv2
                                            # where rows are distances above and columns represent different combinations of PRC to remove so 3 columns where [1], [2], and [1,2] are removed
    geneDifferencesEnv2ToEnv1::Array{Union{Float64, Missing}} # Euclidean distances for: 7) IntactEnv2 compared to BrokenEnv2ToEnv1; 8) IntactEnv1 compared to BrokenEnv2ToEnv1; 9) IntactEnv2 compared to IntactEnv2toEnv1; 10) IntactEnv1 compared to IntactEnv2toEnv1; 11) BrokenEnv2ToEnv1 compared to IntactEnv2ToEnv1
                                                # where rows are distances above and columns represent different combinations of PRC to remove so 3 columns where [1], [2], and [1,2] are removed
    phenotypicPliancyScore1to2::Array{Union{Float64, Missing}} # phenotypic pliancy when switch from env 1 to env 2, where columns represent different combinations of PRC to remove so 3 columns where [1], [2], and [1,2] are removed
    phenotypicPliancyScore2to1::Array{Union{Float64, Missing}} # phenotypic pliancy when switch from env 2 to env 1, where columns represent different combinations of PRC to remove so 3 columns where [1], [2], and [1,2] are removed
    # Number of unstable individuals when remove PRC 1, 2 or 1&2 when switch environments:
    unstableDisrupted1to2::Array{Bool} # whether unstable (== true) or not when switch from env 1 to env 2, where columns represent different combinations of PRC to remove so 3 columns where [1], [2], and [1,2] are removed
    unstableDisrupted2to1::Array{Bool} # whether unstable (== true) or not when switch from env 2 to env 1, where columns represent different combinations of PRC to remove so 3 columns where [1], [2], and [1,2] are removed
    ## Network Characteristics:
    ### For popualtion:
    connectivityIngoing::Float64 # Average connectivity of ingoing connections of all genes (same in both environments)
    connectivityOutgoing::Float64 # Average connectivity of outgoing connections of all genes (same for both environments)
    btwCentrality::Float64 # Average Betweenness Centrality of all genes (same for both environments)
    ### For all genes ("effective network", i.e. remove pcg repressed genes from the network when calculating these three below):
    effectiveConnectivityIngoing::Vector{Float64} # Average effective network connectivity of ingoing connections of all genes for the population where column one represents environment 1 and column 2 represents environment 2 results
    effectiveConnectivityOutgoing::Vector{Float64} # Average effective network connectivity of outgoing connections of all genes where column one represents environment 1 and column 2 represents environment 2 results
    effectiveBtwCentrality::Vector{Float64} # Average effective network Betweenness Centrality of all genes where column one represents environment 1 and column 2 represents environment 2 results
    ### For PRC target genes:
    connectivityIngoingTargets::Float64 # Average connectivity of ingoing connections of all target PRC genes (same for both environments)
    connectivityOutgoingTargets::Float64 # Average connectivity of outgoing connections of all target PRC genes (same for both environments)
    btwCentralityTargets::Float64 # Average Betweenness Centrality of all target PRC genes (same for both environments)
    ### For target genes repressed by a PRC:
    connectivityIngoingTargetsRepressed::Vector{Float64} # Average connectivity of ingoing connections of all target PRC genes that are repressed by a PRC where column one represents environment 1 and column 2 represents environment 2 results
    connectivityOutgoingTargetsRepressed::Vector{Float64} # Average connectivity of outgoing connections of all target PRC genes that are repressed by a PRC where column one represents environment 1 and column 2 represents environment 2 results
    btwCentralityTargetsRepressed::Vector{Float64} # Average Betweenness Centrality of all target PRC genes that are repressed by a PRC where column one represents environment 1 and column 2 represents environment 2 results
    # Target Gene Characteristics:
    numPcgTargets::Float64 # number of targets of any PRC (same for both environments)
    numTargetsRepressed::Vector{Float64} # number of targets that are repressed by any PRC where column one represents environment 1 and column 2 represents environment 2 results
    ## Fitness results:
    totalFitness::Float64 # fitness of individual under both environments
    fitnessEnv1::Float64 # fitness in just env 1
    fitnessEnv2::Float64 # fitness in just env 2
    # Percent of Founder Genes Repressed that also are Repressed by Polycomb during Evolution:
    founderRepressedGenes::Vector{Float64}
end

mutable struct Measurements{I <: Measures}
    individuals::Vector{I}
    # Population Averages for all the measurements:
    phenotypicPliancyScore1to2::Array{Union{Float64, Missing}} # phenotypic pliancy when switch from env 1 to env 2, where columns represent different combinations of PRC to remove so 3 columns where [1], [2], and [1,2] are removed
                                                # row 1 is mean/average, row 2 is standard deviation of the individuals in the population
    phenotypicPliancyScore2to1::Array{Union{Float64, Missing}} # phenotypic pliancy when switch from env 2 to env 1, where columns represent different combinations of PRC to remove so 3 columns where [1], [2], and [1,2] are removed
                                                # row 1 is mean/average, row 2 is standard deviation of the individuals in the population & columns represent different combinations of PRC to remove so 3 columns where [1], [2], and [1,2] are removed
    unstableDisrupted1to2::Array{Int64}  # vector of the number of the individuals in the population that are unstable for when either PRC 1, PRC 2, or PRC 1&2 disrupted when move from env 1 to 2
    unstableDisrupted2to1::Array{Int64} # vector of the number of the individuals in the population that are unstable for when either PRC 1, PRC 2, or PRC 1&2 disrupted when move from env 2 to 1
    connectivityIngoing::Float64 # Average connectivity of ingoing connections of all genes (will be same for every individual because same W matrix, so one number for population)
    connectivityOutgoing::Float64 # Average connectivity of outgoing connections of all genes (will be same for every individual because same W matrix, so one number for population)
    btwCentrality::Float64 # Average Betweenness Centrality of all genes (will be same for every individual because same W matrix, so one number for population)
    ### Effective network characteristics (i.e. remove pcg repressed genes from network and calculate measures)
    effectiveConnectivityIngoing::Array{Float64} # Average effective network connectivity of ingoing connections of all genes for the population where row 1 is for environment 1 and row 2 is for environment 2; and column 1 is average over population and column 2 is std
    effectiveConnectivityOutgoing::Array{Float64} # Average effective network connectivity of outgoing connections of all genes where row 1 is for environment 1 and row 2 is for environment 2; and column 1 is average over population and column 2 is std
    effectiveBtwCentrality::Array{Float64} # Average effective network Betweenness Centrality of all genes where row 1 is for environment 1 and row 2 is for environment 2; and column 1 is average over population and column 2 is std
    ### For PRC target genes:
    connectivityIngoingTargets::Vector{Float64} # Average connectivity of ingoing connections of all target PRC genes where column 1 is average and column 2 is std
    connectivityOutgoingTargets::Vector{Float64} # Average connectivity of outgoing connections of all target PRC genes where column 1 is average and column 2 is std
    btwCentralityTargets::Vector{Float64} # Average Betweenness Centrality of all target PRC genes where column 1 is average and column 2 is std
    ### For target genes repressed by a PRC:
    connectivityIngoingTargetsRepressed::Array{Float64} # Average connectivity of ingoing connections of all target PRC genes that are repressed by a PRC where row 1 is for environment 1 and row 2 is for environment 2; and column 1 is average over population and column 2 is std
    connectivityOutgoingTargetsRepressed::Array{Float64} # Average connectivity of outgoing connections of all target PRC genes that are repressed by a PRC where row 1 is for environment 1 and row 2 is for environment 2; and column 1 is average over population and column 2 is std
    btwCentralityTargetsRepressed::Array{Float64} # Average Betweenness Centrality of all target PRC genes that are repressed by a PRC where row 1 is for environment 1 and row 2 is for environment 2; and column 1 is average over population and column 2 is std
    ### Target Gene Characteristics:
    numPcgTargets::Vector{Float64} # number of targets of any PRC where column 1 is average and column 2 is std
    numTargetsRepressed::Array{Float64} # number of targets that are repressed by any PRC where row 1 is for environment 1 and row 2 is for environment 2; and column 1 is average over population and column 2 is std
    ## Fitness results where first element is average of population and 2nd is standard deviation of population:
    totalFitness::Vector{Float64} # fitness of individual under both environments where column 1 is average and column 2 is std
    fitnessEnv1::Vector{Float64} # fitness in just env 1 where column 1 is average and column 2 is std
    fitnessEnv2::Vector{Float64} # fitness in just env 2 where column 1 is average and column 2 is std
    ##
    founderRepressedGenes::Array{Float64} # row 1 is mean and row 2 is std, and column 1 is env1 develstate/optstate, column 2 is env2 develstate and column 3 is env2 optstate
end

#const Time = Int64
