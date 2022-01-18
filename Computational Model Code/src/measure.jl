function Measure()
    Measure(Array{Union{Float64,Missing}}(missing, 6, 3), Array{Union{Float64,Missing}}(missing, 5, 3),
     Array{Union{Float64,Missing}}(missing, 3), Array{Union{Float64,Missing}}(missing, 3),
      zeros(Bool, 3), zeros(Bool, 3), 0., 0., 0.,
       zeros(Float64, 2), zeros(Float64, 2), zeros(Float64, 2), 0., 0., 0., zeros(Float64, 2), zeros(Float64, 2), zeros(Float64, 2),
        0., zeros(Float64, 2), 0., 0., 0., zeros(Float64, 3))
    # set Measure type with the initial values above
    # Measure(connectivityIngoing, connectivityOutgoing,
    #             btwCentrality, cliqueSize, phenotypicPliancyScore1,
    #                 phenotypicPliancyScore2, newEnvPliancyScore)
end

function genMeasures()
    measvect1 = Array{Measure{Float64}}(undef, N) # Julia assigns Measure type as whole as Float64
    for i = 1:N
        measvect1[i] = Measure()
    end
    return measvect1
end

## Network Architecture Measures:
# 1) Calculate Adjacency Matrix:
# 2) Connectivity of ingoing and outgoing connections:
# 3) Betweenness Centrality
# 4) Relative Clique Size
function convertNetworkToAdjacencyMat(indNetwork::Matrix{Float64})
    adjacencyMat = copy(indNetwork)
    adjacencyMat[adjacencyMat.!=0] .= 1
    return adjacencyMat
end

function connectivityIngoing(ind::Individual, indMeasures::Measure)
# loop through every gene and calculate its incoming connections and return
# a vector with incoming connectivity of each gene:
    incomingConnections = zeros(Float64,G)
    adjacencyMatrix = convertNetworkToAdjacencyMat(ind.network)
    for i = 1:G
        # sum over the ith row of adjacency matrix giving you the total number of connections onto node i (Wij is the effect of gene j on the product of gene i --> so a row is everything coming into gene i to affect it)
        incomingConnections[i] = sum(adjacencyMatrix[i,:])
    end
    indMeasures.connectivityIngoing = mean(incomingConnections)
    indMeasures.connectivityIngoingTargets = mean(avgForAllPrcTargets(ind, incomingConnections))
    # Connectivity in env 1 vs env 2 because target genes repressed can be different:
    indMeasures.connectivityIngoingTargetsRepressed = mean.(avgForRepressedPrcTargets(ind, incomingConnections))
end

function connectivityOutgoing(ind::Individual, indMeasures::Measure)
# loop through every gene and calculate its outgoing connections and return
# a vector with outgoing connectivity of each gene:
    outgoingConnections = zeros(Float64,G)
    adjacencyMatrix = convertNetworkToAdjacencyMat(ind.network)
    for i = 1:G
        # sum over the ith column to give you the total number of connections coming from node i to node j
        outgoingConnections[i] = sum(adjacencyMatrix[:,i])
    end
    indMeasures.connectivityOutgoing = mean(outgoingConnections)
    indMeasures.connectivityOutgoingTargets = mean(avgForAllPrcTargets(ind, outgoingConnections))
    # Connectivity in env 1 vs env 2 because target genes repressed can be different:
    indMeasures.connectivityOutgoingTargetsRepressed = mean.(avgForRepressedPrcTargets(ind, outgoingConnections))
end

function betweennessCentrality(ind::Individual, indMeasures::Measure)
# loop through every gene and calculate its incoming connections and return
# a vector with connectivity of each gene:
    adjacencyMatrix = convertNetworkToAdjacencyMat(ind.network)
    g = SimpleWeightedDiGraph(adjacencyMatrix)
    btwCentrality = Graphs.betweenness_centrality(g, normalize = true)
    indMeasures.btwCentrality = mean(btwCentrality)
    indMeasures.btwCentralityTargets = mean(avgForAllPrcTargets(ind, btwCentrality))
    # Connectivity in env 1 vs env 2 because target genes repressed can be different:
    indMeasures.btwCentralityTargetsRepressed = mean.(avgForRepressedPrcTargets(ind, btwCentrality))
end

## Effective network measures:
function effectiveConnectivityIngoing(ind::Individual, indMeasures::Measure)
# loop through every gene and calculate its incoming connections and return
# a vector with incoming connectivity of each gene:
    incomingConnections = zeros(Float64,G)
    # Remove pcg repressed genes from network (i.e. make column and row all zeros when gene is repressed by a pcg)
    # Different effective network for each environment:
    for env = 1:DIFENVS
        repressedGenes = findall(==(0.0), ind.polycombstate[:,env])
        effectiveNetwork = deepcopy(ind.network)
        effectiveNetwork[repressedGenes,:] = zeros(length(repressedGenes),G)
        effectiveNetwork[:,repressedGenes] = zeros(G,length(repressedGenes))
        adjacencyMatrix = convertNetworkToAdjacencyMat(effectiveNetwork)
        for i = 1:G
            # sum over the ith row of adjacency matrix giving you the total number of connections onto node i (Wij is the effect of gene j on the product of gene i --> so a row is everything coming into gene i to affect it)
            incomingConnections[i] = sum(adjacencyMatrix[i,:])
        end
        indMeasures.effectiveConnectivityIngoing[env] = mean(incomingConnections)
    end
end

function effectiveConnectivityOutgoing(ind::Individual, indMeasures::Measure)
# loop through every gene and calculate its incoming connections and return
# a vector with incoming connectivity of each gene:
    outgoingConnections = zeros(Float64,G)
    # Remove pcg repressed genes from network (i.e. make column and row all zeros when gene is repressed by a pcg)
    for env = 1:DIFENVS
        repressedGenes = findall(==(0.0), ind.polycombstate[:,env])
        effectiveNetwork = deepcopy(ind.network)
        effectiveNetwork[repressedGenes,:] = zeros(length(repressedGenes),G)
        effectiveNetwork[:,repressedGenes] = zeros(G,length(repressedGenes))
        adjacencyMatrix = convertNetworkToAdjacencyMat(effectiveNetwork)
        for i = 1:G
            # sum over the ith column to give you the total number of connections coming from node i to node j
            outgoingConnections[i] = sum(adjacencyMatrix[:,i])
        end
        indMeasures.effectiveConnectivityOutgoing[env] = mean(outgoingConnections)
    end
end

function effectiveBetweennessCentrality(ind::Individual, indMeasures::Measure)
# loop through every gene and calculate its incoming connections and return
# a vector with connectivity of each gene:
    # Remove pcg repressed genes from network (i.e. make column all zeros when gene is repressed by a pcg)
    for env = 1:DIFENVS
        repressedGenes = findall(==(0.0), ind.polycombstate[:,env])
        effectiveNetwork = deepcopy(ind.network)
        effectiveNetwork[repressedGenes,:] = zeros(length(repressedGenes),G)
        effectiveNetwork[:,repressedGenes] = zeros(G,length(repressedGenes))
        adjacencyMatrix = convertNetworkToAdjacencyMat(effectiveNetwork)
        g = SimpleWeightedDiGraph(adjacencyMatrix)
        btwCentrality = Graphs.betweenness_centrality(g, normalize = true)
        indMeasures.effectiveBtwCentrality[env] = mean(btwCentrality)
    end
end

## Polycomb-like mechanism targets and repressed targets:
function avgForAllPrcTargets(ind::Individual, indMeasurement::Vector{Float64})
    # Function that finds the average for a given measurment for just the target genes of any PRC
    # Genes are a target to any PRC are given in ind.polycombvec when marked at 1, where each column represents each PRC and each row represents each gene
    # Sum over the rows for genes with sum > 0 meaning they are a target gene of at least one PRC
    targetGenes = findall(>(0.0), sum(ind.polycombvec, dims=2))
    return indMeasurement[targetGenes]
end

function avgForRepressedPrcTargets(ind::Individual, indMeasurement::Vector{Float64})
    # Function that finds the average for a given measurment for just the *represseed* target genes of any PRC
    # Find cases when polycombstate equals 0 in each environemnt
    repressedGenesEnv1 = findall(==(0.0), ind.polycombstate[:,1])
    repressedGenesEnv2 = findall(==(0.0), ind.polycombstate[:,2])
    return [indMeasurement[repressedGenesEnv1], indMeasurement[repressedGenesEnv2]]
end

function numRepressedTargets(ind::Individual, indMeasures::Measure)
    # Find number of repressed polycomb target genes based on which environment in
    indMeasures.numTargetsRepressed[1] = length(findall(==(0.0), ind.polycombstate[:,1]))
    indMeasures.numTargetsRepressed[2] = length(findall(==(0.0), ind.polycombstate[:,2]))
end

function numTargetGenes(ind::Individual, indMeasures::Measure)
    # Find number of target genes for both environments where the targets are the same in both environments, which get repressed can differ
    indMeasures.numPcgTargets = length(findall(>(0.0), sum(ind.polycombvec, dims=2)))
end

function founderGenesRepressedVsPolycomb(ind::Individual, founder::Individual, indMeasures::Measure)
    # Find which genes have expression that is less than the gene threshold for polycomb for the founder's develstate, and compare these genes to the genes that end up as targets throughout evolution:
    # Env 1:
    founderGenesRepressed1 = findall(x->x < GENETHRESH, founder.develstate[:,1])
    polycombGenesRepressed1 = findall(x->x == 0, ind.polycombstate[:,1])
    genesReprssedInBoth1 = length(intersect(founderGenesRepressed1, polycombGenesRepressed1))
    # Env 2:
    founderGenesRepressed2 = findall(x->x < GENETHRESH, founder.develstate[:,2])
    polycombGenesRepressed2 = findall(x->x == 0, ind.polycombstate[:,2])
    genesReprssedInBoth2 = length(intersect(founderGenesRepressed2, polycombGenesRepressed2))
    # Env 2 optstate:
    founderGenesRepressed3 = findall(x->x < GENETHRESH, founder.optstate[:,2])
    genesReprssedInBoth3 = length(intersect(founderGenesRepressed3, polycombGenesRepressed2))

    indMeasures.founderRepressedGenes[1] = round(genesReprssedInBoth1/length(founderGenesRepressed1), digits = 2) * 100
    indMeasures.founderRepressedGenes[2] = round(genesReprssedInBoth2/length(founderGenesRepressed2), digits = 2) * 100
    indMeasures.founderRepressedGenes[3] = round(genesReprssedInBoth3/length(founderGenesRepressed3), digits = 2) * 100
    # ***Genes repressed by polycomb that are less than genethres in one env but greater than in the other environment
    # Compare for optimum state:


    # Genes repressed by polycomb in env 1 different than env 2 genes that are repressed?

    # Later: Potentially compare for that individual's development states at the given generation:


    # Later: Compare to final individuals with highest fitness:
end


function setDistances(ind::Individual, founder::Individual, indMeasures::Measure)
    # Measure euclidean distances between 11 cases for "geneDifferencesEnv1ToEnv2", where: 1) IntactEnv1 compared to BrokenEnv1ToEnv2;
    # 2) IntactEnv2 compared to BrokenEnv1ToEnv2; 3) IntactEnv1 compared to IntactEnv1toEnv2;
    # 4) IntactEnv2 compared to IntactEnv1toEnv2; 5) BrokenEnv1ToEnv2 compared to IntactEnv1ToEnv2;
    # and 6) IntactEnv1 compared to IntactEnv2
    # and for "geneDifferencesEnv2ToEnv1", where: 7) IntactEnv2 compared to BrokenEnv2ToEnv1;
    # 8) IntactEnv1 compared to BrokenEnv2ToEnv1; 9) IntactEnv2 compared to IntactEnv2toEnv1;
    # 10) IntactEnv1 compared to IntactEnv2toEnv1; and 11) BrokenEnv2ToEnv1 compared to IntactEnv2ToEnv1
    prcCombinations = collect(combinations((1:2)))
    envs = ((1,2),(2,1))
    for i = 1:DIFENVS # Collect measure when move from env 1 to env 2, and vice versus
        for j = 1:length(prcCombinations) # Collect measure when break different combinations of PRC complexes
            envMoveTo = envs[i][1]
            envMoveFrom = envs[i][2]
            disruptPcGResults = disruptPolycomb(ind, envMoveTo, envMoveFrom, founder, prcCombinations[j])
            indWithBrokenPolycombReplaceEnvs = disruptPcGResults[1]
            indWithIntactPolycombReplaceEnvs = disruptPcGResults[2]
            # Moving from env2 --> env1:
            if i == 1
                # check for stability:
                if ((indWithBrokenPolycombReplaceEnvs.stable[envMoveTo] == true) & (indWithIntactPolycombReplaceEnvs.stable[envMoveTo] == true))
                    # 1) IntactEnv2 (ind.) compared to BrokenEnv2ToEnv1
                    indMeasures.geneDifferencesEnv2ToEnv1[1,j] = euclidean(ind.develstate[:,envMoveFrom],indWithBrokenPolycombReplaceEnvs.develstate[:,envMoveTo])
                    # 2) IntactEnv1 compared to BrokenEnv2ToEnv1
                    indMeasures.geneDifferencesEnv2ToEnv1[2,j] = euclidean(ind.develstate[:,envMoveTo],indWithBrokenPolycombReplaceEnvs.develstate[:,envMoveTo])
                    # 3) IntactEnv2 compared to IntactEnv2toEnv1
                    indMeasures.geneDifferencesEnv2ToEnv1[3,j] = euclidean(ind.develstate[:,envMoveFrom],indWithIntactPolycombReplaceEnvs.develstate[:,envMoveTo])
                    # 4) IntactEnv1 compared to IntactEnv2toEnv1
                    indMeasures.geneDifferencesEnv2ToEnv1[4,j] = euclidean(ind.develstate[:,envMoveTo],indWithIntactPolycombReplaceEnvs.develstate[:,envMoveTo])
                    # 5) BrokenEnv2ToEnv1 compared to IntactEnv2ToEnv1
                    indMeasures.geneDifferencesEnv2ToEnv1[5,j] = euclidean(indWithBrokenPolycombReplaceEnvs.develstate[:,envMoveTo],indWithIntactPolycombReplaceEnvs.develstate[:,envMoveTo])
                else
                    indMeasures.geneDifferencesEnv2ToEnv1[:,j] = [missing, missing, missing, missing, missing]
                end
            # Moving from env1 --> env2:
            else
                if ((indWithBrokenPolycombReplaceEnvs.stable[envMoveTo] == true) & (indWithIntactPolycombReplaceEnvs.stable[envMoveTo] == true))
                    # 1) IntactEnv1 (ind.) compared to BrokenEnv1ToEnv2
                    indMeasures.geneDifferencesEnv1ToEnv2[1,j] = euclidean(ind.develstate[:,envMoveFrom],indWithBrokenPolycombReplaceEnvs.develstate[:,envMoveTo])
                    # 2) IntactEnv2 compared to BrokenEnv1ToEnv2
                    indMeasures.geneDifferencesEnv1ToEnv2[2,j] = euclidean(ind.develstate[:,envMoveTo],indWithBrokenPolycombReplaceEnvs.develstate[:,envMoveTo])
                    # 3) IntactEnv1 compared to IntactEnv1toEnv2
                    indMeasures.geneDifferencesEnv1ToEnv2[3,j] = euclidean(ind.develstate[:,envMoveFrom],indWithIntactPolycombReplaceEnvs.develstate[:,envMoveTo])
                    # 4) IntactEnv2 compared to IntactEnv1toEnv2
                    indMeasures.geneDifferencesEnv1ToEnv2[4,j] = euclidean(ind.develstate[:,envMoveTo],indWithIntactPolycombReplaceEnvs.develstate[:,envMoveTo])
                    # 5) BrokenEnv1ToEnv2 compared to IntactEnv1ToEnv2
                    indMeasures.geneDifferencesEnv1ToEnv2[5,j] = euclidean(indWithBrokenPolycombReplaceEnvs.develstate[:,envMoveTo],indWithIntactPolycombReplaceEnvs.develstate[:,envMoveTo])
                    # and 6) IntactEnv1 compared to IntactEnv2
                    indMeasures.geneDifferencesEnv1ToEnv2[6,j] = euclidean(ind.develstate[:,envMoveFrom],ind.develstate[:,envMoveTo])
                else
                    indMeasures.geneDifferencesEnv1ToEnv2[:,j] = [missing, missing, missing, missing, missing, missing]
                end
            end
        end
    end
end


# Set functional pliancy measure per individual:
function measureFunctionalPliancy(ind::Individual, indMeasures::Measure, founder::Individual)
    prcCombinations = collect(combinations((1:2)))
    envs = ((1,2),(2,1))
    for i = 1:DIFENVS # Collect measure when move from env 1 to env 2, and vice versus
        for j = 1:length(prcCombinations) # Collect measure when break different combinations of PRC complexes
            if i == 1
                indMeasures.phenotypicPliancyScore2to1[j] = functionalPliancyPerIndividual(ind, founder, envs[i][1], envs[i][2], prcCombinations[j])
                indMeasures.unstableDisrupted2to1[j] = ismissing(indMeasures.phenotypicPliancyScore2to1[j])
            else
                indMeasures.phenotypicPliancyScore1to2[j] = functionalPliancyPerIndividual(ind, founder, envs[i][1], envs[i][2], prcCombinations[j])
                indMeasures.unstableDisrupted1to2[j] = ismissing(indMeasures.phenotypicPliancyScore1to2[j])
            end
        end
    end
end

# Collect fitness of each individual:
function collectFitness(ind::Individual, indMeasures::Measure)
    indMeasures.totalFitness = copy(ind.fitness)
    indMeasures.fitnessEnv1 = copy(ind.fitnessUnderEachEnv[1])
    indMeasures.fitnessEnv2 = copy(ind.fitnessUnderEachEnv[2])
end


# Population averages and std measures:
function popAvgForMeasurement(popMeasurements::Measurements, measurementName::String)
    measurementVals = map(x -> getproperty(popMeasurements.individuals[x], Symbol(measurementName)), 1:N)
    avg = mean(measurementVals)
    standDev = std(measurementVals)
    return [avg'; standDev']
end

# function for population value for connectivity and betweenness centrality for network of every individual in population (network is the same so only 1 value and std is 0)
function networkMeasurementResults(popMeasurements::Measurements, measurementName::String)
    return getproperty(popMeasurements.individuals[1], Symbol(measurementName)) #all individuals will be the same
end

function pliancyPopAverages(popMeasurements::Measurements, measurementName::String)
    measurementVals = map(x -> getproperty(popMeasurements.individuals[x], Symbol(measurementName)), 1:N)
    temp = hcat(measurementVals...)
    avg = map(x -> mean(skipmissing(temp[x,:])), 1:size(temp)[1])
    standDev = map(x -> std(skipmissing(temp[x,:])), 1:size(temp)[1])
    return [avg'; standDev']
end

# Function to collect the population averages and standard deviations for the various measurments:
function collectPopMeasures(popMeasurements::Measurements)
    popMeasurements.phenotypicPliancyScore1to2 = pliancyPopAverages(popMeasurements, "phenotypicPliancyScore1to2")
    popMeasurements.phenotypicPliancyScore2to1 = pliancyPopAverages(popMeasurements, "phenotypicPliancyScore2to1")
    popMeasurements.unstableDisrupted1to2 = map(y -> length(findall(map(x -> popMeasurements.individuals[x].unstableDisrupted1to2[y], 1:N))), 1:3)
    popMeasurements.unstableDisrupted2to1 = map(y -> length(findall(map(x -> popMeasurements.individuals[x].unstableDisrupted2to1[y], 1:N))), 1:3)
    popMeasurements.numPcgTargets = popAvgForMeasurement(popMeasurements, "numPcgTargets")
    popMeasurements.numTargetsRepressed = popAvgForMeasurement(popMeasurements, "numTargetsRepressed")
    popMeasurements.totalFitness = popAvgForMeasurement(popMeasurements, "totalFitness")
    popMeasurements.fitnessEnv1 = popAvgForMeasurement(popMeasurements, "fitnessEnv1")
    popMeasurements.fitnessEnv2 = popAvgForMeasurement(popMeasurements, "fitnessEnv2")
    # Measurements for network (same for all individuals in population throughout evolution)
    popMeasurements.connectivityIngoing = networkMeasurementResults(popMeasurements, "connectivityIngoing")
    popMeasurements.connectivityOutgoing = networkMeasurementResults(popMeasurements, "connectivityOutgoing")
    popMeasurements.btwCentrality = networkMeasurementResults(popMeasurements, "btwCentrality")
    # Measurements for effective network (i.e. when take into account genes repressed by PRC(s) in the network)):
    popMeasurements.effectiveConnectivityIngoing = popAvgForMeasurement(popMeasurements, "effectiveConnectivityIngoing")
    popMeasurements.effectiveConnectivityOutgoing = popAvgForMeasurement(popMeasurements, "effectiveConnectivityOutgoing")
    popMeasurements.effectiveBtwCentrality = popAvgForMeasurement(popMeasurements, "effectiveBtwCentrality")
    # Measurements for targets and repressed genes:
    popMeasurements.connectivityIngoingTargets = popAvgForMeasurement(popMeasurements, "connectivityIngoingTargets")
    popMeasurements.connectivityOutgoingTargets = popAvgForMeasurement(popMeasurements, "connectivityOutgoingTargets")
    popMeasurements.btwCentralityTargets = popAvgForMeasurement(popMeasurements, "btwCentralityTargets")
    popMeasurements.connectivityIngoingTargetsRepressed = popAvgForMeasurement(popMeasurements, "connectivityIngoingTargetsRepressed")
    popMeasurements.connectivityOutgoingTargetsRepressed = popAvgForMeasurement(popMeasurements, "connectivityOutgoingTargetsRepressed")
    popMeasurements.btwCentralityTargetsRepressed = popAvgForMeasurement(popMeasurements, "btwCentralityTargetsRepressed")
    popMeasurements.founderRepressedGenes = popAvgForMeasurement(popMeasurements, "founderRepressedGenes")
end
## Collect all the measurements:
function collectMeasures(ind::Individual, indMeasures::Measure, founder::Individual)
    # Incoming and outgoing measurements:
    connectivityIngoing(ind, indMeasures)
    connectivityOutgoing(ind, indMeasures)
    betweennessCentrality(ind, indMeasures)
    effectiveConnectivityIngoing(ind, indMeasures)
    effectiveConnectivityOutgoing(ind, indMeasures)
    effectiveBetweennessCentrality(ind, indMeasures)
    numRepressedTargets(ind, indMeasures)
    numTargetGenes(ind, indMeasures)
    collectFitness(ind, indMeasures)
    measureFunctionalPliancy(ind, indMeasures, founder)
    setDistances(ind, founder, indMeasures)
    founderGenesRepressedVsPolycomb(ind, founder, indMeasures)
end
