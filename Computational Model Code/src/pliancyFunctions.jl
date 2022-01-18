# Phenotypic Pliancy Functions written by Maryl Lambros on November 2021

## Functional Phenotypic pliancy measure:
# 1) Break polycomb-like mechanism
# 2) Remove different PRC's at a time up until remove both 2
# 3) Alter environment by switching between environments 1 and 2, and vice versus
# 4) Compare phenotype to final evolved phenotype in environment moved to to phenotype when polycomb broken and placed in environment moved to to make sure when broken it is moving closer to environment moved to
# 5) Compare phenotype in environment moved to when polycomb broken and left intact to make sure broken moves closer than intact to new environment
# 6) Calculate phenotypic pliancy score (number of genes that pass criteria 4 & 5)

function breakPRC(me::Individual, prcToBreak::Array{Int64}, envState::Int64)
    # prcToBreak can be one or multiple Array element(s), where just one many PRC to break
    # (prcToBreak = [1]), OR can give multiple PRC's to break (prcToBreak = [1,3])
    polycombstateVec = copy(me.polycombstate[:,envState])
    polycombpos = findall(x->x == 0, polycombstateVec) # find polycomb response element on genes that are active (gene is supressed by polycomb), meaning polycombstate = 0
    polycombVec = copy(me.polycombvec)
    if length(polycombpos) > 0 # check that at least one polycomb response element is active before breaking one
        prcToBreakMat = polycombVec[:,prcToBreak] # pull out only polycombvec for prcToBreak
        prcToBreakPos = findall(prcToBreakMat.>0) # find all genes (positions) that are susceptible to give prcToBreak (polycombVec == 1 if susceptible)
        if length(prcToBreakPos) > 0 # If there is at least one gene susceptible to prcToBreak
            prcToBreakMat[prcToBreakPos] .= 0 # set these susceptible genes to no longer be susceptible, so set polycombvec to 0
            polycombVec[:,prcToBreak] = prcToBreakMat # update polycombvec for prcToBreak so that all susceptible genes are set to 0 so that no longer susceptible
            prcStillPresentVec = sum(polycombVec, dims = 2) # sum rows of polycombvec to find genes that are still susceptible to other active PRCs (prcStillPresentVec == 1 then gene susceptible to another active PRC)
            prcNotPresent = findall(x -> x == 0, prcStillPresentVec) # find position of genes that are NOT still susceptible to other active PRCs
            presToBreak = map(x -> prcNotPresent[x][1], 1:length(prcNotPresent)) # turn prcNotPresent from Cartisean Index to list of indexes in vector
            polycombstateVec[intersect(presToBreak, polycombpos)] .= 1 # set genes that were ONLY susceptible to prcToBreak as NOT being suppressed by that PRC (genes that were not suppressed don't need to have their polycombstate updated because they are already 1)
            me.polycombstate[:,envState] = polycombstateVec # update polycombstate for given envState with prcToBreak breakdown gene results
        #else
            #print("PRC(s) choosen to break not active ")
        end
    #else
        #print("no PRC that is active ")
    end
end

function breakpoly(me::Individual, envState::Int64)
# function takes an individual's existing polycomb state vector
# and abrogates an entry from 0 to 1 to make exactly one polycomb
# response element inactive (break it's polycomb such that
# polycomb can no longer supress polycomb susceptible gene)

# Original way to breakdown polycomb when only one universal PRC
    polycombstateVec = copy(me.polycombstate[:,envState])
    polycombpos = findall(x->x == 0, polycombstateVec) # find polycomb response element on genes that are active (gene is supressed by polycomb), meaning polycombstate = 0
    if length(polycombpos) > 0 # check that at least one polycomb response element is active before breaking one
        breakidx = rand(1:length(polycombpos))  # break exactly 1 polycomb response element, i.e. make the polycomb response element inactive by stating polycombstate = 1
        polycombstateVec[polycombpos[breakidx]] = 1
        me.polycombstate[:,envState] = polycombstateVec
    end
end

function iterateindForStability(me::Individual, envState::Int64, envStateVecToTest::Vector{Float64}, oldEnvState::Int64, indInOldEnv::Individual)
# iterate individuals from their development state to
# their later state and store the results
# in the Individual object - version to be used with an already set polycomb state vector
# i.e. no iteration to TIMECRIT needed --> so NOT going through development, just testing for stability
# Testing for cell pliancy (so change EnvState but not polycomb state (so typically test
# using individual after went through development state)
    currstate = copy(indInOldEnv.develstate[:,oldEnvState])
    stateupdate = zeros(Float64,G)# Define stateupdate & stateupdateenv before for loop because for
        # loop defines a new scope
    stateupdateenv = zeros(Float64,G)
    vectAD = zeros(Float64,MAXCONV)
    for i=1:MAXCONV
        stateupdate = me.network*currstate # W*Sw; gene interaction matrix times gene state vector
        stateupdateenv = me.EnvIndInter*envStateVecToTest # (W x E)*Se: environment interaction matrix times environment state vector
        stateupdate = stateupdate + stateupdateenv # W*Sw + (W x E)*Se
        stateupdate = 1.0./(1.0.+exp.(-SIGSTR*stateupdate)) # sigma s in eqs
        stateupdate = stateupdate .* me.polycombstate[:,oldEnvState] # turn off polycomb supressed genes that were turned off in old environment individual originally evolved in
        # save vector of length TAU of stateupdate and currstate to input to slidingWindow function
        actualDist = sum(abs.(currstate - stateupdate))/G # normalized Hamming distance between the previous time point state, Sw t (currstate),
                                                         # and the current time point, Sw t+1 (stateupdate)
        # Sliding window - added 7/5/18 - SG
        # if number of iterations once testing for stability (so after TIMECRIT) is >= TAU,
        # then start calling sliding window function to test from (i-(TAU-1)) to i
        vectAD[i] = actualDist
        if i >= TAU
            distBool = (vectAD[i] <= 10.0^(-6.))
            if JACOBIANFLAG
                jacobianMat = jacobian(me,stateupdate)
                eigenVals = eigvals(jacobianMat)
                largestEigVal = maximum(real(eigenVals))
                stableBool = (largestEigVal < 0.0)
            else # else run slidingWindow to test for stability
                stableBool = slidingWindow(i,vectAD)
            end
            if distBool && stableBool # returns true if
                me.stable[envState] = true
                me.develstate[:,envState] = copy(stateupdate) #changed to copy() - SG 6/21/18
                me.pathlength[envState] = i #changed to copy() - SG 6/21/18
                break
            elseif i==MAXCONV
                me.stable[envState] = false
                me.develstate[:,envState] = copy(stateupdate) #changed to copy() - SG 6/21/18
                me.pathlength[envState] = MAXCONV #changed to copy() - SG 6/21/18
            end
        end

        currstate = copy(stateupdate)
    end
end


function disruptPolycomb(ind::Individual, envMoveTo::Int64, evolvedEnv::Int64, founder::Individual, prcsToRemove::Vector{Int64})
    # Generate different individual variables to manipulate in different ways for taking different measurements:
    indWithBrokenPrcs = deepcopy(ind)
    indWithoutBrokenPrcs = deepcopy(ind)
    if BREAKENTIREPRC != true # If break only 1 target element when only 1 PRC:
        breakpoly(indWithBrokenPrcs, 1) # break PRC  for env 1
        breakpoly(indWithBrokenPrcs, 2) # break PRC  for env 2
    else # If multiple PRCs:
        breakPRC(indWithBrokenPrcs, prcsToRemove, 1) # break PRC (prcsToRemove) in env 1
        breakPRC(indWithBrokenPrcs, prcsToRemove, 2) # break PRC (prcsToRemove) in env 2
    end
    # Populations to test in totally different environment than ones evolved in:
    # popBrokenPrcsAndTestNewEnv = deepcopy(indWithBrokenPrcs);
    # popIntactPrcsAndTestNewEnv = deepcopy(indWithoutBrokenPrcs);
    # # Populations to test for stability when start with develstate instead of initstate:
    # indWithBrokenPolycomb = deepcopy(indWithBrokenPrcs);
    # indWithIntactPolycomb = deepcopy(indWithoutBrokenPrcs);
    # Populations to test stability when start with develstate instead of initstate and replace env2 with env1 and vice versus:
    indWithBrokenPolycombReplaceEnvs = deepcopy(indWithBrokenPrcs);
    indWithIntactPolycombReplaceEnvs = deepcopy(indWithoutBrokenPrcs);
    # Iterate for stability when stay in same environemnt evolved in (i.e. start with develstate instead of initstate and stay in same env)
    ### Test for stability when use final develstate as starting point of iterations when testing for phenotypic pliancy instead of initstate:
    ## Testing when stay in same enviroment when polycomb intact and when broken:
    # # In Env 1:
    # iterateindForStability(indWithIntactPolycomb, 1, popRead.founder.EnvState[:,1], 1, ind);
    # iterateindForStability(indWithBrokenPolycomb, 1, popRead.founder.EnvState[:,1], 1, ind);
    # distanceWhenPolycombIntactInEnv1 = meanDistanceBtwPops(indWithIntactPolycomb, pop1, 1, 1);
    # distanceWhenPolycombBrokenInEnv1 = meanDistanceBtwPops(indWithBrokenPolycomb, pop1, 1, 1);
    # Iterate for stability when move to other environment (i.e. start with develstate instead of initstate and replace env2 with env1 and vice versus)
    iterateindForStability(indWithIntactPolycombReplaceEnvs, envMoveTo, founder.EnvState[:,envMoveTo], evolvedEnv, ind) # Intact
    iterateindForStability(indWithBrokenPolycombReplaceEnvs, envMoveTo, founder.EnvState[:,envMoveTo], evolvedEnv, ind) # Broken
    return (indWithBrokenPolycombReplaceEnvs, indWithIntactPolycombReplaceEnvs)
end


function distanceMeasureForPliancy(disruptedInd::Individual, normalInd::Individual, envMoveTo::Int64, oldEnv::Int64)
     if ((disruptedInd.stable[envMoveTo] == true) & (normalInd.stable[envMoveTo] == true) & (normalInd.stable[oldEnv] == true))
         distanceBtwIndForEnvMoveTo = sum(abs.(disruptedInd.develstate[:,envMoveTo] - normalInd.develstate[:,envMoveTo]))/G
         distanceBtwIndAndEnvMoveFrom = sum(abs.(disruptedInd.develstate[:,envMoveTo] - normalInd.develstate[:,oldEnv]))/G
         distanceBtwIndEvolvedEnv1and2 = sum(abs.(normalInd.develstate[:,oldEnv] - normalInd.develstate[:,envMoveTo]))/G
     else
         distanceBtwIndForEnvMoveTo = missing
         distanceBtwIndAndEnvMoveFrom = missing
         distanceBtwIndEvolvedEnv1and2 = missing
     end
     return (distanceBtwIndForEnvMoveTo, distanceBtwIndEvolvedEnv1and2, distanceBtwIndAndEnvMoveFrom)
end


function distBrokenVsIntact(distanceBtwPopsForEnvMoveToIntact::Float64,distanceBtwPopsForEnvMoveToBroken::Float64)
    percentBrokenCloserThanIntactToEnvMoveTo = round(((distanceBtwPopsForEnvMoveToIntact - distanceBtwPopsForEnvMoveToBroken)/distanceBtwPopsForEnvMoveToIntact)*100, digits = 2)
    return percentBrokenCloserThanIntactToEnvMoveTo
end

function functionalPliancy(distanceBtwIndForEnvMoveToBrokenNum::T, distanceBtwIndForEnvMoveToIntactNum::T, distanceBtwIndEvolvedEnv1and2Num::T) where T # could be type Float64 or type Missing
    # Calculate at individual level
    # Remove "missing" entries (when not stable after iterating for stability)
    if !(ismissing(distanceBtwIndForEnvMoveToBrokenNum) & ismissing(distanceBtwIndForEnvMoveToIntactNum))
        if (distanceBtwIndForEnvMoveToBrokenNum < (distanceBtwIndEvolvedEnv1and2Num/2)) & (distanceBtwIndForEnvMoveToBrokenNum < distanceBtwIndForEnvMoveToIntactNum)
            return true
        else
            return false
        end
    else
        return missing
    end
end

function perGeneDistForPliancy(gene::Int64, disruptedInd::Individual, normalInd::Individual, envMoveTo::Int64, oldEnv::Int64)
     distanceBtwIndForEnvMoveTo = abs(disruptedInd.develstate[gene,envMoveTo] - normalInd.develstate[gene,envMoveTo])
     distanceBtwIndAndEnvMoveFrom = abs(disruptedInd.develstate[gene,envMoveTo] - normalInd.develstate[gene,oldEnv])
     distanceBtwIndEvolvedEnv1and2 = abs(normalInd.develstate[gene,oldEnv] - normalInd.develstate[gene,envMoveTo])
     return (distanceBtwIndForEnvMoveTo, distanceBtwIndEvolvedEnv1and2, distanceBtwIndAndEnvMoveFrom)
end

function functionalPliancyPerGene(distanceBtwIndForEnvMoveToBrokenNum::T, distanceBtwIndForEnvMoveToIntactNum::T, distanceBtwIndEvolvedEnv1and2Num::T) where T # could be type Float64 or type Missing
    # Calculate at gene level for given individual
    # If individual is stable then find fraction of genes that are functionally pliant
        # Meaning that when PcG-like mechanism is broken then it moves closer to
        # final state after evolution in environment move to, and also broken moves
        # closer to this state than when PcG-like mechanism is left intact but
        # move to this environment
    if (distanceBtwIndForEnvMoveToBrokenNum < (distanceBtwIndEvolvedEnv1and2Num/2)) & (distanceBtwIndForEnvMoveToBrokenNum < distanceBtwIndForEnvMoveToIntactNum)
        return true
    else
        return false
    end
end

function functionalPliancyPerIndividual(ind::Individual, founder::Individual, envMoveTo::Int64, oldEnv::Int64, prcsToRemove::Vector{Int64})
    disruptPcGResults = disruptPolycomb(ind, envMoveTo, oldEnv, founder, prcsToRemove)
    disruptedBrokenInd = disruptPcGResults[1]
    disruptedIntactInd = disruptPcGResults[2]
    #geneFuncPliancyResults = fill(false, G)
    geneFuncPliancyResults = Array{Union{Missing, Bool}}(missing, G)
    indStability = true
    if ((disruptedBrokenInd.stable[envMoveTo] == true) & (disruptedIntactInd.stable[envMoveTo] == true) & (ind.stable[envMoveTo] == true) & (ind.stable[oldEnv] == true))
        for gene = 1:G
            distOfGeneBroken = perGeneDistForPliancy(gene, disruptedBrokenInd, ind, envMoveTo, oldEnv)
            distOfGeneBrokenToEnvMoveTo = distOfGeneBroken[1]
            distOfGeneEnv1ToEnv2 = distOfGeneBroken[2]
            indStability = true
            distOfGeneIntactToEnvMoveTo = perGeneDistForPliancy(gene, disruptedIntactInd, ind, envMoveTo, oldEnv)[1]
            geneFuncPliancyResults[gene] = functionalPliancyPerGene(distOfGeneBrokenToEnvMoveTo, distOfGeneIntactToEnvMoveTo, distOfGeneEnv1ToEnv2)
        end
    else
        indStability = false
    end
    if indStability
        indFunctionalPliancyResults = length(findall(geneFuncPliancyResults))/length(geneFuncPliancyResults)
        return indFunctionalPliancyResults
    else
        return missing
    end
end


function functionalPliancyForPopulation(distanceBtwPopsForEnvMoveTo::Vector{Float64}, distanceBtwPopsEvolvedEnv1and2::Vector{Float64})
    # Remove "missing" entries (when not stable after iterating for stability)
    distanceBtwEnvMoveTo = distanceBtwPopsForEnvMoveTo[findall(!ismissing, distanceBtwPopsForEnvMoveTo)]
    distanceBtwEnv1and2 = distanceBtwPopsEvolvedEnv1and2[findall(!ismissing, distanceBtwPopsEvolvedEnv1and2)]
    if length(distanceBtwEnvMoveTo) >= 1 # meaning at least 1 individual was stable in both environments after removing missing entries
        return round(length(findall(y -> distanceBtwEnvMoveTo[y] < (distanceBtwEnv1and2[y]/2), 1:length(distanceBtwEnv1and2)))/length(distanceBtwEnv1and2)*100, digits = 5)
    else
        return 0.0
    end
end
