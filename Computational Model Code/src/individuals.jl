# ---------------------------
# CONSTRUCTORS
# ---------------------------
function Individual(network::Matrix{Float64}, initstate::Vector{Float64})
    Individual(copy(network), copy(initstate), zeros(Float64,(G,ENVS)), zeros(Float64,ENVS,DIFENVS),
                zeros(Float64,G,DIFENVS), zeros(Float64,G,DIFENVS), falses(DIFENVS), 0., zeros(Float64,DIFENVS), zeros(Int64,DIFENVS),
                    zeros(Float64,G,DIFENVS), ones(Float64,G,DIFENVS))#, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.)
    # Individual composed of:
    # Individual(network, initstate, EnvIndInter, EnvState,
    #         develstate, optstate, stable, fitness, fitnessUnderEachEnv,
    #           pathlength, polycombvec, polycombstate)
    #           , connectivityIngoing,
    #           connectivityOutgoing, btwCentrality, cliqueSize,
    #           connectivityIngoingTargets, connectivityOutgoingTargets,
    #           btwCentralityTargets, cliqueSizeTargets, connectivityIngoingTargetsRepressed,
    #           connectivityOutgoingTargetsRepressed, btwCentralityTargetsRepressed,
    #           cliqueSizeTargetsRepressed, phenotypicPliancyScore1,
    #           phenotypicPliancyScore2, newEnvPliancyScore1, newEnvPliancyScore2)

end

function Individual(network::Matrix{Float64})
    Individual(copy(network), zeros(Float64,G))
end

function Individual()
    Individual(zeros(G,G))
end


# ---------------------------------------------------------
# METHODS: Generate Individuals for the initial Population:
# ---------------------------------------------------------
function geninds(founder::Individual)
# If given an individual to use as the founder:
# Initialize a Population by creating
# N individuals that are all a copy of founder
    indvect1 = Array{Individual{Float64}}(undef, N) # Julia assigns Individual type as whole as Float64
    for i = 1:N
        if RANDPOP # generate random population that is NOT a copy of the founder
            indvect1[i] = stableind(founder.initstate)
        else # All individuals in population are copy of a founder
            indvect1[i] = deepcopy(founder)
        end
    end
    return indvect1
end


function geninds()
# If not given an individual to use as the founder:
# Generate N individuals by generating a
# stable founder first and then copy this
# founder N times
    founder = genfounder()

    indvect1 = Array{Individual{Float64}}(undef, N) # Julia assigns Individual type as whole as Float64
    for i = 1:N
        if RANDPOP # generate random population that is NOT a copy of the founder
            indvect1[i] = stableind(founder.initstate)
        else # All individuals in population are copy of a founder
            indvect1[i] = deepcopy(founder)
        end
    end
    return indvect1
end

# ----------------------------------------------
# METHODS: Generate Founder Step
# ----------------------------------------------
function genfounder()
# generate a founding individual whose
# developmentally stable state is equivalent to
# its optimal state (however rand
# assumes user not putting in optimal state vector)
# and thereby determines the optimal state for the population,
# i.e. the optimal state is set as the final state of founder
# once stability is reached within MAXCONV iterations
    if isempty(INP) & isempty(OPT)# **this is what currently runs
        founder = randstableind()
        # print(" initstate & optstate empty")
    elseif isempty(OPT) #this what used to run bc INP was set in constants (before 7/25/18)
        founder = stableind(INP) # founder should be stable at this point
        # print(" inistate NOT empty")
    elseif isempty(INP)
        # print(" optstate NOT empty")
        founder = randstableind()
    else
        # print(" initstate & optstate NOT empty")
        founder = stableind(INP)
    end
    return founder
end

function randstableind() # set Sw vector
# generate a random Individual that is
# developmentally stable by initializing the state vector as rand(0:1,G)*2-1 and passing through stableind function
# initialzing the initial state vector and generating a stable individual
    randind = stableind(convert(Vector{Float64},rand(0:1,G))) # assumes user not inputting optstate into stableind and then generates initial state
    return randind
end

function stableind(initstate::Vector{Float64}) # used if only do not input optstate into function; optstate will be empty vector
    if isempty(OPT) & isempty(WMAT) #this is what usually runs as of Aug 2019 because W & optstate are empty in constants
        ind = stableind(initstate,(Float64)[],Array{Float64,2}(undef,0,0))
        # print(" wmat & optstate NOT given")
    elseif isempty(OPT)
        ind = stableind(initstate,(Float64)[],WMAT)
        # print(" optstate NOT given, w mat given")
    elseif isempty(WMAT)
        ind = stableind(initstate,OPT,Array{Float64,2}(undef,0,0))
        # print(" wmat NOT given, optstate given")
    else
        ind = stableind(initstate,OPT,WMAT)
        # print(" given optstate & wmatrix")
    end
    return ind
end


function stableind(initstate::Vector{Float64},optstate::Vector{Float64},wmatrix::Matrix{Float64})
# This function is used to generate a founder that is stable:
# generate a random Individual that is developmentally stable by
# giving an initial individual state vector and can optionally give
# optimal state vector if desire. This function generates an individual's
# Gene Regulatory Network (W), Environment State vector (Se) and
# Environment Gene interaction matrix (W x E).
    if isempty(wmatrix)
        randind = Individual() # generating object of type individual
        randind.initstate = copy(initstate)
        # print(" wmat NOT given in stableind(inistate,optstate,wmat) function")
    else
        # print(" wmat given in stableind(inistate,optstate,wmat) function")
        randind = Individual() # generating object of type individual
        randind.initstate = copy(initstate)
        randind.network = copy(wmatrix)
    end
    # what happens if this wmatrix and initstate isn't stable for saurabh version but for extended?? 3/21/19 -ML
    # so maybe always want to generate Saurabh population first and use this Wmatrix to find stable extended population founder
    # b/c can alter WxE interaction matrix and/or Se vector when searching for stability with given initstate vector and W matrix
    countsForInitstates = 1
    while (randind.stable!=trues(DIFENVS)) & (countsForInitstates < MAXINITSTATES)
        while ((randind.stable[1]!=true) | (if isempty(optstate); false; else; randind.develstate!=optstate; end))
            # right in OR statment is only for when user inputs an optstate,
            # but in PcG model optstate is always empty because set as such in constants
            if isempty(wmatrix) # if not W matrix (network) given, then generate random one
                randind.network = zeros(Float64,G,G)
                if IPANETWORKS == false
                    for i=1:G^2
                        if rand()<C # connectivity probability for individuals' W (Gaussian matrices)
                            if INDWEIGHTS == "discrete" # individual weights; if "discrete" then sample {1, -1, 0} via rand(-1:1)
                                #randind.network[i] = [-1,1][rand(1:2)]
                                randind.network[i] = rand(-1:1)
                            elseif INDWEIGHTS == "gaussian" # if "gaussian" then sample -inf to inf (is this correct??? 2/28/18 -ML --> yes) via randn()
                                randind.network[i] = randn()
                            elseif INDWEIGHTS == "grnInference"
                                randind.network[i] = rand([-1,1])
                            else
                                error(string("wrong specification of individual weights:",
                                       " use discrete or gaussian"))
                            end
                        end
                    end
                    if SYMMETRIC_STARTWMAT == true
                        for i=1:G
                            for j=1:G
                                # make diagonal not zero so will not end up with 0.5 expression after 1 iteration through dynamics in grn inference network building to tests
                                if j==i
                                    randind.network[i,j] = rand([-1,1])
                                    # make second layer diagonal to make sure network is fully connected:
                                    # if j < G
                                    #     randind.network[i,(j+1)] = rand([-1,1])
                                    # end
                                end
                                # make symmetric
                                randind.network[j,i] = randind.network[i,j]
                            end
                        end
                    end
                else
                    randind.network = generateWmatUsingIpaNetwork(IPAFILENAME)
                end
            end

            # Initialize Environment Gene Network Interaction matrix (W x E)
            randind.EnvIndInter = zeros(Float64, G, ENVS)
            if !SAURABVERSION # if running Saurab's older version of PcG, then do not have environmentalStateVec & environmentInteractionMatrix (thus equal to zero, which are set in Individual function above)
                for i=1:G
                    if rand()<CENMAT # proportion of environments that can possibly affect a gene(s)
                        for j=1:ENVS
                            if rand()<CENGENE # proportion of genes that can be affected by a given environment
                                randind.EnvIndInter[i, j] = randn() # effect of environment j on gene i
                            end
                        end
                    end
                end
            end

            # Initialize Environment State vector (environment present or not)
            if SAURABVERSION # if running Saurab's older version of PcG, then do not have environmentalStateVec (thus equal to zero)
                randind.EnvState = zeros(Float64,ENVS,DIFENVS)
            else
                envState1 = convert(Vector{Float64},rand(0:1,ENVS))
                envStatesVectors = map(x -> envState1, 1:DIFENVS)
                randind.EnvState = hcat(envStatesVectors...) # ***check what the ... does
            end
            iterateindfound(randind, 1)
        end

        # For different environmental states:
        # First find other environmental inputs (EnvState) that are stable for randind
        # but 70% different than envState of randind.envState[:,1]
        # assuming that initstate are the same for both environmental state cases
        if DIFENVS > 1
            for i = 2:DIFENVS
                counts = 1
                while (randind.stable[i] == false) & (counts < MAXENVINPUT)
                    EnvState2ToChange = randperm(ENVS)[1:convert(Int64,PERCENTDIFENVS*ENVS)]
                    randind.EnvState[EnvState2ToChange, i] = 1 .- randind.EnvState[EnvState2ToChange, 1]
                    iterateindfound(randind, i) # develstate, stable, pathlength, and polycombstate get set in this function for this particular environmental state, i
                    counts += 1
                end
            end
        end
        countsForInitstates += 1
    end
    if countsForInitstates >= MAXINITSTATES
        randind = stableind(convert(Vector{Float64},rand(0:1,G)))
    end
    # set optstate
    if isempty(optstate)
        randind.optstate[:,1] = copy(randind.develstate[:,1]) # set optstate of first environment so this environment is one that is already optimum
    else # if user inputs optstate then develstate and optstate will not be equal
        randind.optstate = copy(optstate)
    end
    # Find and set optstate(s) in other environmental states that is at least 40% different than founder opstate in environment state 1
    if DIFENVS > 1
        for j = 2:DIFENVS
            distanceOptStates = 0.0
            countNum = 1
            while (distanceOptStates < PERCENTDIFOPTSTATES) & (countNum < MAXOPTSTATES)
                randind.optstate[:,j] = copy(randind.optstate[:,1])
                optState2ToChange = randperm(G)[1:convert(Int64,PERCENTDIFGENES*G)]
                randind.optstate[optState2ToChange,j] = abs.(rand() .- randind.optstate[optState2ToChange,1])
                distanceOptStates = sum(abs.(randind.optstate[:,1] - randind.optstate[:,j]))/G
                countNum += 1
            end
            if countNum > MAXOPTSTATES
                randind = stableind(initstate)
            end
        end
    end
    # calculate fitness
    fitnesseval(randind)
    fitnessUnderEachEnvEval(randind)
    return randind
end


function iterateindfound(me::Individual, envState::Int64)
# iterate founder individual from its initial state to
# their developmental state and store the results
# in the Individual object under only the first environmental state
# This function used to find a developmentally stable founder to use for initial population
    currstate = copy(me.initstate)
    vectAD = zeros(Float64,MAXCONV)
    for i=1:MAXCONV
        stateupdate = me.network*currstate
        stateupdateenv = me.EnvIndInter*me.EnvState[:,envState]
        stateupdate = stateupdate + stateupdateenv
        stateupdate = 1.0./(1.0 .+ exp.(-SIGSTR.*stateupdate))
        tempdiff = currstate - stateupdate
        actualDist = sum(abs.(tempdiff))/G

        #sliding window - added 7/5/18 - SG
        vectAD[i] = actualDist
        if i >= TAU
            distBool = (vectAD[i] <= 10.0^(-6.0))
            stableBool = slidingWindow(i,vectAD)
            if distBool && stableBool
                me.stable[envState] = true
                me.develstate[:,envState] = copy(stateupdate) #changed to copy() - SG 6/21/18
                me.pathlength[envState] = copy(i) #changed to copy() - SG 6/21/18
                break
            elseif i==MAXCONV
                me.stable[envState] = false
                me.develstate[:,envState] = copy(stateupdate) #changed to copy() - SG 6/21/18
                me.pathlength[envState] = copy(MAXCONV) #changed to copy() - SG 6/21/18
            end
        end
        currstate = copy(stateupdate)
    end
end


function slidingWindow(t::Int64, x::Array{Float64,1})
# Measures if individual has reached stability using a sliding window of t-(TAU-1) to t and
# taking the normalized Hamming distance of the measurements recorded and tests if they
# fulfill the criteria of stability - i.e. the average Hamming distances in a sliding
# window of length TAU is less than epsilon, which equals 10^-4.
# Written by Sarah Garaff July 5, 2018
    distances = 0.0
    for i = (t-(TAU-1)):t
        distances += x[i] #add up each of the distances between measurements
    end
    std = distances/TAU #find std by divding by TAU
    if std <= EPSI # stable if normalized variance is less than epsilon (EPSI), which equals 10^(-6)
         return true # it's stable
    else
        return false
    end
end



# -------------------------------------------
# METHODS: Evolution Step
# -------------------------------------------
function update(mes::Vector{Individual{T}}) where T
# This function updates the state of a vector (Sw) of individuals
# ****One run of "update" function corresponds to one generation in evolution****

    # pmap runs in parallel if julia is invoked with multiple threads
    #pmap(develop, mes) # Run the "develop" function for the given input individual
    #*** Why running develop function for every generation (each iteration of update
    # function) because at end of while loop, will have replaced every individual
    # in the population and these individuals have undergo development during
    # while loop, so essentially developing individuals twice??? Except in the case
    # of the first generation development step - ML question 11/13/17
    oldinds = deepcopy(mes)

    newind = 1 # new individual who is stable and with
    # high enough fitness will replace the old individual

    while newind <= length(mes)
        #tempind = Individual()

        z = rand(1:N) # Pick a random individual from population
                    # with N individuals total to undergo asexual
                    # or sexual reproduction

        if SEXUALREPRO
            r = rand(1:N) # SEXUALREPRO is false for polycomb code and analysis --> not after 2018
            tempind = deepcopy(oldinds[z])
            reproduce(oldinds[z], oldinds[r], tempind)
        else
            tempind = deepcopy(oldinds[z])
        end

        mutate(tempind) # possibly mutate individuals gene network (W)

        if !SAURABVERSION & MUTATEENVINTERACTIONFLAG
            mutateEnvIndInter(tempind) # possibly mutate individuals environment gene network interaction matrix (W x E)
        end

        if POLYCOMBFLAG # True for polycomb code and analysis
            polymut(tempind)
        end

        # Individual goes through development in environment 1
        iterateIndForOneEnvState(tempind, 1)

        if tempind.stable[1] # select for developmental stability for individual in first environment because if not stable in the first environemnt then don't need to run in the second since needs to be stable in both environments to be selected for the next generation
            if DIFENVS > 1
                for j = 2:DIFENVS
                    iterateIndForOneEnvState(tempind, j)
                end

                if tempind.stable == trues(DIFENVS) # select for developmental stability for individual under all different environmental states
                    # If any is not stable, then don't keep
                    if MEASUREFIT
                        fitnesseval(tempind) # set to true in polycomb analysis
                        fitnessUnderEachEnvEval(tempind) # measure fitness under each environmental condition separately
                    else
                        tempind.fitness = 1.
                        tempind.fitnessUnderEachEnv = ones(DIFENVS)
                    end

                    if tempind.fitness >= rand() # Stabilizing selection step:
                            # higher fitness means more likely to
                            # have reproductive success because more likely to get picked
                            # out of population --> rand() produces numbers 0 to 1 with
                            # Uniform distribution, therefore the larger the number, the
                            # more likely it will be greater than the rand() number
                        mes[newind] = deepcopy(tempind)
                        newind += 1
                    end
                end
            end
        end
    end
end


# ------------------------------------------
# METHODS: Development Step
# ------------------------------------------
function iterateIndForOneEnvState(me::Individual, envState::Int64)
# iterate individuals from their initial state to
# their developmental state and store the results
# in the Individual object
    currstate = copy(me.initstate) # currstate = Sw vector; initstate does NOT get changed
    stateupdate = zeros(Float64,G)# Define stateupdate & stateupdateenv before for loop because for
        # loop defines a new scope
    stateupdateenv = zeros(Float64,G)
    vectAD = zeros(Float64,MAXCONV) #vector to 'save' the actualDist
    if POLYCOMBFLAG # If zero then not running with no polycomb
        for i=1:TIMECRIT
            stateupdate = me.network*currstate # W*Sw; gene interaction matrix times gene state vector
            stateupdateenv = me.EnvIndInter*me.EnvState[:,envState] # (W x E)*Se: environment interaction matrix times environment state vector
            stateupdate = stateupdate + stateupdateenv # W*Sw + (W x E)*Se
            stateupdate = 1.0./(1.0.+exp.(-SIGSTR*stateupdate)) # converting resulting state vector into 0's or 1's --> if vector contains only -1's and 1's then will not alter vector
            currstate = copy(stateupdate)
            actualDist = sum(abs.(currstate - stateupdate))/G # normalized Hamming distance between the previous time point state, Sw t (currstate),
                                                             # and the current time point, Sw t+1 (stateupdate)
            vectAD[i] = actualDist
        end

        # run after timecrit
        polytemp = findall(x-> x < GENETHRESH, stateupdate) # find genes that are less than or equal to a gene threshold, GENETHRESH
        prcPresent = sum(me.polycombvec, dims = 2) # find if gene is susceptible to *any* PRCs
        prcPresentPos = findall(x -> x > 0, prcPresent) # returns Cartisean coordinates so must transform into vector of gene numbers that are susceptible; do this is step below --> polys
        # to turn from Cartisean Index to list of indexes in vector:
        polys = map(x -> prcPresentPos[x][1], 1:length(prcPresentPos)) # polys tells which genes are susceptible to at least 1 PRC
        # Theta vector == polycombvec. Genes possibly susceptible to polycomb when
        # polycombvec = 1, 2, ..., PRCTOTAL: for example, gene has polycomb response element for either PRC1, PRC2 or both.
        # Thus, polycombstate entry needs to be set to 0 if this gene's expression is less than a gene threshold and
        # polycomb protein is active (assumed to be active until explicitly break PcG protein).
        polycombind = intersect(polytemp, polys) # Genes to be suspressed by PRC --> Find index of genes that are susceptible to
        # polycomb elements (have a PRE & thus polycombvec >= 1) and less than or equal to a gene threshold after
        # critial time in development (TIMECRIT)
        me.polycombstate[:,envState] = ones(Int64,G) # need this step to reset the past polycomb control and reset for the new polycomb control for this generation based on evolution of genes' susceptiblility to PRCs and PRC evolution for this particular envState (environment)
        # (meaning polycomb response element (PcG protein) is inactive) bc set to 1
        me.polycombstate[polycombind, envState] .= 0 # polycombstate entry is 0 if polycomb susceptible
        # gene is turned off by polycomb response element, i.e. the gene's corresponding
        # polycomb response element is active for given environmental state, envState (i.e. column of polycombstate)
    end

    for i=(1+TIMECRIT):MAXCONV
        stateupdate = me.network*currstate # W*Sw; gene interaction matrix times gene state vector
        stateupdateenv = me.EnvIndInter*me.EnvState[:,envState] # (W x E)*Se: environment interaction matrix times environment state vector
        stateupdate = stateupdate + stateupdateenv # W*Sw + (W x E)*Se
        stateupdate = 1.0./(1.0.+exp.(-SIGSTR*stateupdate)) # sigma s in eqs
        stateupdate = stateupdate .* me.polycombstate[:,envState] # turn off polycomb supressed genes
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


function reproduce(me::Individual, you::Individual, us::Individual) # Do not do sexual reproduction in polycomb stuff because individual cells behave more like asexual reproduction
# sexual reproduction via independent row segregation

    reproindxs = rand(0:1,G)

    for i in findall(x->x==1,reproindxs)
        us.network[i,:] = deepcopy(me.network[i,:])
    end

    for j in findall(x->x==0,reproindxs)
        us.network[j,:] = deepcopy(you.network[j,:])
    end

    return reproindxs
end


function mutate(me::Individual)
# mutate nonzero elements of an individuals
# network according to a rate parameter normalized
# by the size of the nonzero entries in an
# individual network
# done for each reproduction step for each generation

    # Find the non-zero entries as potential mutation sites
    nzindx = findall(x->x!=0,me.network)
    # Determine the connectivity of the matrix
    # by counting the number of non-zeros:
    cnum = length(nzindx) # note cnum=C*G^2, where c is the connectivity of network (W)
    mutflag = 0
    mutatedElement = 0
    for i=1:cnum # chance to mutate any of the non-zero elements in W matrix
    # with a total of 10% chance of a *single* mutation per individual
        if rand() < MUTRATE/cnum # about MUTRATE/cnum percent chance that
        # will mutate this particular non-zero element of W matrix
            me.network[nzindx[i]] = MUTMAG*randn() # Mutate a non-zero entry by a number chosen from a gaussian
            mutflag =+ 1
            mutatedElement = nzindx[i]
        end
    end
    return [mutflag, mutatedElement]
end


function mutateEnvIndInter(me::Individual)
# Function mutates nonzero elements of an individuals
# Environment Individual Interaction matrix (W x E)
# according to a rate parameter normalized
# by the size of the nonzero entries in an
# individual environment Interaction matrix.
# Function is ran for each reproduction step
# for each generation

    # Find the non-zero entries as potential mutation sites
    nzindx = findall(x->x!=0,me.EnvIndInter)

    # Determine the connectivity of the matrix
    # by counting the number of non-zeros
    cnum = length(nzindx)
    envmutflag=0
    mutatedElement=0
    for i=1:cnum
        # For each non-zero entry:
        # With probability Rate/C*G^2, note cnum=(C*(G^2))
        if  rand() < ENVMUTRATE/cnum
            # Mutate a non-zero entry by a number chosen from a gaussian
            me.EnvIndInter[nzindx[i]] = ENVMUTRATE*randn()
            envmutflag =+ 1
            mutatedElement = nzindx[i]
        end
    end
    return [envmutflag, mutatedElement]
end


function polymut(me::Individual)
# perform a one-site mutation of a given
# individual's polycomb vector
# to get population with individuals that can have
# different polycomb response elements inactivated
    mutmat = copy(me.polycombvec)
    for i = 1:G
        if rand() < POLYMUTRATE/G #10% of polycomb susceptible genes can be mutated and there are G polycomb susceptible genes
            prcToMutate = rand(1:DIFENVS) # select with PRC to mutate, so column to mutate because columns of polycombvec are PRCs when have more than one PRC
            mutmat[i,prcToMutate] = 1-mutmat[i,prcToMutate] # change zero entry to 1, or 1 entry to zero;
        end
    end
    me.polycombvec = mutmat
    return mutmat
end




#--------------------------------------------------------------
# MEASUREMENTS FOR FITNESS, ROBUSTNESS, AND MODULARITY
#--------------------------------------------------------------
# function fitnessevalWithoutAdjustment(me::Individual)
# # measure fitness according to the distance
# # between the developmental state determined
# # by running iterateIndForOneEnvState(me) and the optimum
# # state (founder development state) for the
# # two individuals without adjusting for fact that fitness
# # will be higher in first environment to start and adjusting
# # for this (for fitness with adjustment, see fitnesseval() function below)
#     if me.stable == trues(DIFENVS)
#         distanceForFitness = sum(map(x -> sum(abs.(me.develstate[:,x] - me.optstate[:,x])), 1:DIFENVS))/(G*DIFENVS) # ranges from 0 to infinity
#         fitnessResult = exp.(-(distanceForFitness/SELSTR)) # me.fitness ranges from 1 to
#         # zero, with 1 being the most fit (i.e. distance = 0)
#         me.fitness = fitnessResult
#     else
#         me.fitness=0 # fitness set to zero when individual does not reach
#         # developmental stability (steady state by MAXCONV iterations)
#     end
# end


function fitnesseval(me::Individual)
# measure fitness according to the distance
# between the developmental state determined
# by running iterateIndForOneEnvState(me) and the optimum
# state (founder development state) for the
# two individuals
    if me.stable == trues(DIFENVS)
        sumOfDistancesForFitness = sum(map(x -> sum(abs.(me.develstate[:,x] - me.optstate[:,x])), 1:DIFENVS))
        distanceForFitness = sumOfDistancesForFitness/(G*DIFENVS) # ranges from 0 to infinity
        fitnessResult = (exp.(-(distanceForFitness/SELSTR)) + maximum([0,(PEAK-(PEAK/MINFIT)*(distanceForFitness))]))/2 # me.fitness ranges from 1 to
        # zero, with 1 being the most fit (i.e. distance = 0)
        me.fitness = fitnessResult
    else
        me.fitness=0 # fitness set to zero when individual does not reach
        # developmental stability (steady state by MAXCONV iterations)
    end
end


function fitnessUnderEachEnvEval(me::Individual)
    if me.stable == trues(DIFENVS)
        for i = 1:DIFENVS
            distanceForFitness = sum(abs.(me.develstate[:,i] - me.optstate[:,i]))/G # ranges from 0 to infinity
            fitnessResult = exp.(-(distanceForFitness/SELSTR)) # me.fitness ranges from 1 to
            # zero, with 1 being the most fit (i.e. distance = 0)
            me.fitnessUnderEachEnv[i] = fitnessResult
        end
    else
        me.fitnessUnderEachEnv=zeros(Float64,DIFENVS) # fitness set to zero when individual does not reach
        # developmental stability (steady state by MAXCONV iterations)
    end
end
