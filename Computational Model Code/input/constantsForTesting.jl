# Parameters for polycomb model (Maryl's version)
const INDTRIALS = 2 # the number of different evolutionary trajectories for the same starting population
# ***update below***
# const INDTRIALSFOLDER = "savedPopsForIndTrialsWithSlidingWinAndActualDistAndSigstrOf6" # name of folder to use that contain the .jld
#     # files for the independent populations previously generated using generateSavePopulations() function to do data analysis on full dataset generated for the paper
const G = 50 # number of genes per individual, default 10, small 3 # use 50 for paper
const N = 500 # population size, default 500, small 10 (N individuals) # use 1000 for paper
const C = 0.1 # connectivity probability for individuals' W (Gaussian matrices); network where each gene is on average connected to C of the other genes
const GENS = 1000 # total number of generations, default 400, small 10; use 1000 for paper
const ENVS = 50 # total number of possible environmental states an individual can encounter during evolution and lifetime (length of EnvState vector) # use 50 for paper
const MAXCONV = 100 # max number of iterations to test for convergence for individuals and founder
const MUTRATE = 0.1 # mutation rate used in MUTRATE/(cG^2)
const ENVMUTRATE = 0.1 # environment individual interaction mutation rate used in "mutateEnvIndInter" function
const MUTMAG = 1 # magnitude of mutations; also used when generating founder; standard deviation of guassian drawing random numbers from
const SELSTR = 0.5 # selection strength (sigma in fitness equation), highest selection when close to 0, no selection as SELSTR approaches infinity
const P = 0.5 # connectivity probability for population structure graph
const TAU = 10 # Sliding window length to take average of expression levels to compare to find stability
                #look-behind depth for convergence testing with sigmoidal function, default 10
const SIGSTR = 1.0 # sigmoid stretch ("a" in "sigma s" equation)
const POLYMUTRATE = 0.1 # "T" in Saurabh paper if divide this by G; mutation rate for vector of genes that are susceptible to polycomb
const GENETHRESH = 0.15 # critical threshold to measure if polycomb susceptible genes are expressed or not (gamma); 0.05 originally
const RANDSEEDNUM = 6 # Number to use to set random seed for the starting population
const PRCS = 2 # number of different types of PcG (Polycomb group protein) repressive complexes (PRCs) that can function in each individual
randNum = 1 # number to create RANDSEEDVEC from below # can also make rand(1:1000) if want it to be randomly selected; for data generation we had all the different starting populations run with the same random seed vector (so go through same evolutionary trajectory but different starting populations)
const RANDSEEDVEC = collect(randNum:(INDTRIALS+randNum+1))# If want random vector of random seeds to use then use: collect(randNum:(INDTRIALS+randNum+1)). If want to directly input random seed vector used in past run of code then manually input that vector that is saved in "resultsWriteUp.txt". # random seed to use for each independent trial
const PRCTOTAL = sum(map(x -> binomial(PRCS,x), 1:PRCS)) # number of different combinations of PRC(s) that can affect one gene
const MINDISTFORPLIANCY = 0.05 # minimum distance between gene expression patterns for individual to be considered phenotypically pliant
const MINDISTFORPLIANCYNEWENV = 0.1 # minimum distance between gene expression patterns for individual to be considered phenotypically pliant when introduce to totally new environment
const DIFENVS = 2 # number of different environments to run during evolution
# Extended Model Environmental constants:
const CENVEC = 0.4 # proportion of environment states set to 1 (meaning on and affecting individual)
const CENMAT = 0.4 # proportion of gene's able to be affected by any environments (meaning all environment elements have potential to affect particular gene) --> probably should be even lower
const CENGENE = 0.4 # proportion of environments that affect a given gene: connectivity of row in W x E -->
# # Robustness testing:
const ROBIT = 100 # number of iterations to run in genetic robustness testing
const ALTERENV = 0.04 # when mutating an environment, alter this fraction of entries of environment state vector, Se --> when G & ENVS = 50, want to alter 2 entries then ALTERENVS = 0.04
const ENVCHANGES = 100 #binomial(G,convert(Int64,ALTERENV*G)) # total number of different mutated environments want to test with intact and broken polycomb analysis
#                         # AND number of different environments to test when measuring the average environmental robustness for one individual during evolution
#                         # Total number of possible combinations of G choose ALTERENV*G = binomial(G,convert(Int64,ALTERENV*G))
const DIST = 0.05 # maximum distance two development state vectors, Sw, to be considered the same (stable)
const EPSI = 10.0^(-6.0) # epsilon to use when testing for stability with the slidingwindow() function,
    # such that individual is considered stable if normalized variance is less than epsilon
const WMAT = Array{Float64,2}(undef,0,0) # to set W matrix for founder to generate population
const INP = (Float64)[] # initial state input --> to set initstate vector for founder to generate population
const OPT = (Float64)[] # initial opt state of founder --> if need to set optimum state for founder

# For two input - output relationships where evolve to desired output:
const PERCENTDIFENVS = 0.70 # percent of entries to change in envState vector to then use to set envState2 vector
const PERCENTDIFGENES = 0.50 # percent of entries to change in optstate vector to then use to set optState2 vector such that optState2 is at least 40% different than optstate
const PERCENTDIFOPTSTATES = 0.40 # percentage of minimum difference between optstate and optState2 vectors when randomly setting optState2 vector to evolve to during evolution
const MAXENVINPUT = 1000 # maximum number of times to test for finding a stable individual with a randomly choosen envState2 vector that is 70% different than envState input
const MAXOPTSTATES = 1000 # maximum number of times to test for finding randomly choosen optState2 vector that is at least 40% different than optstate vector
const MAXINITSTATES = 3
const PEAK = 1 # peak of sum of distances between 2 environment develstates and their optstates; used when measuring fitness and adjusting for low fitness of individual in env 2 at start of evolution
const MINFIT = 0.45 # where distances cross the x-axis, so want to be a little less than minimum of distances between 2 environment develstates and their optstates; used when measuring fitness and adjusting for low fitness of individual in env 2 at start of evolution

# For measurements:
const SAVEDATA = (0, 25, 50, 75, 100, 500, 1000) # Generations to save population data; If testing 10,000 Generations: (0, 25, 100, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000, 6500, 7000, 7500, 8000, 8500, 9000, 9500, 10000)
const TOTALMEAS = 20 # Total number of times going to make measurements
const FRACMEAS = TOTALMEAS/GENS # parameters for when to make measurements of population properties; function called "measure" that uses these.
const MEASPERIOD = 1/FRACMEAS # tells how many generations to skip before making new measurements of population properties (i.e. how often to make measurements of population properties)
const NUMMEAS = convert(Int64,round((GENS/MEASPERIOD) + 1)) # Total number of times going to make measurements of population properties plus 1: total number of generations divided by how many generations to skip between making measurements for population properties plus 1 # Used to keep track of measurements making in main.jl

# Flags:
const RANDSEEDFLAG = true # if want to use a random.seed or not
const RUNWITHANDWITHOUTPOLYCOMB = false # if running bash script where generate results with and without polycomb to include in same dataOutput folder to then use to plot both results on same graph in R
const SAVEDPOPSFLAG = false # if true, then use previously created and saved populations when run "runMultipleIndependentTrails.jl"
const SAURABVERSION = false # if true, then Saurab's version of code (same methods used in paper -> Epigenetics Decouples Mutational from Environmental Robustness. Did it also facilitate multicellularity?) is used such that environmentalStateVec is set to zero, and change the environment by altering initstate vector
const MUTATEENVINTERACTIONFLAG = true # if true, then the environment interaction with network matrix can be mutated during evolution. False, if it is fixed throughout evolution.
const SAMESTARTINGPOPFORINDTRIALSFLAG = true # if true, then each independent trial starts with the same population but has a different random seed, so the evolution will be different but same starting population. This is better when comparing different independent trials (populations) in R using PCA analysis
const SYMMETRIC_STARTWMAT = false # if true, then generates symmetric W matrix for founder individual
# Polycomb flag:
const POLYCOMBFLAG = true #if true, use polycomb procedure
if POLYCOMBFLAG # if true, then run code with evolution of Polycomb Repressive Complexes (PRCs), and then set TIMECRIT to be 2
    const TIMECRIT = 2
else
    const TIMECRIT = 0
end
const JACOBIANFLAG = false
const INDWEIGHTS = "gaussian" # individual weights; if "discrete" then sample {1, -1, 0} via rand(-1:1) if "gaussian" use randn()
const SETFOUNDERFLAG = false # if true, then imports founder from previous run of code to use as the founder for this run of the code as well
const RANDPOP = false # if true generate initial population with random interactions rather than a homogeneous one
const SEXUALREPRO = true # if true interleave rows from two individuals each generation
# For using IPA networks:
const IPANETWORKS = false # boolean that tells you if you should generate W matrix for an individual using a given IPA generated file containing gene regulatory network relationships
const IPAFILENAME = "TFAP2Anetwork.txt" # file name of IPA generated gene regulatory network relationships file
# # For Phenotypic Pliancy
const BREAKPOLY = false #breakdown the polycomb mechanism after the evolution process
const BREAKENTIREPRC = true # if want to breakdown entire an entire given PRC (==true), OR just want to break one gene that is controlled by given PRC to break (==false)
# # Measurement Flags:
const MEASUREROBUST = true # measure robustness
const MEASUREFIT = true # measure fitness
