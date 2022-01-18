#-------------------------------------------------------------------------------
# PCA plot for a single population:
function plotPcaForSinglePopAndEvolTraj(popSeed::Int64, evolSeed::Int64, genNum::Int64, envMoveTo::Int64, envMoveFrom::Int64, dataDir::String)
    # Pull in population.jld file to get the gene expression results
    folderName = joinpath(dataDir,string("Pop_Seed_",popSeed,"\\Pop_Seed_",popSeed,"EvolSeed_",evolSeed,"\\Pop_Seed_",popSeed,"EvolSeed_",evolSeed,"_Gen",genNum))
    popFile = joinpath(folderName, string("Pop_Seed_",popSeed,"EvolSeed_",evolSeed,"_Gen",genNum,".jld"))
    popRead = jldopen(popFile, "r") do file
        read(file, "Population")
    end

    # Initialize Matrix to Store Results to use for PCA:
    matForPCA = zeros(Float64, G, N*4)
    namesVec =  fill("", N*4)
    # Pull out resulting stable gene expression pattern for every individual
    # for each of the 4 cases below to build a matrix with gene expression
    # results for these for categories for every individual, where rows are individuals:
    # 1) Immediately after evolution in Env 1
    matForPCA[:,1:N] = hcat(map(x -> popRead.individuals[x].develstate[:,envMoveFrom], 1:N)...)
    namesVec[1:N] = fill(string("Env",envMoveFrom), N)
    # 2) Immediately after evolution in Env2
    matForPCA[:,N+1:N*2] = hcat(map(x -> popRead.individuals[x].develstate[:,envMoveTo], 1:N)...)
    namesVec[N+1:N*2] = fill(string("Env",envMoveTo), N)
    for i = 1:N
        indWithBrokenPolycombReplaceEnvs, indWithIntactPolycombReplaceEnvs = disruptPolycomb(popRead.individuals[i], envMoveTo, envMoveFrom, popRead.founder, [1, 2])
        # 3) When Break PRC1and2 and move from Env 1 --> Env 2 (or vice versus)
        matForPCA[:,(N*2+1)+(i-1)] = indWithBrokenPolycombReplaceEnvs.develstate[:,envMoveTo]
        # 4) When PRC1and2 remain intact and move from Env 1 --> Env 2 (or vice versus)
        matForPCA[:,(N*3+1)+(i-1)] = indWithIntactPolycombReplaceEnvs.develstate[:,envMoveTo]
    end
    namesVec[(N*2+1):N*3] = fill(string("BrokenEnv",envMoveFrom,"To",envMoveTo), N)
    namesVec[(N*3+1):N*4] = fill(string("IntactEnv",envMoveFrom,"To",envMoveTo), N)

    M = fit(PCA, matForPCA; pratio=1, maxoutdim=2)
    components1 = MultivariateStats.transform(M, matForPCA)
    components1 = transpose(components1)
    df = DataFrame(components1, [:PC1, :PC2])
    df[!, :case] = namesVec
    scatter1 = @df df scatter(:PC1, :PC2, group=:case, legend=:topright, markersize=[3 3 3 3], c=[:firebrick2 :lime :blue :cyan], shape=[:o :o :x :o], markerstrokewidth=[0.1 0.1 0 0.1], title = string("PCA for Pop ", popSeed, " Evol Traj ", evolSeed, " at Gen ", genNum, " for Env ", envMoveFrom," To ",envMoveTo), xlabel = "PC1", ylabel = "PC2")
    savefig(scatter1, joinpath(folderName,string("PcaPop", popSeed, "EvolTraj", evolSeed, "Gen", genNum, "Env", envMoveFrom,"To",envMoveTo,".pdf")))
    display(scatter1)
end

# UMAP plot for a single population:
function plotUmapForSinglePopAndEvolTraj(popSeed::Int64, evolSeed::Int64, genNum::Int64, envMoveTo::Int64, envMoveFrom::Int64, dataDir::String)
    # Pull in population.jld file to get the gene expression results
    folderName = joinpath(dataDir,string("Pop_Seed_",popSeed,"\\Pop_Seed_",popSeed,"EvolSeed_",evolSeed,"\\Pop_Seed_",popSeed,"EvolSeed_",evolSeed,"_Gen",genNum))
    popFile = joinpath(folderName, string("Pop_Seed_",popSeed,"EvolSeed_",evolSeed,"_Gen",genNum,".jld"))
    popRead = jldopen(popFile, "r") do file
        read(file, "Population")
    end

    # Initialize Matrix to Store Results to use for PCA:
    matForPCA = zeros(Float64, G, N*4)
    namesVec =  fill("", N*4)
    # Pull out resulting stable gene expression pattern for every individual
    # for each of the 4 cases below to build a matrix with gene expression
    # results for these for categories for every individual, where rows are individuals:
    # 1) Immediately after evolution in Env 1
    matForPCA[:,1:N] = hcat(map(x -> popRead.individuals[x].develstate[:,envMoveFrom], 1:N)...)
    namesVec[1:N] = fill(string("Env",envMoveFrom), N)
    # 2) Immediately after evolution in Env2
    matForPCA[:,N+1:N*2] = hcat(map(x -> popRead.individuals[x].develstate[:,envMoveTo], 1:N)...)
    namesVec[N+1:N*2] = fill(string("Env",envMoveTo), N)
    for i = 1:N
        indWithBrokenPolycombReplaceEnvs, indWithIntactPolycombReplaceEnvs = disruptPolycomb(popRead.individuals[i], envMoveTo, envMoveFrom, popRead.founder, [1, 2])
        # 3) When Break PRC1and2 and move from Env 1 --> Env 2 (or vice versus)
        matForPCA[:,(N*2+1)+(i-1)] = indWithBrokenPolycombReplaceEnvs.develstate[:,envMoveTo]
        # 4) When PRC1and2 remain intact and move from Env 1 --> Env 2 (or vice versus)
        matForPCA[:,(N*3+1)+(i-1)] = indWithIntactPolycombReplaceEnvs.develstate[:,envMoveTo]
    end
    namesVec[(N*2+1):N*3] = fill(string("BrokenEnv",envMoveFrom,"To",envMoveTo), N)
    namesVec[(N*3+1):N*4] = fill(string("IntactEnv",envMoveFrom,"To",envMoveTo), N)

    res_umap = umap(matForPCA)

    df = DataFrame(transpose(res_umap), [:UMAP1, :UMAP2])
    df[!, :case] = namesVec
    scatter1 = @df df scatter(:UMAP1, :UMAP2, group=:case, legend=:topright, markersize=[3 3 3 3], c=[:firebrick2 :lime :blue :cyan], shape=[:o :o :x :o], markerstrokewidth=[0.1 0.1 0 0.1], title = string("UMAP for Pop ", popSeed, " Evol Traj ", evolSeed, " at Gen ", genNum, " for Env ", envMoveFrom," To ",envMoveTo), xlabel = "UMAP Dim 1", ylabel = "UMAP Dim 2")
    savefig(scatter1, joinpath(folderName,string("UmapPop", popSeed, "EvolTraj", evolSeed, "Gen", genNum, "Env", envMoveFrom,"To",envMoveTo,".pdf")))
    display(scatter1)
end


#-------------------------------------------------------------------------------
# PCA plot for supplemental material when transfer from env 1 to 2 and then
# transfer back to env 1 to make sure switching is sustainable:
# PCA plot for a single population:
function suppSwitchBackPcaForSinglePopAndEvolTraj(popSeed::Int64, evolSeed::Int64, genNum::Int64, envMoveTo::Int64, envMoveFrom::Int64, dataDir::String)
    # Pull in population.jld file to get the gene expression results
    folderName = joinpath(dataDir,string("Pop_Seed_",popSeed,"\\Pop_Seed_",popSeed,"EvolSeed_",evolSeed,"\\Pop_Seed_",popSeed,"EvolSeed_",evolSeed,"_Gen",genNum))
    popFile = joinpath(folderName, string("Pop_Seed_",popSeed,"EvolSeed_",evolSeed,"_Gen",genNum,".jld"))
    popRead = jldopen(popFile, "r") do file
        read(file, "Population")
    end

    # Initialize Matrix to Store Results to use for PCA:
    matForPCA = zeros(Float64, G, N*5)
    namesVec =  fill("", N*5)
    # Pull out resulting stable gene expression pattern for every individual
    # for each of the 4 cases below to build a matrix with gene expression
    # results for these for categories for every individual, where rows are individuals:
    # 1) Immediately after evolution in Env 1
    matForPCA[:,1:N] = hcat(map(x -> popRead.individuals[x].develstate[:,envMoveFrom], 1:N)...)
    namesVec[1:N] = fill(string("Env",envMoveFrom), N)
    # 2) Immediately after evolution in Env2
    matForPCA[:,N+1:N*2] = hcat(map(x -> popRead.individuals[x].develstate[:,envMoveTo], 1:N)...)
    namesVec[N+1:N*2] = fill(string("Env",envMoveTo), N)
    for i = 1:N
        # Break Polycomb-like mechanisms and switch environment
        indWithBrokenPolycombReplaceEnvs, indWithIntactPolycombReplaceEnvs = disruptPolycomb(popRead.individuals[i], envMoveTo, envMoveFrom, popRead.founder, [1, 2])
        # 3) When Break PRC1and2 and move from Env 1 --> Env 2 (or vice versus)
        matForPCA[:,(N*2+1)+(i-1)] = indWithBrokenPolycombReplaceEnvs.develstate[:,envMoveTo]
        # # 4) When PRC1and2 remain intact and move from Env 1 --> Env 2 (or vice versus)
        # matForPCA[:,(N*3+1)+(i-1)] = indWithIntactPolycombReplaceEnvs.develstate[:,envMoveTo]

        #
        indWithBrokenPcgMoveBackToOrgEnv = deepcopy(indWithBrokenPolycombReplaceEnvs)
        indWithIntactPcgMoveBackToOrgEnv = deepcopy(indWithIntactPolycombReplaceEnvs)
        iterateindForStability(indWithBrokenPcgMoveBackToOrgEnv, envMoveFrom, popRead.founder.EnvState[:,envMoveFrom], envMoveTo, indWithBrokenPcgMoveBackToOrgEnv) # Intact
        iterateindForStability(indWithIntactPcgMoveBackToOrgEnv, envMoveFrom, popRead.founder.EnvState[:,envMoveFrom], envMoveTo, indWithIntactPcgMoveBackToOrgEnv)
        # 5) When Break PRC1and2 and move from Env 1 --> Env 2 --> Env 1 (or vice versus)
        matForPCA[:,(N*3+1)+(i-1)] = indWithBrokenPcgMoveBackToOrgEnv.develstate[:,envMoveFrom]
        # 6) When PRC1and2 remain intact and move from Env 1 --> Env 2 --> Env 1 (or vice versus)
        # matForPCA[:,(N*5+1)+(i-1)] = indWithIntactPcgMoveBackToOrgEnv.develstate[:,envMoveFrom]

        # 7) When Break PRC1and2 but stay in EnvMoveFrom (so do not move environment)
        indBrokenPcgStayInSameEnv, indIntactPcgStayInSameEnv = disruptPolycomb(popRead.individuals[i], envMoveFrom, envMoveFrom, popRead.founder, [1, 2])
        matForPCA[:,(N*4+1)+(i-1)] = indBrokenPcgStayInSameEnv.develstate[:,envMoveFrom]
        # matForPCA[:,(N*7+1)+(i-1)] = indIntactPcgStayInSameEnv.develstate[:,envMoveFrom]
    end
    namesVec[(N*2+1):N*3] = fill(string("BrokenEnv",envMoveFrom,"To",envMoveTo), N)
    # namesVec[(N*3+1):N*4] = fill(string("IntactEnv",envMoveFrom,"To",envMoveTo), N)
    namesVec[(N*3+1):N*4] = fill(string("BrokenEnvTo",envMoveFrom,"To",envMoveTo,"To",envMoveFrom), N)
    # namesVec[(N*5+1):N*6] = fill(string("IntactEnvBackTo",envMoveFrom), N)
    namesVec[(N*4+1):N*5] = fill(string("BrokenEnv",envMoveFrom), N)
    # namesVec[(N*7+1):N*8] = fill(string("IntactStayEnv",envMoveFrom), N)

    M = fit(PCA, matForPCA; pratio=1, maxoutdim=2)
    components1 = MultivariateStats.transform(M, matForPCA)
    components1 = transpose(components1)
    df = DataFrame(components1, [:PC1, :PC2])
    df[!, :case] = namesVec
    scatter1 = @df df scatter(:PC1, :PC2, group=:case, legend=:bottomleft, markersize=[3 3 3 3 3], c=[:tan1 :firebrick2 :purple2 :lime :blue], shape=[:o :o :x :o :x], markerstrokewidth=[0.1 0.1 0 0.1 0], title = string("PCA for Pop ", popSeed, " Evol Traj ", evolSeed, " at Gen ", genNum, " for Env ", envMoveFrom," To ",envMoveTo, "and Back To", envMoveFrom), xlabel = "PC1", ylabel = "PC2")
    savefig(scatter1, joinpath(folderName,string("SuppPcaPop", popSeed, "EvolTraj", evolSeed, "Gen", genNum, "Env", envMoveFrom,"To",envMoveTo,"BackTo",envMoveFrom,".pdf")))
    display(scatter1)
end


# PCA plot for supplemental material when don't switch environments (so "switch" to environment currently in)
function suppStayInEnvPcaForSinglePopAndEvolTraj(popSeed::Int64, evolSeed::Int64, genNum::Int64, envMoveTo::Int64, envMoveFrom::Int64, dataDir::String)
    # Pull in population.jld file to get the gene expression results
    folderName = joinpath(dataDir,string("Pop_Seed_",popSeed,"\\Pop_Seed_",popSeed,"EvolSeed_",evolSeed,"\\Pop_Seed_",popSeed,"EvolSeed_",evolSeed,"_Gen",genNum))
    popFile = joinpath(folderName, string("Pop_Seed_",popSeed,"EvolSeed_",evolSeed,"_Gen",genNum,".jld"))
    popRead = jldopen(popFile, "r") do file
        read(file, "Population")
    end

    # Initialize Matrix to Store Results to use for PCA:
    matForPCA = zeros(Float64, G, N*4)
    namesVec =  fill("", N*4)
    # Pull out resulting stable gene expression pattern for every individual
    # for each of the 4 cases below to build a matrix with gene expression
    # results for these for categories for every individual, where rows are individuals:
    # 1) Immediately after evolution in Env 1
    matForPCA[:,1:N] = hcat(map(x -> popRead.individuals[x].develstate[:,envMoveFrom], 1:N)...)
    namesVec[1:N] = fill(string("Env",envMoveFrom), N)
    # 2) Immediately after evolution in Env2
    matForPCA[:,N+1:N*2] = hcat(map(x -> popRead.individuals[x].develstate[:,envMoveTo], 1:N)...)
    namesVec[N+1:N*2] = fill(string("Env",envMoveTo), N)
    for i = 1:N
        # 7) When Break PRC1and2 but stay in EnvMoveFrom (so do not move environment)
        indBrokenPcgStayInSameEnv, indIntactPcgStayInSameEnv = disruptPolycomb(popRead.individuals[i], envMoveFrom, envMoveFrom, popRead.founder, [1, 2])
        matForPCA[:,(N*2+1)+(i-1)] = indBrokenPcgStayInSameEnv.develstate[:,envMoveFrom]
        matForPCA[:,(N*3+1)+(i-1)] = indIntactPcgStayInSameEnv.develstate[:,envMoveFrom]
    end
    namesVec[(N*2+1):N*3] = fill(string("BrokenStayEnv",envMoveFrom), N)
    namesVec[(N*3+1):N*4] = fill(string("IntactStayEnv",envMoveFrom), N)

    M = fit(PCA, matForPCA; pratio=1, maxoutdim=2)
    components1 = MultivariateStats.transform(M, matForPCA)
    components1 = transpose(components1)
    df = DataFrame(components1, [:PC1, :PC2])
    df[!, :case] = namesVec
    scatter1 = @df df scatter(:PC1, :PC2, group=:case, legend=:topright, markersize=[3 3 3 3], c=[:lime :blue :tan1 :purple4], shape=[:o :x :o :o], markerstrokewidth=[0.1 0.1 0.1 0.1], title = string("PCA for Pop ", popSeed, " Evol Traj ", evolSeed, " at Gen ", genNum, " for Env ", envMoveFrom," To ",envMoveTo, "and Back To", envMoveFrom), xlabel = "PC1", ylabel = "PC2")
    savefig(scatter1, joinpath(folderName,string("SuppPcaPop", popSeed, "EvolTraj", evolSeed, "Gen", genNum, "StayInEnv", envMoveFrom,".pdf")))
    display(scatter1)
end
