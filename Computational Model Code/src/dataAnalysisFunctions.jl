# Data Analysis Functions written by Maryl Lambros November 2021
#-------------------------------------------------------------------------------
# Function to correct naming mistake of StdEnv1 and AvgEnv2 (should be swapped)
function generateDfForSinglePopAvgEvolTraj(measurementFileName::String, dataOutputFolderName::String)
    measureDf = DataFrame(CSV.File(joinpath(dataOutputFolderName, measurementFileName),
        types=Dict(:Pop=>Int64, :Gen=>Int64)))
    CSV.write(joinpath(dataOutputFolderName,measurementFileName), rename!(measureDf, ["Pop", "Evol", "Gen", "AvgEnv1", "AvgEnv2", "StdEnv1", "StdEnv2"]))
end

#-------------------------------------------------------------------------------
# Plot functions:
function plotMeasurements(measurementFileName::String, dataOutputFolderName::String, popSeed::Int64, evolSeed::Int64)
    if popSeed == 0 # if want to plot average of all populations and evolutionary trajectories
        popTitle = " for all Pops and Evol Traj"
        df = generateDfForAllPopAndEvolTraj(measurementFileName, dataOutputFolderName)
    elseif popSeed > 0 && evolSeed == 0 # if want to plot single population but average of all evolutionary trajectories
        popTitle = string(" for All Evol Traj for Pop ", popSeed)
        df = generateDfForSinglePopAvgEvolTraj(measurementFileName, dataOutputFolderName, popSeed)
    elseif popSeed > 0 && evolSeed > 0 # if want to plot single population and evolutionary trajectories separately
        popTitle = string(" for Evol Traj ", evolSeed, " for Pop ", popSeed)
        df = generateDfForSinglePopSingleEvolTraj(measurementFileName, dataOutputFolderName, popSeed, evolSeed)
    end

    numMeasures = size(df)[2] - 1
    if numMeasures == 6 && !(occursin("df_unstable",measurementFileName))
        if occursin("founder",measurementFileName)
            plotWhenFounderRepressedGenes(measurementFileName, dataOutputFolderName, df, popTitle, "Generations")
        else
            plotWhen6MeasurementResults(measurementFileName, dataOutputFolderName, df, popTitle, "Generations")
        end
    elseif numMeasures == 4
        plotWhen4MeasurementResults(measurementFileName, dataOutputFolderName, df, popTitle, "Generations")
    elseif occursin("df_unstable",measurementFileName) #numMeasures is actually 6 but only 3 measurements taken but generate both mena and std
        plotWhen3MeasurementResults(measurementFileName, dataOutputFolderName, df, popTitle, "Generations")
    elseif numMeasures == 2
        plotWhen2MeasurementResults(measurementFileName, dataOutputFolderName, df, popTitle, "Generations")
    end
end

function plotEuclideanDistResultsTogether(measurementFileName1::String, measurementFileName2::String, measurementFileName3::String, measurementFileName4::String, dataOutputFolderName::String)
    popTitle = " for all Pops and Evol Traj"
    df1 = generateDfForAllPopAndEvolTraj(measurementFileName1, dataOutputFolderName)
    df2 = generateDfForAllPopAndEvolTraj(measurementFileName2, dataOutputFolderName)
    df3 = generateDfForAllPopAndEvolTraj(measurementFileName3, dataOutputFolderName)
    df4 = generateDfForAllPopAndEvolTraj(measurementFileName4, dataOutputFolderName)

end

function plotWhen6MeasurementResults(measurementFileName::String, dataOutputFolderName::String, df::DataFrame, popTitle::String, xLab::String)
    # Plotting average of all populations data generated and collected when measurements are:
    # Measurements: 1) "AvgPRC1", 2) "StdPRC1", 3) "AvgPRC2", 4) "StdPRC2", 5) "AvgPRC1and2", 6) "StdPRC1and2"
    # for the following measurement types:

    ## Functional pliancy:
    # 1) df_phenotypicPliancyScore1to2, and 2) df_phenotypicPliancyScore2to1

    ## Gene euclidean distances:
    # 1) df_difIntact1_BrokenEnv1ToEnv2, 2) df_difIntact2_BrokenEnv1ToEnv2
    # 3) df_difIntact1_IntactEnv1toEnv2, 4) df_difIntact2_IntactEnv1toEnv2,
    # 5) df_difBrokenEnv1ToEnv2_IntactEnv1ToEnv2, 6) df_difIntactEnv1_IntactEnv2
    # 7) df_difIntact1_BrokenEnv2ToEnv1, 8) df_difIntact2_BrokenEnv2ToEnv1
    # 9) df_difIntact1_IntactEnv2toEnv1, 10) df_difIntact2_IntactEnv2toEnv1
    # and 11) df_difBrokenEnv2ToEnv1_IntactEnv2ToEnv1
    prcVec = ["PRC1" "PRC2" "PRC1and2"]
    if (occursin("df_dif",measurementFileName)) # if euclidean distance measurement
        meas1, meas2 = split(match(r"(?<=df_dif)(\w+)", measurementFileName).match, "_")
        plotTitle = string("Euclidean Dist Btw ",meas1, " and ", meas2, popTitle)
        ylab = string("Avg. Euclidean Distance Btw ", meas1, " and ", meas2)
    else # functional pliancy measurement
        meas1, meas2 = parse.(Int64, [m.match for m in eachmatch(r"\d", measurementFileName)])
        plotTitle = string("Functional Pliancy when move from Env ", meas1," to Env ", meas2, popTitle)
        ylab = string("Avg. Func Pliancy Score for Env ", meas1, " to Env ", meas2)
    end

    plot1 = @df df plot(cols(1), [:AvgPRC1_mean :AvgPRC2_mean :AvgPRC1and2_mean]
         , yerror=[:StdPRC1_mean./sqrt(N) :StdPRC2_mean./sqrt(N) :StdPRC1and2_mean./sqrt(N)]
         , linewidth=2, linestyle=:solid, legend=:bottomright, c=[:blue :cyan :firebrick2]
         , label=prcVec, title=plotTitle)
    xlabel!(xLab)
    ylabel!(ylab)
    savefig(plot1, joinpath(dataOutputFolderName, string(plotTitle, ".pdf")))
    display(plot1)
end

function plotWhenFounderRepressedGenes(measurementFileName::String, dataOutputFolderName::String, df::DataFrame, popTitle::String, xLab::String)

    plotTitle = string("Repressed genes", popTitle)
    ylab = string("% Avg Genes Repressed by PRCs & Founder")

    plot1 = @df df plot(cols(1), [:AvgEnv1_mean :AvgEnv2Dev_mean :AvgEnv2Opt_mean]
         , yerror=[:StdEnv1_mean./sqrt(N) :StdEnv2Dev_mean./sqrt(N) :StdEnv2Opt_mean./sqrt(N)]
         , linewidth=2, linestyle=:solid, label=["Env1" "Env2Dev" "Env2Opt"], legend=:bottomright, title=plotTitle)
    xlabel!(xLab)
    ylabel!(ylab)
    savefig(plot1, joinpath(dataOutputFolderName, string(plotTitle, ".pdf")))
    display(plot1)
end

function plotWhen4MeasurementResults(measurementFileName::String, dataOutputFolderName::String, df::DataFrame, popTitle::String, xLab::String)
    # df_effectiveConnectivityIngoing, df_effectiveConnectivityOutgoing, df_effectiveBtwCentrality
    # df_connectivityIngoingTargetsRepressed, df_connectivityOutgoingTargetsRepressed
    # df_btwCentralityTargetsRepressed, df_numTargetsRepressed
    meas1 = match(r"(?<=df_)(\w+)", measurementFileName).match
    plotTitle = string(meas1, " in Env 1 and Env 2", popTitle)
    ylab = string("Avg. ", meas1)
    plot1 = @df df plot(cols(1), [:AvgEnv1_mean :AvgEnv2_mean]
         , yerror=[:StdEnv1_mean./sqrt(N) :StdEnv2_mean./sqrt(N)]
         , linewidth=2, linestyle=:solid, legend=:topleft, c=[:firebrick2 :blue], label=["Env1" "Env2"], title=plotTitle)
    xlabel!(xLab)
    ylabel!(ylab)
    savefig(plot1, joinpath(dataOutputFolderName, string(plotTitle, ".pdf")))
    display(plot1)
end

function plotWhen3MeasurementResults(measurementFileName::String, dataOutputFolderName::String, df::DataFrame, popTitle::String, xLab::String)
    # df_unstableDisrupted1to2, df_unstableDisrupted2to1

    meas1 = match(r"(?<=df_)(\w+)", measurementFileName).match
    plotTitle = string(meas1, " when PRC1, PRC2, and PRC1 and 2 broken", popTitle)
    ylab = string("Avg. ", meas1)

    plot1 = @df df plot(cols(1), [:PRC1_mean :PRC2_mean :PRC1and2_mean]
         , yerror=[:PRC1_std./sqrt(N) :PRC2_std./sqrt(N) :PRC1and2_std./sqrt(N)]
         , linewidth=2, linestyle=:solid, legend=:topright, c=[:blue :cyan :firebrick2]
         , label=["PRC1" "PRC2" "PRC1and2"], title=plotTitle)
    xlabel!(xLab)
    ylabel!(ylab)
    savefig(plot1, joinpath(dataOutputFolderName, string(plotTitle, ".pdf")))
    display(plot1)
end


function plotWhen2MeasurementResults(measurementFileName::String, dataOutputFolderName::String, df::DataFrame, popTitle::String, xLab::String)
    # df_connectivityIngoingTargets, df_connectivityOutgoingTargets, df_btwCentralityTargets
    # df_numPcgTargets, df_totalFitness, df_fitnessEnv1, df_fitnessEnv2

    meas1 = match(r"(?<=df_)(\w+)", measurementFileName).match
    plotTitle = string(meas1, popTitle)
    ylab = string("Avg. ", meas1)

    plot1 = @df df plot(cols(1), :Avg_mean
         , yerror=:Std_mean./sqrt(N)
         , linewidth=2, linestyle=:solid, label=["Avg"], title=plotTitle)
    xlabel!(xLab)
    ylabel!(ylab)
    savefig(plot1, joinpath(dataOutputFolderName, string(plotTitle, ".pdf")))
    display(plot1)
end


function plotPopLevelNetworkArchMeasures(measurementFileName::String, dataOutputFolderName::String)
    # df_connectivityIngoingTargets, df_connectivityOutgoingTargets, df_btwCentralityTargets
    # df_numPcgTargets, df_totalFitness, df_fitnessEnv1, df_fitnessEnv2
    df = DataFrame(CSV.File(joinpath(dataOutputFolderName, measurementFileName),
        types=Dict(:Pop=>Int64)))

    meas1 = match(r"(?<=df_)(\w+)", measurementFileName).match
    plotTitle = string("Avg. ", meas1, " for each population")
    ylab = string("Avg. ", meas1)

    plot1 = @df df scatter(:Pop, cols(2), title=plotTitle)
    xlabel!("Population")
    ylabel!(ylab)
    savefig(plot1, joinpath(dataOutputFolderName, string(plotTitle, ".pdf")))
    display(plot1)
end

# Function to plot multiple different dataframe results in same plot:
function plotMultipleDfsInOne(measurementFileNames::Vector{String}, dataOutputFolderName::String, prc::String)
    # Plotting average of all populations data generated and collected when measurements are:
    # Measurements: 1) "AvgPRC1", 2) "StdPRC1", 3) "AvgPRC2", 4) "StdPRC2", 5) "AvgPRC1and2", 6) "StdPRC1and2"
    # for the following measurement types when want to plot multiple together from the following:

    ## Gene euclidean distances:
    # 1) df_difIntact1_BrokenEnv1ToEnv2, 2) df_difIntact2_BrokenEnv1ToEnv2
    # 3) df_difIntact1_IntactEnv1toEnv2, 4) df_difIntact2_IntactEnv1toEnv2,
    # 5) df_difBrokenEnv1ToEnv2_IntactEnv1ToEnv2, 6) df_difIntactEnv1_IntactEnv2
    # 7) df_difIntact1_BrokenEnv2ToEnv1, 8) df_difIntact2_BrokenEnv2ToEnv1
    # 9) df_difIntact1_IntactEnv2toEnv1, 10) df_difIntact2_IntactEnv2toEnv1
    # and 11) df_difBrokenEnv2ToEnv1_IntactEnv2ToEnv1
    colorVec = [:orange, :blue, :purple, :red]

    df1 = generateDfForAllPopAndEvolTraj(measurementFileNames[1], dataOutputFolderName)
    meas1, meas2 = split(match(r"(?<=df_dif)(\w+)", measurementFileNames[1]).match, "_")
    plotTitle = string("Euclidean Dist for ", prc, " for all Pops and Evol Traj")
    ylab = string("Avg. Euclidean Distance")

    avgColToPlot = findall( x -> occursin(string("Avg",prc,"_mean"), x), names(df1))[1]
    stdColToPlot = findall( x -> occursin(string("Std",prc,"_mean"), x), names(df1))[1]
    plot1 = @df df1 plot(cols(1), cols(avgColToPlot)
         , yerror=cols(stdColToPlot)./sqrt(N)
         , linewidth=2, linestyle=:solid, label=string(meas1, "-", meas2), legend=:left, color = colorVec[1], title=plotTitle)
    xlabel!("Generations")
    ylabel!(ylab)

    for i = 2:length(measurementFileNames)
        df = generateDfForAllPopAndEvolTraj(measurementFileNames[i], dataOutputFolderName)
        meas1, meas2 = split(match(r"(?<=df_dif)(\w+)", measurementFileNames[i]).match, "_")
        plot1 = @df df plot!(plot1, cols(1), cols(avgColToPlot)
             , yerror=cols(stdColToPlot)./sqrt(N)
             , linewidth=2, linestyle=:solid, label=string(meas1, "-", meas2), color = colorVec[i])
    end
    savefig(plot1, joinpath(dataOutputFolderName, string(plotTitle, ".pdf")))
    display(plot1)
end


# Function to plot 2 different dataframe results in same plot for connectivity
# and betweenness centrality:
function plotMultipleNetworkArchsInOne(measurementFileNames::Vector{String}, dataOutputFolderName::String, measurementString::String)
    #
    df1 = generateDfForAllPopAndEvolTraj(measurementFileNames[1], dataOutputFolderName)
    plotTitle = string(measurementString, " in Env 1 and Env 2 for all Pops & Evol Trajs")
    ylab = string("Avg. ", measurementString)
    plot1 = @df df1 plot(:Gen, [:AvgEnv1_mean :AvgEnv2_mean]
         , yerror=[:StdEnv1_mean./sqrt(N) :StdEnv2_mean./sqrt(N)]
         , linewidth=2, linestyle=[:dot :solid], c=[:cyan :magenta], label=["Env1 All Genes" "Env2 All Genes"], title=plotTitle, legend=:bottomleft)
    xlabel!("Generations")
    ylabel!(ylab)

    df2 = generateDfForAllPopAndEvolTraj(measurementFileNames[2], dataOutputFolderName)
    plot1 = @df df2 plot!(plot1, :Gen, [:AvgEnv1_mean :AvgEnv2_mean]
         , yerror=[:StdEnv1_mean./sqrt(N) :StdEnv2_mean./sqrt(N)]
         , linewidth=2, linestyle=[:dot :solid], c=[:blue :firebrick2], label=["Env1 Repressed" "Env2 Repressed"])

    df3 = generateDfForAllPopAndEvolTraj(measurementFileNames[3], dataOutputFolderName)
    plot1 = @df df3 plot!(plot1, :Gen, :Avg_mean
         , yerror=:Std_mean./sqrt(N)
         , linewidth=2, linestyle=:dashdot, c=:lime, label="Target Genes")

    savefig(plot1, joinpath(dataOutputFolderName, string(plotTitle, ".pdf")))
    display(plot1)
end

#-------------------------------------------------------------------------------
# Functions when looking at averages across ALL populations and evolutionary trajectories:
function generateDfForAllPopAndEvolTraj(measurementFileName::String, dataOutputFolderName::String)
    measureDf = DataFrame(CSV.File(joinpath(dataOutputFolderName, measurementFileName),
        types=Dict(:Pop=>Int64, :Gen=>Int64)))
    numMeasures = size(measureDf)[2] - 3
    if numMeasures == 6 # 6 measurements
        if occursin("founder",measurementFileName)
            # Combine the different evolutionary trajectory results and plot results when break PRC1, PRC2, and PRC1and2:
            dfPops = combine(groupby(measureDf, [:Gen]), :AvgEnv1 => mean, :StdEnv1 => mean, :AvgEnv2Dev => mean, :StdEnv2Dev => mean, :AvgEnv2Opt => mean, :StdEnv2Opt => mean)
        else
            # Combine the different evolutionary trajectory results and plot results when break PRC1, PRC2, and PRC1and2:
            dfPops = combine(groupby(measureDf, [:Gen]), :AvgPRC1 => mean, :StdPRC1 => mean, :AvgPRC2 => mean, :StdPRC2 => mean, :AvgPRC1and2 => mean, :StdPRC1and2 => mean)
        end
    elseif numMeasures == 4
        measureDf.AvgEnv1 = map(x -> isnan(x) ? missing : x, measureDf.AvgEnv1)
        measureDf.StdEnv1 = map(x -> isnan(x) ? missing : x, measureDf.StdEnv1)
        measureDf.AvgEnv2 = map(x -> isnan(x) ? missing : x, measureDf.AvgEnv2)
        measureDf.StdEnv2 = map(x -> isnan(x) ? missing : x, measureDf.StdEnv2)
        dfPops = combine(groupby(measureDf, [:Gen]), :AvgEnv1 => (x -> mean(skipmissing(x))), :StdEnv1 => (x -> mean(skipmissing(x))), :AvgEnv2 => (x -> mean(skipmissing(x))), :StdEnv2 => (x -> mean(skipmissing(x))))
        dfPops = rename!(dfPops, ["Gen", "AvgEnv1_mean", "StdEnv1_mean", "AvgEnv2_mean", "StdEnv2_mean"])
    elseif numMeasures == 3
        dfPops = combine(groupby(measureDf, [:Gen]), :PRC1 => mean, :PRC2 => mean, :PRC1and2 => mean, :PRC1 => std, :PRC2 => std, :PRC1and2 => std)
    elseif numMeasures == 2
        dfPops = combine(groupby(measureDf, [:Gen]), :Avg => mean, :Std => mean)
    end
    return dfPops
end

# End Functions for looking at averages of ALL populations and evolutionary trajectories


#-------------------------------------------------------------------------------
# Functions for when want to look at single population with average of all 10 evolutionary trajectories:
function generateDfForSinglePopAvgEvolTraj(measurementFileName::String, dataOutputFolderName::String, popSeed::Int64)
    measureDf = DataFrame(CSV.File(joinpath(dataOutputFolderName, measurementFileName),
        types=Dict(:Pop=>Int64, :Gen=>Int64)))
    df = measureDf[(measureDf.Pop .== popSeed),:]
    numMeasures = size(df)[2] - 3
    if numMeasures == 6 # 6 measurements
        # Combine the different evolutionary trajectory results and plot results when break PRC1, PRC2, and PRC1and2:
        dfEvols = combine(groupby(df, :Gen), :AvgPRC1 => mean, :StdPRC1 => mean, :AvgPRC2 => mean, :StdPRC2 => mean, :AvgPRC1and2 => mean, :StdPRC1and2 => mean)
    elseif numMeasures == 4
        dfEvols = combine(groupby(df, :Gen), :AvgEnv1 => mean, :StdEnv1 => mean, :AvgEnv2 => mean, :StdEnv2 => mean)
    elseif numMeasures == 3
        dfEvols = combine(groupby(df, :Gen), :PRC1 => mean, :PRC2 => mean, :PRC1and2 => mean, :PRC1 => std, :PRC2 => std, :PRC1and2 => std)
    elseif numMeasures == 2
        dfEvols = combine(groupby(df, :Gen), :Avg => mean, :Std => mean)
    end
    return dfEvols
end

# End Functions for when want to look at single population with average of all 10 evolutionary trajectories:


#-------------------------------------------------------------------------------
# Functions for when want to plot single starting population and single evolutionary trajectory (or all 10 as own line?)
function generateDfForSinglePopSingleEvolTraj(measurementFileName::String, dataOutputFolderName::String, popSeed::Int64, evolSeed::Int64)
    measureDf = DataFrame(CSV.File(joinpath(dataOutputFolderName, measurementFileName),
        types=Dict(:Pop=>Int64, :Gen=>Int64)))
    df = measureDf[(measureDf.Pop .== popSeed) .& (measureDf.Evol .== evolSeed),:]
    numMeasures = size(df)[2] - 3
    if numMeasures == 6 # 6 measurements
        # Combine the different evolutionary trajectory results and plot results when break PRC1, PRC2, and PRC1and2:
        dfEvol = combine(groupby(df, :Gen), :AvgPRC1 => mean, :StdPRC1 => mean, :AvgPRC2 => mean, :StdPRC2 => mean, :AvgPRC1and2 => mean, :StdPRC1and2 => mean)
    elseif numMeasures == 4
        dfEvol = combine(groupby(df, :Gen), :AvgEnv1 => mean, :StdEnv1 => mean, :AvgEnv2 => mean, :StdEnv2 => mean)
    elseif numMeasures == 3
        dfEvol = combine(groupby(df, :Gen), :PRC1 => mean, :PRC2 => mean, :PRC1and2 => mean, :PRC1 => std, :PRC2 => std, :PRC1and2 => std)
    elseif numMeasures == 2
        dfEvol = combine(groupby(df, :Gen), :Avg => mean, :Std => mean)
    end
    return dfEvol
end

# End of Functions for when want to plot single starting population and single evolutionary trajectory (or all 10 as own line?)

#-------------------------------------------------------------------------------
function generateDfForAllPopsButPlotPopsSeparately(measurementFileName::String, dataOutputFolderName::String, popSeed::Int64, evolSeed::Int64)
    measureDf = DataFrame(CSV.File(joinpath(dataOutputFolderName, measurementFileName),
        types=Dict(:Pop=>Int64, :Gen=>Int64)))
    numMeasures = size(df)[2] - 3
    if numMeasures == 6 # 6 measurements
        # Combine the different evolutionary trajectory results and plot results when break PRC1, PRC2, and PRC1and2:
        dfEvol = combine(groupby(measureDf, [:Pop,:Gen]), :AvgPRC1 => mean, :StdPRC1 => mean, :AvgPRC2 => mean, :StdPRC2 => mean, :AvgPRC1and2 => mean, :StdPRC1and2 => mean)
    elseif numMeasures == 4
        dfEvol = combine(groupby(measureDf, [:Pop,:Gen]), :AvgEnv1 => mean, :StdEnv1 => mean, :AvgEnv2 => mean, :StdEnv2 => mean)
    elseif numMeasures == 3
        dfEvol = combine(groupby(measureDf, [:Pop,:Gen]), :PRC1 => mean, :PRC2 => mean, :PRC1and2 => mean, :PRC1 => std, :PRC2 => std, :PRC1and2 => std)
    elseif numMeasures == 2
        dfEvol = combine(groupby(measureDf, [:Pop,:Gen]), :Avg => mean, :Std => mean)
    end
    return dfEvol
end

#-------------------------------------------------------------------------------
# Plots for comparing across each population avg across all evolutionary
# trajectories for the final generation of evolution
function generateDfForFinalGeneration(measurementFileName::String, dataOutputFolderName::String, genNum::Int64)
    measureDf = DataFrame(CSV.File(joinpath(dataOutputFolderName, measurementFileName),
        types=Dict(:Pop=>Int64, :Gen=>Int64)))
    df = measureDf[(measureDf.Gen .== genNum),:]
    numMeasures = size(df)[2] - 3
    if numMeasures == 6 # 6 measurements
        # Combine the different evolutionary trajectory results and plot results when break PRC1, PRC2, and PRC1and2:
        dfEvols = combine(groupby(df, :Pop), :AvgPRC1 => mean, :StdPRC1 => mean, :AvgPRC2 => mean, :StdPRC2 => mean, :AvgPRC1and2 => mean, :StdPRC1and2 => mean)
    elseif numMeasures == 4
        dfEvols = combine(groupby(df, :Pop), :AvgEnv1 => mean, :StdEnv1 => mean, :AvgEnv2 => mean, :StdEnv2 => mean)
    elseif numMeasures == 3
        dfEvols = combine(groupby(df, :Pop), :PRC1 => mean, :PRC2 => mean, :PRC1and2 => mean, :PRC1 => std, :PRC2 => std, :PRC1and2 => std)
    elseif numMeasures == 2
        dfEvols = combine(groupby(df, :Pop), :Avg => mean, :Std => mean)
    end
    return dfEvols
end

function plotFinalGenMeasurements(measurementFileName::String, dataOutputFolderName::String, genNum::Int64)
    popTitle = string(" for Gen ", genNum ,"of all Pops")
    df = generateDfForFinalGeneration(measurementFileName, dataOutputFolderName, genNum)

    numMeasures = size(df)[2] - 1
    if numMeasures == 6 && !(occursin("df_unstable",measurementFileName))
        plotWhen6MeasurementResults(measurementFileName, dataOutputFolderName, df, popTitle, "Population")
    elseif numMeasures == 4
        plotWhen4MeasurementResults(measurementFileName, dataOutputFolderName, df, popTitle, "Population")
    elseif occursin("df_unstable",measurementFileName) #numMeasures is actually 6 but only 3 measurements taken but generate both mena and std
        plotWhen3MeasurementResults(measurementFileName, dataOutputFolderName, df, popTitle, "Population")
    elseif numMeasures == 2
        plotWhen2MeasurementResults(measurementFileName, dataOutputFolderName, df, popTitle, "Population")
    end
end
