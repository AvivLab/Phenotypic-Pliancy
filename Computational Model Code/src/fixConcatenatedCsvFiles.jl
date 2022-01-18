# Code to crrect CSV dataframe mistake that swaps StdEnv1 and AvgEnv2 for the following measurement results:
# "df_btwCentralityTargetsRepressed.csv","df_connectivityOutgoingTargetsRepressed.csv",
# "df_connectivityIngoingTargetsRepressed.csv", "df_effectiveConnectivityIngoing.csv",
# "df_effectiveConnectivityOutgoing.csv", "df_effectiveBtwCentrality.csv", and "df_numTargetsRepressed.csv"

#---------------------------------------------------------------------------
# Load packagees and files needed: (Note: install by Pkg.add("Distributions") in julia)
using CSV
using DataFrames
using StatsPlots

include("dataAnalysisFunctions.jl")

#--------------------------------------------------------------------------------
#dataOutputFolderName = joinpath("..", "dataOutputForPaper", "Summary1_10000")
dataOutputFolderName = joinpath("..", "dataOutputForPaper")

generateDfForSinglePopAvgEvolTraj("df_btwCentralityTargetsRepressed.csv", dataOutputFolderName)
generateDfForSinglePopAvgEvolTraj("df_connectivityOutgoingTargetsRepressed.csv", dataOutputFolderName)
generateDfForSinglePopAvgEvolTraj("df_connectivityIngoingTargetsRepressed.csv", dataOutputFolderName)
generateDfForSinglePopAvgEvolTraj("df_effectiveConnectivityIngoing.csv", dataOutputFolderName)
generateDfForSinglePopAvgEvolTraj("df_effectiveConnectivityOutgoing.csv", dataOutputFolderName)
generateDfForSinglePopAvgEvolTraj("df_effectiveBtwCentrality.csv", dataOutputFolderName)
generateDfForSinglePopAvgEvolTraj("df_numTargetsRepressed.csv", dataOutputFolderName)
