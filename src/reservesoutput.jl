# module reservesoutput
######################	using
using ReservesTypes
VERSION<VersionNumber("0.7.0") ? using JLD : using JLD2 #Needed for solvereservesmodel! to be able to save intermediae iterations
import Cairo # need to add
using Gadfly # need to add
using DataFrames # need to add
using Colors, FixedPointNumbers # need to add. Define colors
using LinearAlgebra, Statistics # For dot product and correlations
######################
export basemodel, basesimul, basemoments, reservesfigures
# 0. Load solved model
# Loads basemodel
workdir=pwd()
#cd(homedir()"\\dropbox\\U-penn\\research\\ReservesProject\\Julia\\Results")
# cd(homedir()"\\dropbox\\U-penn\\research\\ReservesProject\\Julia\\Results\\Drafts")
#@load "basenewsolved.jld"

cd(workdir)
# 1. Make figures
include("reservesfigures.jl")
outputframe=reservesfigures(basemodel, basesimul, basemoments)
# 2. Make tables
include("reservestables.jl")
reservestables(basemodel, basesimul, basemoments)
#########################
#end # module end