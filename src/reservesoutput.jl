# module reservesoutput
######################	
using ReservesTypes
using JLD
using Gadfly
using DataFrames
using Colors
######################
export basemodel, basesimul, basemoments, reservesfigures
# 0. Load solved model
# Loads basemodel
workdir=pwd()
cd(homedir()"\\dropbox\\U-penn\\research\\ReservesProject\\Julia\\Results")
@load "basenewsolved.jld"
cd(workdir)
# 1. Make figures
include("reservesfigures.jl")
reservesfigures(basemodel, basesimul, basemoments)
# 2. Make tables
include("reservestables.jl")
reservestables(basemodel, basesimul, basemoments)
#########################
#end # module end