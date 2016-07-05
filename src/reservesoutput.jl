module reservesoutput
######################	
using ReservesTypes
using JLD
using Gadfly
using DataFrames

######################
export basemodel, basesimul, basemoments, reservesfigures
# 0. Load solved model
# Loads basemodel
@load "firstsolved.jld"
# 1. Make figures
include("reservesfigures.jl")
reservesfigures(basemodel, basesimul, basemoments)

#########################
end # module end