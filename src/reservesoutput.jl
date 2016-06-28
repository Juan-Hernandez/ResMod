module reservesoutput
######################	
using ReservesTypes
using JLD
using Gadfly
using DataFrames

######################
export basemodel, basesimul, basemoments
# 0. Load solved model
# Loads basemodel
@load "simulatedmodel.jld"
# 1. Make figures
include("reservesfigures.jl")
reservesfigures(basemodel, basesimul, basemoments)

#########################
end # module end