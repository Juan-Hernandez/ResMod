""" mainscript """
lamierda=pwd()
push!(LOAD_PATH, lamierda)
using ReservesTypes
using JLD
# 1. Load current model
@load "debugmodels.jld"
# 2. Continue solving
flag=solvereservesmodel!(basemodel,false,iternum,18000)
# 2.1 Save solved
jldopen("solvedmodel.jld", "w") do file
				write(file, "basemodel", basemodel)
				write(file, "flag", flag)
end
# 3. Simulate model
basesimul=ModelSimulation(100000)
simulatemodel!(basesimul,basemodel)
# 4. Obtain moments
basemoments=ModelMoments()
flag=getmoments!(basemoments, basesimul, basemodel.grids, 1000) # burnin 1000 observations
jldopen("simulatedmodel.jld", "w") do file
	write(file, "basemodel", basemodel)
	write(file, "basesimul", basesimul)
	write(file, "basemoments", basemoments)
end
