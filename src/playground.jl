
#currdir=pwd()
#cd("..")
#push!(LOAD_PATH,pwd())
#cd(currdir)

# 0. Load packages and definitions
using ReservesTypes
VERSION<VersionNumber("0.7.0") ? using JLD : using JLD2
filename="playground.jld"	

# 1. Load parameters (in separate file for gitignore)
include("testparameters.jl")
# This defines basecompuparams, baseeconparams, basesolverparams

# 2. Create and initialize model
basemodel=ReservesModel(basecompuparams,baseeconparams)
modelinitialize!(basemodel)

# 3. Solve model
# 3.1 Solve routine and extract gaps 
solveroutvec=solvereservesmodel!(basemodel, basesolverparams)	
# # 3.1.1 Profiling
# basesolverparams.itermax=3
# Profile.@profile solveroutvec=solvereservesmodel!(basemodel, basesolverparams)	

# 3.2 Extract gaps
resiternum = floor(Int64, solveroutvec[1])
valuegap = solveroutvec[2]
pricegap = solveroutvec[3]
defaultgap = solveroutvec[4]
 
# 3.3 Save solved model
jldopen(filename,"w") do file
	write(file, "basemodel", basemodel)
	write(file, "basesolverparams", basesolverparams)
end	

# 4. Simulate and get moments
# 4.1. Simulate model
basesimul=ModelSimulation(100000)
simulatemodel!(basesimul,basemodel,true)
# 4.2 Save simulated
jldopen(filename,"r+") do file
	write(file, "basesimul", basesimul)
end
# 4.3. Obtain moments
basemoments=ModelMoments()
flag=getmoments!(basemoments, basesimul, basemodel.grids, 1000) # burnin 1000 observations
# 4.4 Print Moments
println("	debt	|	reserves	|	spravg		|	sprvar		|	defstat		|	defchoice	| sprXgrowth |   maxgap   |")
println("----------------------------------------------------------------------------------------------------------------------")
show(IOContext(stdout, :compact => true), basemoments.debtmean)
print("  | ")
show(IOContext(stdout, :compact => true), basemoments.reservesmean)
print("  | ")
show(IOContext(stdout, :compact => true), basemoments.spreadmean)
print(" | ")
show(IOContext(stdout, :compact => true), basemoments.spreadsigma)
print(" | ")
show(IOContext(stdout, :compact => true), basemoments.defaultstatemean)
print("  | ")
show(IOContext(stdout, :compact => true), basemoments.defaultchoicemean)
print("  | ")
show(IOContext(stdout, :compact => true), basemoments.spreadXgdp)
print("  | ")
show(IOContext(stdout, :compact => true), basemoments.spreadXgrowth)
print("  | ")
show(IOContext(stdout, :compact => true), basemoments.deltaspreadXgdp)
print("  | ")
show(IOContext(stdout, :compact => true), basemoments.deltaspreadXgrowth)
print("  | ")
show(IOContext(stdout, :compact => true), maximum([valuegap,pricegap,defaultgap]) )
println("  |")
println("======================================================================================================================")
# Profile.print()
# 4.5. Save moments
jldopen(filename,"r+") do file
	write(file, "basemoments", basemoments)
end	