function reservestables(model::ReservesModel, simulated::ModelSimulation, moments::ModelMoments)

detnum=model.compuparams.debtnum
resnum=model.compuparams.resnum
ynum=model.compuparams.ynum

# 1. Set directory
currdir=pwd()
# Inside ReservesMod
cd("..\\..\\tables")
# Latex directory (uncomment to put output on JMP directly) 
# TEST BEFORE
# cd("..\\..\\..\\Latex\\Tables")

# 2. momentsparams table

momparamtab=open("tab_momentsparams.tex","w+")
# Initialize table
println(momparamtab, "\\begin{tabular}{lll}")
# Header
println(momparamtab, "	\\toprule")
println(momparamtab, "	Parameter & Value & Source \\\\")
println(momparamtab, "	\\midrule")
# Data
# risk free rate
println(momparamtab, "	\$r^f\$ & $(model.econparams.rfree) &  \\\\" )
# bond maturity
println(momparamtab, "	\$\\lambda\$ & $(model.econparams.llambda) &  \\\\" )
# bond coupon
println(momparamtab, "	\$z\$ & $(model.econparams.coupon) &  EMBI avg spread \\\\" )
# Output Autocorrelation
println(momparamtab, "	\$\\rho\$ & $(model.econparams.logoutputrho) &  \\\\" )
# Output Innovation Variance
println(momparamtab, "	\$\\sigma_{\\nu}\$ & $(model.econparams.logoutputrho) &  \\\\" )

# Lower rule end table and close
println(momparamtab, "	\\bottomrule")
println(momparamtab, "\\end{tabular}")
close(momparamtab)
cd(currdir)

end # Function end
