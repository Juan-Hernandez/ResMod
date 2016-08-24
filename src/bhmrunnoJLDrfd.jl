module bhmrunnoJLDbiggrid
#
addprocs(23)


push!(LOAD_PATH,pwd())
export basemodel

using ReservesTypes

curriter=0
sizegrid=140
entrygrid=1./collect(linspace(8.0, 12.0, 5))
freqgrid=collect(linspace(16.0, 40.0, 7))
durgrid=collect(linspace(5.0, 8.0, 4))


basecompuparams=ComputationParams(
	# Output Grid Lenght
	25,		# ynum::Int64
	# Debt grid parameters
	0.0,	# debtmin::Float64 
	0.16,	# debtmax::Float64
	21,		# debtnum::Int64
	# Reserves grid parameters
	0.0,		# resmin::Float64
	0.95, 	# resmax::Float64
	20,		# resnum::Int64
	# Temporary (smoothing shock parameters)
	0.011, 	# msigma::Float64
	2.0,	# msdwidth::Float64 (For convergence)
	13,		# mnum::Int64
	-100.0,	# thrmin::Float64
	)

# Temporary arrays to store last model solution. Helps with speed
tempvaluepay=Array{Float64}(basecompuparams.debtnum, basecompuparams.resnum, basecompuparams.mnum, basecompuparams.ynum, 2)
tempvaluedefault=Array{Float64}(basecompuparams.resnum, basecompuparams.ynum, 2)
tempbondprice=Array{Float64}(basecompuparams.debtnum, basecompuparams.resnum, basecompuparams.ynum, 2)



for i=1:sizegrid
	
	baseeconparams=EconParams(
		# Consumer
		0.9832,	# bbeta::Float64
		2,		# ggamma::Int64;  # HAS TO BE EQUAL TO 2. This cannot change. Will destroy threshold solution.
		# Risk free rate
		0.01, 	# rfree::Float64 # 4% yearly
		# Bond Maturity and coupon
		0.033, 	# llambda::Float64 # 5 year avg maturity
		1.0, 	# coupon:: Float64 
		# Expected Output grid parameters
		0.94, 	# logOutputRho::Float64
		0.045, 	# logOutputSigma::Float64
		# Default output cost parameters
		-0.6,	# defcost1::Float64
		0.751, 	# defcost2::Float64
		entrygrid[mod1(i,5)], # reentry::Float64
		# Sudden Stop Probability
		freqgrid[mod1(i,7)], 	# panicfrequency::Float64 -- One every 16 quarters
		durgrid[mod1(i,4)]   # panicduration::Float64 -- 8 quarters
		)
	
	basemodel=ReservesModel(basecompuparams,baseeconparams)
	# if i==1
		modelinitialize!(basemodel)
	# else
	# 	setindex!(basemodel.valuepay, tempvaluepay, :)
	# 	setindex!(basemodel.valuedefault, tempvaluedefault, :)
	# 	setindex!(basemodel.bondprice, tempbondprice, :)
	# end

	basesolverparams=SolverParams(0.25, 0, 200, 1800, 5000, false, 2e-04)
	solvereservesmodel!(basemodel, basesolverparams)	

	# basesolverparams=SolverParams(0.1, 0, 100, 200, 5000, false, 5e-05)
	# solvereservesmodel!(basemodel, basesolverparams)	

	# setindex!(tempvaluepay,basemodel.valuepay,:)
	# setindex!(tempvaluedefault,basemodel.valuedefault,:)
	# setindex!(tempbondprice,basemodel.bondprice,:)

	# Simulate model
	basesimul=ModelSimulation(100000)
	simulatemodel!(basesimul,basemodel,true)
	# 4. Obtain moments
	basemoments=ModelMoments()
	flag=getmoments!(basemoments, basesimul, basemodel.grids, 1000) # burnin 1000 observations
	println("---------------------------------------------------------------------------------------------")
	println("   index   |    debt    |  reserves  |   spravg   |   sprvar   |    defstat  |  defchoice  |")
	print("     ")
	showcompact(i)
	print("    |  ")
	showcompact(basemoments.debtmean)
	print("  | ")
	showcompact(basemoments.reservesmean)
	print("  | ")
	showcompact(basemoments.spreadmean)
	print(" | ")
	showcompact(basemoments.spreadsigma)
	print(" | ")
	showcompact(basemoments.defaultstatemean)
	print("  | ")
	showcompact(basemoments.defaultchoicemean)
	println("  |")
	println("=============================================================================================")
end	
	
end # Module end

