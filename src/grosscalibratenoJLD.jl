module grosscalibratenoJLD
#
addprocs(23)


push!(LOAD_PATH,pwd())
export basemodel

using ReservesTypes


df1=Array{Float64}(9)
df1=collect(linspace(-0.375,-0.075,9))

basecompuparams=ComputationParams(
	# Output Grid Lenght
	21,		# ynum::Int64
	# Debt grid parameters
	0.0,	# debtmin::Float64
	1.02,	# debtmax::Float64
	52,		# debtnum::Int64
	# Reserves grid parameters
	0.0,		# resmin::Float64
	1.0, 	# resmax::Float64
	51,		# resnum::Int64
	# Temporary (smoothing shock parameters)
	0.008, 	# msigma::Float64
	2.0,	# msdwidth::Float64 (For convergence)
	13,		# mnum::Int64
	-100.0,	# thrmin::Float64
	)

# Temporary arrays to store last model solution. Helps with speed
tempvaluepay=Array{Float64}(basecompuparams.debtnum, basecompuparams.resnum, basecompuparams.mnum, basecompuparams.ynum, 2)
tempvaluedefault=Array{Float64}(basecompuparams.resnum, basecompuparams.ynum, 2)
tempbondprice=Array{Float64}(basecompuparams.debtnum, basecompuparams.resnum, basecompuparams.ynum, 2)

basesolverparams=SolverParams(0.25, 0, 25, 3000, 5000, false, 1e-05)

for i=1:18

	baseeconparams=EconParams(
		# Consumer
		0.985,	# bbeta::Float64
		2,		# ggamma::Int64;  # HAS TO BE EQUAL TO 2. This cannot change. Will destroy threshold solution.
		# Risk free rate
		0.01, 	# rfree::Float64 # 4% yearly
		# Bond Maturity and coupon
		0.05, 	# llambda::Float64 # 5 year avg maturity
		0.0185, 	# coupon:: Float64 
		# Expected Output grid parameters
		0.8, 	# logOutputRho::Float64
		0.0716, 	# logOutputSigma::Float64
		# Default output cost parameters
		df1[i],# defcost1::Float64
		0.25, 	# defcost2::Float64
		0.1, # reentry::Float64
		# Sudden Stop Probability
		24.0, 	# panicfrequency::Float64 -- One every 16 quarters
		8.0   # panicduration::Float64 -- 8 quarters
		)
	
	basemodel=ReservesModel(basecompuparams,baseeconparams)
	if i==1
		modelinitialize!(basemodel)
	else
		setindex!(basemodel.valuepay, tempvaluepay, :)
		setindex!(basemodel.valuedefault, tempvaluedefault, :)
		setindex!(basemodel.bondprice, tempbondprice, :)
	end

	solvereservesmodel!(basemodel, basesolverparams)	
	
	setindex!(tempvaluepay,basemodel.valuepay,:)
	setindex!(tempvaluedefault,basemodel.valuedefault,:)
	setindex!(tempbondprice,basemodel.bondprice,:)
	# Simulate model
	basesimul=ModelSimulation(100000)
	simulatemodel!(basesimul,basemodel)
	# 4. Obtain moments
	basemoments=ModelMoments()
	flag=getmoments!(basemoments, basesimul, basemodel.grids, 1000) # burnin 1000 observations
	println(" defcost1 | debt	|  reserves  |	spravg	|	sprvar	|	defstat  |  defchoice  |")
	println("---------------------------------------------------------------------------------------------")
	showcompact(baseeconparams.defcost1)
	print("	|	")
	showcompact(basemoments.debtmean)
	print("	|	")
	showcompact(basemoments.reservesmean)
	print("	|	")
	showcompact(basemoments.spreadmean)
	print("	|	")
	showcompact(basemoments.spreadsigma)
	print("	|	")
	showcompact(basemoments.defaultstatemean)
	print("	|	")
	showcompact(basemoments.defaultchoicemean)
	println("	|	")
	println("=============================================================================================")
end	
	
end # Module end

