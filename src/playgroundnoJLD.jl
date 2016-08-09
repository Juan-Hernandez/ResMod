module playground
	addprocs(23)
	export basemodel
	push!(LOAD_PATH,pwd())
	using ReservesTypes
	# using JLD

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
		2.0,	# msdwidth::Float64 
		13,		# mnum::Int64
		-100.0,	# thrmin::Float64
		)

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
		-0.1875,# defcost1::Float64
		0.25, 	# defcost2::Float64
		0.1, # reentry::Float64
		# Sudden Stop Probability
		24.0, 	# panicfrequency::Float64 -- One every 16 quarters
		8.0   # panicduration::Float64 -- 8 quarters
		)
	
	basemodel=ReservesModel(basecompuparams,baseeconparams)
	
	modelinitialize!(basemodel)

	basesolverparams=SolverParams(
		0.25, 	# updatespeed::Float64 
		0, 		# startiternum::Int64
		25,		# interprint::Int64 
		3000,	# itermax::Int64
		4000, 	# intermediatesave::Int64
		false,	# policiesout::Bool
		1e-5, 	# valtol::Float64 
		)

	solvereservesmodel!(basemodel, basesolverparams)	
	
	# 3 Results
	# 3.1. Simulate model
	basesimul=ModelSimulation(100000)
	simulatemodel!(basesimul,basemodel)
	# 3.2. Obtain moments
	basemoments=ModelMoments()
	flag=getmoments!(basemoments, basesimul, basemodel.grids, 1000) # burnin 1000 observations
	# 3.3 Print Moments
	println("	debt	|	reserves	|	spravg		|	sprvar		|	defstat		|	defchoice	|")
	println("---------------------------------------------------------------------------------------------")
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
	# 3.4. Save solved
	#@save "firstsolved.jld" basemodel basesimul basemoments
	
	
end

