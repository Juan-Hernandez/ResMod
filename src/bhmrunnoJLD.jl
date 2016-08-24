module bhmrunnoJLD
	addprocs(23)
	export basemodel
	push!(LOAD_PATH,pwd())
	using ReservesTypes
	# using JLD

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
		2.0,	# msdwidth::Float64 
		13,		# mnum::Int64
		-100.0,	# thrmin::Float64
		)

	baseeconparams=EconParams(
		# Consumer
		0.98,	# bbeta::Float64
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
		-1.017,# defcost1::Float64
		1.19, 	# defcost2::Float64
		0.083, # reentry::Float64
		# Sudden Stop Probability
		24.0, 	# panicfrequency::Float64 -- One every 16 quarters
		8.0   # panicduration::Float64 -- 8 quarters
		)
	
	basemodel=ReservesModel(basecompuparams,baseeconparams)
	
	modelinitialize!(basemodel)

	basesolverparams=SolverParams(
		0.25, 	# updatespeed::Float64 
		0, 		# startiternum::Int64
		100,	# interprint::Int64 
		1800,	# itermax::Int64
		4000, 	# intermediatesave::Int64
		false,	# policiesout::Bool
		1e-4, 	# valtol::Float64 
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

