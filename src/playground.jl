module playground
	export basemodel
	push!(LOAD_PATH,pwd())
	using ReservesTypes
	using JLD

	basecompuparams=ComputationParams(
		# Output Grid Lenght
		21,		# ynum::Int64
		# Debt grid parameters
		0.0,	# debtmin::Float64
		1.02,	# debtmax::Float64
		35,		# debtnum::Int64
		# Reserves grid parameters
		0.0,		# resmin::Float64
		1.0, 	# resmax::Float64
		26,		# resnum::Int64
		# Temporary (smoothing shock parameters)
		0.01, 	# msigma::Float64
		1.6,	# msdwidth::Float64
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
		0.03, 	# coupon:: Float64 
		# Expected Output grid parameters
		0.93, 	# logOutputRho::Float64
		0.076, 	# logOutputSigma::Float64
		# Default output cost parameters
		-0.1875,# defcost1::Float64
		0.25, 	# defcost2::Float64
		0.0385, # reentry::Float64
		# Sudden Stop Probability
		16.0, 	# panicfrequency::Float64 -- One every 16 quarters
		8.0   # panicduration::Float64 -- 8 quarters
		)
	
	basemodel=ReservesModel(basecompuparams,baseeconparams)
	
	modelinitialize!(basemodel)

	basesolverparam=SolverParams(
		0.25, 	# updatespeed::Float64 
		0, 		# startiternum::Int64
		25,		# interprint::Int64 
		6001,	# itermax::Int64
		1000, 	# intermediatesave::Int64
		false,	# policiesout::Bool
		1e-5, 	# valtol::Float64 
		)

	solvereservesmodel!(basemodel)
		
	# Simulate model
	basesimul=ModelSimulation(100000)
	simulatemodel!(basesimul,basemodel)
	# 4. Obtain moments
	basemoments=ModelMoments()
	flag=getmoments!(basemoments, basesimul, basemodel.grids, 1000) # burnin 1000 observations

	@save "firstsolved.jld" basemodel basesimul basemoments
	
	
end

