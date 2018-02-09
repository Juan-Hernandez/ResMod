module playground
	export basemodel
	push!(LOAD_PATH,pwd())
	using ReservesTypes
	using JLD
	
	growth=0.0065 # Avg quarterly growth
	basecompuparams=ComputationParams(
		# Output Grid Lenght
		50,		# ynum::Int64
		# Debt grid parameters
		0.0,	# debtmin::Float64
		1.0,	# debtmax::Float64
		41,		# debtnum::Int64
		# Reserves grid parameters
		0.0,		# resmin::Float64
		1.0, 	# resmax::Float64
		41,		# resnum::Int64
		# Temporary (smoothing shock parameters)
		0.008, 	# msigma::Float64
		2.0,	# msdwidth::Float64u
		13,		# mnum::Int64
		-100.0,	# thrmin::Float64
		)

		baseeconparams=EconParams(
			# Consumer
			0.9895*(1.0+growth)^(1-2),		# bbeta::Float64
			2,								# ggamma::Int64;  # HAS TO BE EQUAL TO 2. This cannot change. Will destroy threshold solution.
			# Risk free rate
			(0.01-growth)/(1.0+growth), 	# rfree::Float64 # 4% yearly
			# Bond Maturity and coupon
			(0.05+growth)/(1.0+growth), 	# llambda::Float64 # 5 year avg maturity (6% avg quarterly debt service)
			0.01575*(1+growth)-growth,		# coupon:: Float64 
			# Expected Output grid parameters
			0.7584, 	# logOutputRho::Float64
			0.0982, 	# logOutputSigma::Float64
			# Default output cost parameters
			-0.458,		# defcost1::Float64
			0.59655, 	# defcost2::Float64
			0.125, 		# reentry::Float64
			# Sudden Stop Probability
			24.0, #19.5122, 	# panicfrequency::Float64 -- One every 24 quarters (25% of the time in panic)
			8.0   # panicduration::Float64 -- 8 quarters
			)
	
	basemodel=ReservesModel(basecompuparams,baseeconparams)
	
	modelinitialize!(basemodel)

	basesolverparams=SolverParams(
		0.25, 	# updatespeed::Float64 
		0, 		# startiternum::Int64
		20,		# interprint::Int64 
		1000,	# itermax::Int64
		2000, 	# intermediatesave::Int64
		false,	# policiesout::Bool
		1e-5, 	# valtol::Float64 
		)

	(resiternum,valuegap,pricegap,defaultgap)=solvereservesmodel!(basemodel, basesolverparams)	
	
	# 3 Results
	# 3.1. Simulate model
	basesimul=ModelSimulation(100000)
	simulatemodel!(basesimul,basemodel,true)
	# Save simulated
	@save "firstsim.jld" basemodel basesimul

	# 3.2. Obtain moments
	basemoments=ModelMoments()
	flag=getmoments!(basemoments, basesimul, basemodel.grids, 1000) # burnin 1000 observations
	# 3.3 Print Moments
	println("	debt	|	reserves	|	spravg		|	sprvar		|	defstat		|	defchoice	| sprXgrowth |   maxgap   |")
	println("----------------------------------------------------------------------------------------------------------------------")
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
		print("  | ")
		showcompact(basemoments.spreadXgdp)
		print("  | ")
		showcompact(basemoments.spreadXgrowth)
		print("  | ")
		showcompact(basemoments.deltaspreadXgdp)
		print("  | ")
		showcompact(basemoments.deltaspreadXgrowth)
		print("  | ")
		showcompact(maximum([valuegap,pricegap,defaultgap]) )
		println("  |")
	println("======================================================================================================================")
	# 3.4. Save solved
	@save "firstsolved.jld" basemodel basesimul basemoments
		
end

