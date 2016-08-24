module startmodel

	export basemodel, bsecompuparams, baseeconparams
	push!(LOAD_PATH,pwd())
	using ReservesTypes

	basecompuparams=ComputationParams(
		# Output Grid Lenght
		25,		# ynum::Int64
		# Debt grid parameters
		0.0,	# debtmin::Float64
		1.0,	# debtmax::Float64
		41,		# debtnum::Int64
		# Reserves grid parameters
		0.0,	# resmin::Float64
		0.90, 	# resmax::Float64
		31,		# resnum::Int64
		# Temporary (smoothing shock parameters)
		0.01, 	# msigma::Float64
		2.0,	# msdwidth::Float64 
		13,		# mnum::Int64
		-100.0,	# thrmin::Float64
		)

	baseeconparams=EconParams(
		# Consumer
		0.96747,	# bbeta::Float64
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
		-0.3042894,# defcost1::Float64
		0.3042894*1.3251468, 	# defcost2::Float64
		0.1, # reentry::Float64
		# Sudden Stop Probability
		32.0, 	# panicfrequency::Float64 -- One every 16 quarters
		6.0   # panicduration::Float64 -- 8 quarters
		)
	
	basemodel=ReservesModel(basecompuparams,baseeconparams);
	
	modelinitialize!(basemodel);
	
	end