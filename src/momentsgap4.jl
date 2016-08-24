function momentsgap4(params::Array{Float64,1})
	# 0. Unpack
	bbeta=exp(params[1])/(1+exp(params[1]))
	defcost1=-exp(params[2])/(1+exp(params[2]))
	defcost2ratio=1+exp(params[3])/(1+exp(params[3]))
	panicfrequency=1+50*exp(params[4])/(1+exp(params[4]))
	# 1. Check parameters and return if outside

	# 2. Set parameters and solve

	compuparams=ComputationParams(
		# Output Grid Lenght
		21,		# ynum::Int64
		# Debt grid parameters
		0.0,	# debtmin::Float64
		1.0,	# debtmax::Float64
		41,		# debtnum::Int64
		# Reserves grid parameters
		0.0,		# resmin::Float64
		0.90, 	# resmax::Float64
		31,		# resnum::Int64
		# Temporary (smoothing shock parameters)
		0.01, 	# msigma::Float64
		2.0,	# msdwidth::Float64 
		13,		# mnum::Int64
		-100.0,	# thrmin::Float64
		)

	econparams=EconParams(
		# Consumer
		bbeta,		# bbeta::Float64
		2,			# ggamma::Int64;  # HAS TO BE EQUAL TO 2. This cannot change. Will destroy threshold solution.
		# Risk free rate
		0.01, 		# rfree::Float64 # 4% yearly
		# Bond Maturity and coupon
		0.05, 		# llambda::Float64 # 5 year avg maturity
		0.0185, 	# coupon:: Float64 
		# Expected Output grid parameters
		0.8, 		# logOutputRho::Float64
		0.0716, 	# logOutputSigma::Float64
		# Default output cost parameters
		defcost1,	# defcost1::Float64
		-defcost1*defcost2ratio, 	# defcost2::Float64
		0.1, 		# reentry::Float64
		# Sudden Stop Probability
		panicfrequency, 		# panicfrequency::Float64 -- One every 16 quarters
		8.0   		# panicduration::Float64 -- 8 quarters
		)
	

	solverparams=SolverParams(
		0.25, 	# updatespeed::Float64 
		0, 		# startiternum::Int64
		0,		# interprint::Int64 
		1200,	# itermax::Int64
		4000, 	# intermediatesave::Int64
		false,	# policiesout::Bool
		1e-5, 	# valtol::Float64 
		)
	
	# 3. Initialize and Solve the model
	
	# 3.1 Initialize
	basemodel=ReservesModel(compuparams, econparams)
	modelinitialize!(basemodel)

	# 3.2 Solve
	valuegap::Float64=1.0
	pricegap::Float64=1.0
	defaultgap::Float64=1.0
	resiternum::Int64=0
	(resiternum,valuegap,pricegap,defaultgap)=solvereservesmodel!(basemodel, solverparams)	
	
	# 4 Results
	# 4.1. Simulate model
	simul=ModelSimulation(100000)
	simulatemodel!(simul,basemodel,true)
	# 4.2. Obtain moments
	moments=ModelMoments()
	getmoments!(moments, simul, basemodel.grids, 1000) # burnin 1000 observations
	# 4.3 Get error norm
	errorvec=Array{Float64}(4)
	errorvec[1]=(moments.debtmean/0.16-1)^2
	errorvec[2]=(moments.spreadmean/0.0345-1)^2
	errorvec[3]=(moments.defaultchoicemean/0.005-1)^2
	errorvec[4]=(moments.reservesmean/0.09-1)^2
	# Weighted sum (add penalty for no convergence)
	gapnorm::Float64=10.0
	# in log first argument is the logarithmic base
	gapnorm=2*errorvec[2]+sum(errorvec)+max( 0.0, log(10, max(valuegap,pricegap,defaultgap)/solverparams.valtol/10.0) )
	return gapnorm
end

