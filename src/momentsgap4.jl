function momentsgap4(params::Array{Float64,1})
	# 0. Unpack
	bbeta=0.995-0.011*exp(50.0*params[1])/(1+exp(50.0*params[1]))
	par1=0.005+0.11*exp(25.0*params[2])/(1+exp(25.0*params[2]))
	par2=0.15+0.85*exp(50.0*params[3])/(1+exp(50.0*params[3]))
	panicfrequency=36.0-24.0*exp(64.0*params[4])/(1+exp(64.0*params[4]))
	# 1. Check parameters and return if outside
	growth=0.0065 # Avg quarterly growth
	pivot=0.87
	# 2. Set parameters and solve

	compuparams=ComputationParams(
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
		0.0075, 	# msigma::Float64
		2.0,	# msdwidth::Float64 
		13,		# mnum::Int64
		-100.0,	# thrmin::Float64
		)

	econparams=EconParams(
		# Consumer
		bbeta*(1.0+growth)^(1-2),		# bbeta::Float64
		2,							# ggamma::Int64;  # HAS TO BE EQUAL TO 2. This cannot change. Will destroy threshold solution.
		# Risk free rate
		(0.01-growth)/(1.0+growth), 	# rfree::Float64 # 4% yearly
		# Bond Maturity and coupon
		(0.05+growth)/(1.0+growth), 	# llambda::Float64 # 5 year avg maturity (6% avg quarterly debt service)
		0.01575*(1+growth)-growth,		# coupon:: Float64 
		# Expected Output grid parameters
		0.7584, 	# logOutputRho::Float64
		0.0982, 	# logOutputSigma::Float64
		# Default output cost parameters
		2*par1-par2,			# defcost1::Float64
		(par2-par1)/pivot, 	# defcost2::Float64
		0.125, # reentry::Float64
		# Sudden Stop Probability
		panicfrequency, 	# panicfrequency::Float64 -- One every 24 quarters (25% of the time in panic)
		8.0   # panicduration::Float64 -- 8 quarters
		)	

	solverparams=SolverParams(
		0.25, 	# updatespeed::Float64 
		0, 		# startiternum::Int64
		0,		# interprint::Int64 
		1000,	# itermax::Int64
		4000, 	# intermediatesave::Int64
		false,	# policiesout::Bool
		1e-5 	# valtol::Float64 
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
	# 4.3 Store moments
	# if true
		outfilename="momentsgapout.txt"
		calout=open(outfilename,"a")
			showcompact(calout, econparams.bbeta)
			print(calout, " |  ")
			showcompact(calout, econparams.defcost1)
			print(calout, " |  ")
			showcompact(calout, econparams.defcost2)
			print(calout, " |  ")
			showcompact(calout, econparams.panicfrequency)
			print(calout, " |  ")
			showcompact(calout, moments.debtmean)
			print(calout, "  | ")
			showcompact(calout, moments.reservesmean)
			print(calout, "  | ")
			showcompact(calout, moments.spreadmean)
			print(calout, " | ")
			showcompact(calout, moments.spreadsigma)
			print(calout, " | ")
			showcompact(calout, moments.defaultstatemean)
			print(calout, "  | ")
			showcompact(calout, moments.defaultchoicemean)
			print(calout, "  | ")
			showcompact(calout, moments.spreadXgdp)
			print(calout, "  | ")
			showcompact(calout, moments.spreadXgrowth)
			print(calout, "  | ")
			showcompact(calout, moments.deltaspreadXgdp)
			print(calout, "  | ")
			showcompact(calout, moments.deltaspreadXgrowth)
			print(calout, "  | ")
			showcompact(calout, maximum([valuegap,pricegap,defaultgap]) )
			println(calout, "  |")
		close(calout)
	# end

	# 4.3 Get error norm
	errorvec=Array{Float64}(4)
	errorvec[1]=(moments.debtmean/0.158-1)^2
	errorvec[2]=(moments.spreadsigma/0.0105-1)^2
	errorvec[3]=(moments.defaultchoicemean/0.005-1)^2
	errorvec[4]=(moments.deltaspreadXgrowth/0.6144+1)^2
	# Weighted sum (add penalty for no convergence)
	gapnorm::Float64=10.0
	# in log first argument is the logarithmic base
	gapnorm=sum(errorvec)+max( 0.0, log(10, maximum([valuegap,pricegap,defaultgap])/solverparams.valtol/10.0) )
	return gapnorm
end

