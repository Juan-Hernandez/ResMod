	module bhmrun
	export basemodel
	push!(LOAD_PATH,pwd())
	using JLD
	using ReservesTypes
	basecompuparams=ComputationParams(
		# Output Grid Lenght
		21,		# ynum::Int64
		# Debt grid parameters
		0.0,	# debtmin::Float64
		2.4,	# debtmax::Float64
		61,		# debtnum::Int64
		# Reserves grid parameters
		0.0,		# resmin::Float64
		1.4, 	# resmax::Float64
		36,		# resnum::Int64
		# Temporary (smoothing shock parameters)
		0.00625, 	# msigma::Float64 # Now the extremes are at +/- 1.6 sd
		13,		# mnum::Int64
		-100.0,	# thrmin::Float64
		# Tolerances
		5.0*1e-5, 	# valtol::Float64 
		# Price updating step
		1.0/3.0, 	# updatespeed::Float64 
		)
	baseeconparams=EconParams(
		# Consumer
		0.985,	# bbeta::Float64
		2,		# ggamma::Int64;  # HAS TO BE EQUAL TO 2. This cannot change. Will destroy threshold solution.
		# Risk free rate
		0.01, 	# rfree::Float64 # 4% yearly
		# Bond Maturity and coupon
		0.05, 	# llambda::Float64 # 5 year avg maturity
		0.02, 	# coupon:: Float64 
		# Expected Output grid parameters
		0.93, 	# logOutputRho::Float64
		0.05, 	# logOutputSigma::Float64
		# Default output cost parameters
		-1.017,# defcost1::Float64
		1.18, 	# defcost2::Float64
		1.0/16.0, # reentry::Float64 16 quarters in default
		# Sudden Stop Probability
		0.0625, 	# norm2ss::Float64 This is frequency: Panics once every 4 years (16 quarters)
		1.0-1.0/8.0   # ss2ss::Float 64 Average duration of panic 8 quarters (2 years)
		)
	#1. Create and initialize model object
	basemodel=ReservesModel(basecompuparams,baseeconparams)
	modelinitialize!(basemodel)
	# 2. Solve model
	flag=solvereservesmodel!(basemodel)	
	jldopen("solvedmodel.jld", "w") do file
				write(file, "basemodel", basemodel)
				write(file, "flag", flag)
	end
	
	# 3. Simulate model
	basesimul=ModelSimulation(100000)
	simulatemodel!(basesimul,basemodel)
	# 4. Obtain moments
	basemoments=ModelMoments()
	flag=getmoments!(basemoments, basesimul, basemodel.grids, 1000) # burnin 1000 observations
	jldopen("simulatedmodel.jld", "w") do file
		write(file, "basemodel", basemodel)
		write(file, "basesimul", basesimul)
		write(file, "basemoments", basemoments)
	end
		
end

