module playground

export basemodel
push!(LOAD_PATH,pwd())
using ReservesTypes
using JLD

df1=Array{Float64}(9)
df1=collect(linspace(-0.3,0.0,9))

for i=6:9
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
		0.01, 	# msigma::Float64
		13,		# mnum::Int64
		-100.0,	# thrmin::Float64
		# Tolerances
		1e-5, 	# valtol::Float64 
		# Price updating step
		0.25, 	# updatespeed::Float64 
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
		df1[i],# defcost1::Float64
		0.25, 	# defcost2::Float64
		0.0385, # reentry::Float64
		# Sudden Stop Probability
		16.0, 	# panicfrequency::Float64 -- One every 16 quarters
		8.0   # panicduration::Float64 -- 8 quarters
		)
	
	basemodel=ReservesModel(basecompuparams,baseeconparams)
	
	modelinitialize!(basemodel)

	solvereservesmodel!(basemodel)	
	# Simulate model
	basesimul=ModelSimulation(100000)
	simulatemodel!(basesimul,basemodel)
	# 4. Obtain moments
	basemoments=ModelMoments()
	flag=getmoments!(basemoments, basesimul, basemodel.grids, 1000) # burnin 1000 observations
	filename=string("solved",i,".jld")
	save(filename, "basemodel", basemodel, "basesimul", basesimul, "basemoments", basemoments)
end	
	
end # Module end

