module calibratesobolnoJLD
#
addprocs(25)


push!(LOAD_PATH,pwd())
export basemodel

using ReservesTypes
using Sobol

# 1. Define parameters

# 1.1 Sobol points over calibrated parameter box 

# Vector [beta, par1, par2, reentry, 1/freq]
calsequence=SobolSeq(4, [0.95, 0.01, 0.1, 0.025 ], [0.99, 0.09, 1.0, 0.0625])
# par1 is such that proportinoal cost at y=0.9  is (y-y_def)/y= d1+0.9d2 = par1
# par 2 is the derivative loss at y=0.9: d y_def/ dy = 1-d1-1.8d2 = 1 - par2  

# 1.2 Fixed computation parameters
basecompuparams=ComputationParams(
	# Output Grid Lenght
	25,		# ynum::Int64
	# Debt grid parameters
	0.0,	# debtmin::Float64
	1.0,	# debtmax::Float64
	41,		# debtnum::Int64
	# Reserves grid parameters
	0.0,		# resmin::Float64
	0.9, 	# resmax::Float64
	31,		# resnum::Int64
	# Temporary (smoothing shock parameters)
	0.01, 	# msigma::Float64
	2.0,	# msdwidth::Float64 (For convergence)
	13,		# mnum::Int64
	-100.0,	# thrmin::Float64
	)

# 1.3 Solver parameters
basesolverparams=SolverParams(0.25, 0, 0, 1000, 5000, false, 1e-05)

# 1.4 Itereation control and output print
iternum=0
itermax=4096
println("---------------------------------------------------------------------------------------------")
println("        parvec         |    debt    |  reserves  |   spravg   |   sprvar   |    defstat  |  defchoice  |")
	
for parvec in calsequence
	iternum+=1

	# 2. Solve model
	
	# 2.1. Set econ params
	
	baseeconparams=EconParams(
				# Consumer
		parvec[1],							# bbeta::Float64
		2,								# ggamma::Int64;  # HAS TO BE EQUAL TO 2. This cannot change. Will destroy threshold solution.
				# Risk free rate
		0.01, 							# rfree::Float64 # 4% yearly
				# Bond Maturity and coupon
		0.05, 							# llambda::Float64 # 5 year avg maturity
		0.0185, 						# coupon:: Float64 
				# Expected Output grid parameters
		0.8, 							# logOutputRho::Float64
		0.0716, 						# logOutputSigma::Float64
				# Default output cost parameters
		2*parvec[2]-parvec[3],			# defcost1::Float64
		(parvec[3]-parvec[2])/0.9,		# defcost2::Float6
		0.125, 						# reentry::Float64
				# Sudden Stop Probability
		1.0/parvec[4], 					# panicfrequency::Float64 -- One every 32 quarters
		8.0   							# panicduration::Float64 -- 8 quarters
		)
	
	# 2.2 Intialize model object
	
	basemodel=ReservesModel(basecompuparams,baseeconparams)
	modelinitialize!(basemodel)
	
	# 2.3. # Call solver routine

	solvereservesmodel!(basemodel, basesolverparams)	
	
	# 3. Simulate model with seed=true
	
	basesimul=ModelSimulation(100000)
	simulatemodel!(basesimul,basemodel,true)
	
	# 4. Obtain moments
	basemoments=ModelMoments()
	flag=getmoments!(basemoments, basesimul, basemodel.grids, 1000) # burnin 1000 observations
	
	# 5. Print relevans moments
	print("     ")
	showcompact(parvec)
	print("     |  ")
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
	println("  |")
	iternum==itermax && break
end	

println("=============================================================================================")
end # Module end

