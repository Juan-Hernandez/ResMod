module d1d2calibratenoJLD
#
addprocs(25)


push!(LOAD_PATH,pwd())
export basemodel

using ReservesTypes

# par1 is such that proportinoal cost at y=0.9  is (y-y_def)/y= d1+0.9d2 = par1
# par 2 is the derivative loss at y=0.9: d y_def/ dy = 1-d1-1.8d2 = 1 - par2  

par1=collect(linspace(0.01, 0.1, 10))
par2=collect(linspace(0.05, 1.0, 20))'


def1grid=Array{Float64}(10,20)
def2grid=Array{Float64}(10,20)
sizegrid=length(def1grid)

# def1 is just 2*par1-par2
broadcast!(+, def1grid, 2*par1, -par2 )
# def2=(par2-par1)/0.9
broadcast!(+, def2grid, -par1./0.9, par2./0.9 ) 

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

for i=1:sizegrid
	
	# 2. Solve model
	
	# 2.1. Set econ params
	
	baseeconparams=EconParams(
				# Consumer
		0.985,							# bbeta::Float64
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
		def1grid[i],					# defcost1::Float64
		def2grid[i], 					# defcost2::Float6
		0.1, 							# reentry::Float64
				# Sudden Stop Probability
		32.0, 							# panicfrequency::Float64 -- One every 32 quarters
		8.0   							# panicduration::Float64 -- 8 quarters
		)
	
	# 2.2 Intialize model object
	
	basemodel=ReservesModel(basecompuparams,baseeconparams)
	modelinitialize!(basemodel)
	
	# 2.3. Set solver parameters and call solver routine
	basesolverparams=SolverParams(0.25, 0, 200, 1000, 5000, false, 1e-05)

	solvereservesmodel!(basemodel, basesolverparams)	
	
	# 3. Simulate model
	
	basesimul=ModelSimulation(100000)
	simulatemodel!(basesimul,basemodel,true)
	
	# 4. Obtain moments
	basemoments=ModelMoments()
	flag=getmoments!(basemoments, basesimul, basemodel.grids, 1000) # burnin 1000 observations
	
	# 5. Print relevans moments
	println("---------------------------------------------------------------------------------------------")
	println("  index  |    debt    |  reserves  |   spravg   |   sprvar   |    defstat  |  defchoice  |")
	print("     ")
	showcompact(i)
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
	println("=============================================================================================")
end	
	
end # Module end

