# This is now a script
#module serialsobol
#
Pkg.add("Sobol")
using Sobol

# addprocs(3)	# Process allocated outside
# using ClusterManagers
# ClusterManagers.addprocs_sge(59,queue="all.q",qsub_env="LD_LIBRARY_PATH")

cd("..")
push!(LOAD_PATH,pwd())
@everywhere using Base.LinAlg.BLAS
@everywhere using ReservesTypes

cd("\NewCalibration")

# 1. Define parameters
growth=0.0065		# Quarterly growth rate
pivot=0.87 		# Pivot point for default cost 
					# For y=51 sigma=0.098 Point--ergodic density--cumulative 
					# prob: 0.8467--0.0270--0.0595,  0.8706--0.0418--0.1013, 0.8950--0.0598--0.1611

# 1.1 Sobol points over calibrated parameter box 

# par1 is such that proportinoal cost at y=pivot  is (y-y_def)/y= d1+pivot*d2 = par1
# par 2 is the derivative loss at y=pivot: d y_def/ dy = 1-d1-2*pivot*d2 = 1 - par2  

# 1.1.1 One parameter vector [beta]
# calsequence=SobolSeq(1, [0.986], [0.993])

# 1.1.2 Two parameter vector [beta, W]
# calsequence=SobolSeq(2, [0.986, 1.4], [0.993, 3.0])

# 1.1.3 Three parameter vector [duration, par1, par2]
calsequence=SobolSeq(3, [2.0, 0.01, 0.3], [18.0, 0.09, 0.94])

# 1.1.4 Four parameter vector [beta, wealthmean, par1, par2]
# calsequence=SobolSeq(4, [0.986, 3.0, 0.04, 0.3 ], [0.993, 5.0, 0.08, 0.85])


# 1.2 Fixed computation parameters
basecompuparams=ComputationParams(
	# Output Grid Lenght
	51,		# ynum::Int64
	# Debt grid parameters
	0.0,	# debtmin::Float64
	1.0,	# debtmax::Float64
	41,		# debtnum::Int64
	# Reserves grid parameters
	0.0,		# resmin::Float64
	0.75, 	# resmax::Float64
	31,		# resnum::Int64
	# Temporary (smoothing shock parameters)
	0.006, 	# msigma::Float64 	
	2.5,	# msdwidth::Float64 	# This is important for curvature on mmass. 2.0 is very flat and leads to poor convergence 
	13,		# mnum::Int64
	-100.0,	# thrmin::Float64
	)

# 1.3 Solver parameters
basesolverparams=SolverParams(0.2, 0, 0, 800, 5000, false, 1e-05)

# 1.4 Itereation control and output print
iterstart=0
itermax=2048
# skip(calsequence,iterstart)


# 2. Pallalel evaluation

# 2.1 Output file initialization
outfilename="calibtwo.txt"
calout=open(outfilename,"a")
println(calout, "----------------------------------------------------------------------------------------------------------------------------------------------------")
println(calout, "                  parvec                  |   debt    |  reserves  |   spravg   |   sprvar   |    defstat  |  defchoice  | sprXgrowth |   maxgap   |")
close(calout)

# 2.2 Original parameters (pre-growth transformation)
growth=0.0065 		# Avg quarterly growth
bbeta=0.9868		# TO CALIBRATE
ggamma=2 			# HAS TO BE EQUAL TO 2. This cannot change. Will destroy threshold solution.
rfree=0.01 			# 4% yearly
# llambda=0.05 		# 20 quarter year maturity (6% avg quarterly debt service)
coupon=0.016  	# rfree + 240 b.p annualized spread

# 2.3 Parallel call ofer parameter comprehension 
pmap( momentsimulator!, Iterators.repeated(basecompuparams,itermax), 
	[ baseeconparams=EconParams(
		# Consumer
		bbeta, # parvec[1]*(1.0+growth)^(1-ggamma),			# bbeta::Float64
		ggamma,										# ggamma::Int64;  # HAS TO BE EQUAL TO 2. This cannot change. Will destroy threshold solution.
		# Risk free rate, lenders risk aversion and mean wealth
		(rfree-growth)/(1.0+growth),				# rfree::Float64 		
		ggamma,										# ggamalender::Float64 	# Same as borrower. This CANNOT be changed to values different from 2.
		4.6,									# wealthmean::Flotat64	# To calibrate with risk premium
		# Bond Maturity and coupon
		(1.0/parvec[1]+growth)/(1.0+growth), 				# llambda::Float64 		# 5 year avg maturity 
		coupon-growth*(1.0-1.0/parvec[1])/(1.0+growth), 	# coupon:: Float64 
		# Expected Output grid parameters
		0.7584, 									# logOutputRho::Float64 # gdp in USD delfacted by USPCE. HP filtered. 1994Q1-2016Q2
		0.0982,									 	# logOutputSigma::Float64
		0.12,										# govtspending::Float64
		# Default output cost parameters
		# -0.455,										# defcost1::Float64 	# Fixed instead of calibrated
		# 0.59195, 									# defcost2::Float64 	# Fixed instead of calibrated
		2*parvec[2]-parvec[3],					# defcost1::Float64
		(parvec[3]-parvec[2])/pivot,				# defcost2::Float64
		0.125,	 									# reentry::Float64
		# Sudden Stop Probability
		24.0, 										# panicfrequency::Float64 -- One every 16 quarters
		8.0   										# panicduration::Float64 -- 8 quarters
		)
	for parvec in [next(calsequence) for id=1:itermax] ],
	Iterators.repeated(basesolverparams,itermax), Iterators.repeated(outfilename,itermax) )	

calout=open(outfilename,"a")	
println(calout, "=============================================================================================")
close(calout)
#end # Module end

