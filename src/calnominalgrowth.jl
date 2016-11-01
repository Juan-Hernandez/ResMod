# This is now a script
#module serialsobol
#
# addprocs(25)


push!(LOAD_PATH,pwd())
# export basemodel

@everywhere using Base.LinAlg.BLAS
@everywhere using ReservesTypes
using Sobol


@everywhere include("serialmomentsimulator!.jl")
@everywhere include("solvereservesmodelserial!.jl")
# 1. Define parameters
growth=0.0114		# Quarterly growth rate
pivot=0.86 	# Pivot point for default cost
# 1.1 Sobol points over calibrated parameter box 

# Vector [beta, par1, par2, 1/freq]
calsequence=SobolSeq(4, [0.984, 0.0, 0.1, 0.025 ], [1.0, 0.08, 0.9, 0.065])
# par1 is such that proportinoal cost at y=pivot  is (y-y_def)/y= d1+pivot*d2 = par1
# par 2 is the derivative loss at y=pivot: d y_def/ dy = 1-d1-2*pivot*d2 = 1 - par2  

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
	0.75, 	# resmax::Float64
	31,		# resnum::Int64
	# Temporary (smoothing shock parameters)
	0.00625, 	# msigma::Float64
	2.0,	# msdwidth::Float64 (For convergence)
	13,		# mnum::Int64
	-100.0,	# thrmin::Float64
	)

# 1.3 Solver parameters
basesolverparams=SolverParams(0.25, 0, 0, 1000, 5000, false, 1e-05)

# 1.4 Itereation control and output print
iternum=0
itermax=256
outfilename="rho84sigma022new.txt"
calout=open(outfilename,"a")
println(calout, "----------------------------------------------------------------------------------------------------------------------------------------------------")
println(calout, "                  parvec                  |   debt    |  reserves  |   spravg   |   sprvar   |    defstat  |  defchoice  | sprXgrowth |   maxgap   |")
close(calout)

pmap( serialmomentsimulator!, repeated(basecompuparams,itermax), 
	[ baseeconparams=EconParams(
		parvec[1]*(1+growth)^(1-2),		# bbeta::Float64
		2,								# ggamma::Int64;  # HAS TO BE EQUAL TO 2. This cannot change. Will destroy threshold solution.
		(0.01-growth)/(1+growth), 		# rfree::Float64 # 4% nominal yearly
		(0.05+growth)/(1+growth), 		# llambda::Float64 # 5 year avg maturity
		0.0085+(0.01-growth)/(1+growth),# coupon:: Float64 # 340bp divided by 4 plus r
		0.7909,							# logOutputRho::Float64
		0.0747, 						# logOutputSigma::Float64
		2*parvec[2]-parvec[3],			# defcost1::Float64
		(parvec[3]-parvec[2])/pivot,	# defcost2::Float6
		0.125,	 						# reentry::Float64
		1.0/parvec[4], 					# panicfrequency::Float64 -- One every 32 quarters
		8.0   							# panicduration::Float64 -- 8 quarters
		)
	for parvec in [next(calsequence) for id=1:itermax] ],
	repeated(basesolverparams,itermax), repeated(outfilename,itermax) )	

calout=open(outfilename,"a")	
println(calout, "=============================================================================================")
close(calout)
#end # Module end

