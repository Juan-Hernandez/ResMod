# This is now a script
#module serialsobol

using Sobol
using Distributed
# addprocs(3)	# Process allocated outside
# using ClusterManagers
# ClusterManagers.addprocs_sge(59,queue="all.q",qsub_env="LD_LIBRARY_PATH")

@everywhere cd("..")
@everywhere push!(LOAD_PATH,pwd())
@everywhere using LinearAlgebra
@everywhere using ReservesTypes
@everywhere cd(".\\NewCalibration")

# 1. Insrt parameters 

include("calitestparameters.jl")

# 2. Pallalel evaluation

# 2.1 Output file initialization
outfilename="calibzero2.txt"
calout=open(outfilename,"a")
println(calout, "----------------------------------------------------------------------------------------------------------------------------------------------------")
println(calout, "                  parvec                  |   debt    |  reserves  |   spravg   |   sprvar   |    defstat  |  defchoice  | sprXgrowth |   maxgap   |")
close(calout)

# 2.3 Parallel call ofer parameter comprehension 
pmap( momentsimulator!, Iterators.repeated(basecompuparams,itermax), 
	[ baseeconparams=EconParams(
		# Consumer
		parvec[1]*(1.0+growth)^(1-ggamma),			# bbeta::Float64
		ggamma,										# ggamma::Int64;  # HAS TO BE EQUAL TO 2. This cannot change. Will destroy threshold solution.
		# Risk free rate, lenders risk aversion and mean wealth
		(rfree-growth)/(1.0+growth),				# rfree::Float64 		
		parvec[2],									# ggamalender::Float64 	
		wealthmean,									# wealthmean::Flotat64	# To calibrate with risk premium
		# Bond Maturity and coupon
		(llambda+growth)/(1.0+growth), 				# llambda::Float64 		# 5 year avg maturity 
		coupon-growth*(1.0-llambda)/(1.0+growth), 	# coupon:: Float64 
		# Expected Output grid parameters
		logyrho, 									# logOutputRho::Float64 # gdp in USD delfacted by USPCE. HP filtered. 1994Q1-2016Q2
		logysig,								 	# logOutputSigma::Float64
		govspend,									# govtspending::Float64
		# Default output cost parameters
		2*parvec[3]-parvec[4],						# defcost1::Float64
		(parvec[4]-parvec[3])/pivot,				# defcost2::Float64
		reentry,									# reentry::Float64
		# Sudden Stop Probability
		ssfreq,										# panicfrequency::Float64 -- One every 16 quarters
		ssdur  										# panicduration::Float64 -- 8 quarters
		)
	for parvec in [next!(calsequence) for id=1:itermax] ],
	Iterators.repeated(basesolverparams,itermax), Iterators.repeated(outfilename,itermax) )	

calout=open(outfilename,"a")	
println(calout, "=============================================================================================")
close(calout)
#end # Module end

