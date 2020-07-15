# 0. Output file name
outfilename="cali3reg.txt"
# 1. Fixed solver parameters
basesolverparams=SolverParams(0.2, 0, 0, 800, 5000, false, 1e-05, false)
# SolverParams(updatespeed, startiternum, iterprint, itermax, intermediatesave, policiesout, valtol, debugbool)

# 2. Fixed computation parameters
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

# 3. Sobol parameters 

# par1 is such that proportinoal cost at y=pivot  is (y-y_def)/y= d1+pivot*d2 = par1
# par 2 is the derivative loss at y=pivot: d y_def/ dy = 1-d1-2*pivot*d2 = 1 - par2  

# 3.1 One parameter vector [beta]
# calsequence=SobolSeq(1, [0.986], [0.993])

# 3.2 Two parameter vector [beta, W]
# calsequence=SobolSeq(2, [0.986, 1.4], [0.993, 3.0])

# 3.3 Three parameter vector [beta, par1, par2]
# calsequence=SobolSeq(3, [0.9864, 0.0, 0.3], [0.992, 0.08, 0.94])

# 3.4 Four parameter vector [beta, riskdur, par1, par2]
calsequence=SobolSeq(4, [0.986, 2.0, 0.04, 0.3 ], [0.993, 10.0, 0.08, 0.85])

# 3.5 Itereation control and output print
iterstart=0
itermax=4

# 4. Economic parameters
# 4.1 Metaparameters
	growth=0.0065		# Quarterly growth rate
	pivot=0.87 			# Pivot point for default cost 
						# For y=51 sigma=0.098 Point--ergodic density--cumulative 
						# prob: 0.8467--0.0270--0.0595,  0.8706--0.0418--0.1013, 0.8950--0.0598--0.1611
	growth=0.0065 		# Avg quarterly growth


# 4.2 Non sobol Economic parameters (pre-growth transformation)
	# Consumer
		# bbeta=parvec[1]*(1.0+growth)^(1-ggamma)	# bbeta::Float64
		ggamma=2		# ggamma::Int64;  # This cannot change. Will destroy threshold solution. 
	# Risk free rate, lenders risk aversion and mean wealth
		rfree=0.01		# rfree::Float64
		gammalender=0.0	# gammalender::Float64  
		wealthmean=1.0	# wealthmean::Flotat64   
	# Bond Maturity and coupon
		llambda=0.05 	# llambda::Float64	# 20 quarter year maturity (6% avg quarterly debt service)
		coupon=0.01575  # coupon:: Float64	# rfree + 230 b.p annualized spread
	# Expected Output grid parameters
		logyrho=0.7584	# logOutputRho::Float64 # gdp in USD delfact$
		logysig=0.0982	# logOutputSigma::Float64
		govspend=0.12	# govtspending::Float64
	# Default output cost parameters
		# defcost1=2*parvec[3]-parvec[4]			# defcost1::Float64     
		# defcost2=(parvec[4]-parvec[3])/pivot		# defcost2::Float64     
		reentry=0.125	# reentry::Float64
	# Sudden Stop Probability
		safedur=24.0	# safeduration::Float64		# Expected duration of 6 years.
		# riskdur=parvec[2]	# riskduration::Float64		# Expected time in risk regime one and a half years
		panicdur=4.0	# panicduration::Float64 -- 4 quarter
		panicfreq=16.0	# panicfrequency::Float64 -- One every 8 years 


# # 5. Test: one call
# parvec=next!(calsequence)
# testeconparams=EconParams(
# 		# Consumer
# 		parvec[1]*(1.0+growth)^(1-ggamma),			# bbeta::Float64
# 		ggamma,										# ggamma::Int64;  # HAS TO BE EQUAL TO 2. This cannot change. Will destroy threshold solution.
# 		# Risk free rate, lenders risk aversion and mean wealth
# 		(rfree-growth)/(1.0+growth),				# rfree::Float64 		
# 		parvec[2],									# ggamalender::Float64 	
# 		wealthmean,									# wealthmean::Flotat64	# To calibrate with risk premium
# 		# Bond Maturity and coupon
# 		(llambda+growth)/(1.0+growth), 				# llambda::Float64 		# 5 year avg maturity 
# 		coupon-growth*(1.0-llambda)/(1.0+growth), 	# coupon:: Float64 
# 		# Expected Output grid parameters
# 		logyrho, 									# logOutputRho::Float64 # gdp in USD delfacted by USPCE. HP filtered. 1994Q1-2016Q2
# 		logysig,								 	# logOutputSigma::Float64
# 		govspend,									# govtspending::Float64
# 		# Default output cost parameters
# 		2*parvec[3]-parvec[4],						# defcost1::Float64
# 		(parvec[4]-parvec[3])/pivot,				# defcost2::Float64
# 		reentry,									# reentry::Float64
# 		# Sudden Stop Probability
# 		ssfreq,										# panicfrequency::Float64 -- One every 16 quarters
# 		ssdur  										# panicduration::Float64 -- 8 quarters
# 		)
# momentsimulator!(basecompuparams, testeconparams, basesolverparams, outfilename)
# println("test passed")

