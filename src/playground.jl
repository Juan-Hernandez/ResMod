# module playground #script
#	export basemodel
	#currdir=pwd()
	#cd("..")
	#push!(LOAD_PATH,pwd())
	#cd(currdir)
	using ReservesTypes
	using JLD
	filename="playground.jld"	
	basecompuparams=ComputationParams(
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
		0.006, 	# msigma::Float64
		2.5,	# msdwidth::Float64u
		13,		# mnum::Int64
		-100.0,	# thrmin::Float64
		)
	# Original params
	bbeta=0.9895
	ggamma=2 			# HAS TO BE EQUAL TO 2. This cannot change. Will destroy threshold solution.
	growth=0.0065 		# Avg quarterly growth
	rfree=0.01 			# 4% yearly
	llambda=0.05 		# 20 quarter year maturity
	coupon=0.01575  	# rfree + 230 b.p annualized spread 
		baseeconparams=EconParams(
			# Consumer
			bbeta*(1.0+growth)^(1-ggamma),			# bbeta::Float64
			ggamma,									# ggamma::Int64;  # HAS TO BE EQUAL TO 2. This cannot change. Will destroy threshold solution.
			# Risk free rate
			(rfree-growth)/(1.0+growth),		 	# rfree::Float64 # 4% yearly
			# Bond Maturity and coupon
			(llambda+growth)/(1.0+growth), 		# llambda::Float64 # 5 year avg maturity (6% avg quarterly debt service)
			coupon-growth*(1.0-llambda)/(1.0+growth), # coupon:: Float64 
			#  ll + z = ll'+z' = (ll+g)/(1+g)+ z'
 			# z' = z + ll - ll' = z + (ll + ll*g)/(1+g) -(ll + g)/(1+g) = z -g(1-ll)/(1+g) = z - g*(1-ll')
			# OLD z' :  0.01575*(1+growth)-growth,
			#  ll + (1-ll)z = ll'+(1-ll')z' = (ll+g)/(1+g)+ (1-ll)z'/(1+g)
 			# z' = z(1+g) + ll(1+g)/(1-ll) - ll'(1+g)/(1-ll) = z(1+g) + (ll + ll*g)/(1-ll) -(ll + g)/(1-ll) = z(1+g) -g(1-ll)/(1-ll) = z(1+g) - g
 			# TEST: since I changed to coupon in last period, setting the new z as (old z)*(1-llambda) should leave everithing equal
 			# 0.01575*(1-0.05)-growth*(1-0.05)/(1+growth),
			# Expected Output grid parameters
			0.7584, 	# logOutputRho::Float64
			0.0982, 	# logOutputSigma::Float64
			0.0,		# govtspending::Float64
			# Default output cost parameters
			-0.455,		# defcost1::Float64
			0.59195, 	# defcost2::Float64
			0.125, 		# reentry::Float64
			# Sudden Stop Probability
			24.0, #19.5122, 	# panicfrequency::Float64 -- One every 24 quarters (25% of the time in panic)
			8.0   # panicduration::Float64 -- 8 quarters
			)
	
	basemodel=ReservesModel(basecompuparams,baseeconparams)
	
	modelinitialize!(basemodel)

	basesolverparams=SolverParams(
		0.05, 	# updatespeed::Float64 		# 0.25 generates some convergence problems.
		0, 		# startiternum::Int64
		20,		# interprint::Int64 
		801,	# itermax::Int64
		800, 	# intermediatesave::Int64
		false,	# policiesout::Bool
		1e-5, 	# valtol::Float64 
		)

	(resiternum,valuegap,pricegap,defaultgap)=solvereservesmodel!(basemodel, basesolverparams)	
	jldopen(filename,"w") do file
		write(file, "basemodel", basemodel)
	end	
	# 3 Results
	# 3.1. Simulate model
	basesimul=ModelSimulation(100000)
	simulatemodel!(basesimul,basemodel,true)
	# Save simulated
	jldopen(filename,"r+") do file
		write(file, "basesimul", basesimul)
	end
	# 3.2. Obtain moments
	basemoments=ModelMoments()
	flag=getmoments!(basemoments, basesimul, basemodel.grids, 1000) # burnin 1000 observations
	# 3.3 Print Moments
	println("	debt	|	reserves	|	spravg		|	sprvar		|	defstat		|	defchoice	| sprXgrowth |   maxgap   |")
	println("----------------------------------------------------------------------------------------------------------------------")
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
	print("  | ")
	showcompact(basemoments.spreadXgdp)
	print("  | ")
	showcompact(basemoments.spreadXgrowth)
	print("  | ")
	showcompact(basemoments.deltaspreadXgdp)
	print("  | ")
	showcompact(basemoments.deltaspreadXgrowth)
	print("  | ")
	showcompact(maximum([valuegap,pricegap,defaultgap]) )
	println("  |")
	println("======================================================================================================================")
	# 3.4. Save solved
	jldopen(filename,"r+") do file
		write(file, "basemoments", basemoments)
	end	
	flag=errr	
# end

