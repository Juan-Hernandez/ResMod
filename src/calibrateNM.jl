module calibrateNM
	addprocs(25)
	# export basemodel
	push!(LOAD_PATH,pwd())
	using ReservesTypes
	using Optim
	
	include("momentsgap3.jl")
	# Create simplex
	immutable StepSimplexer <: Optim.Simplexer
		a::Array{Float64,1}
	end
	
	StepSimplexer(;a1 = 0.01, a2 = -0.5125, a3=-0.2) = StepSimplexer([a1, a2, a3])

	function Optim.simplexer{T, N}(simplexstep::StepSimplexer, initial_x::Array{T, N})
		n = length(initial_x)
		initial_simplex = Array{T, N}[initial_x for i = 1:n+1]
		for j = 1:n
			initial_simplex[j+1][j] += simplexstep.a[j]
		end
		initial_simplex
	end
	#mystep=StepSimplexer()
	#showcompact(mystep)
	# mysimplex=Optim.simplexer(mystep, [0.965, -0.1, 1.4])
	
	calibration=optimize(momentsgap3, [log(0.97/0.02), log(0.0875/0.8625), log(0.32/0.18)], NelderMead(), OptimizationOptions(show_trace=true) )
	
	showcompact(Optim.minimizer(calibration))
	bbetastar=Optim.minimizer(calibration)[1]
	bbetastar=0.99*exp(bbetastar)/(1+exp(bbetastar))
	defcost1star=Optim.minimizer(calibration)[2]
	defcost1star=-0.1-0.95*exp(defcost1star)/(1+exp(defcost1star))
	df2ratiostar=Optim.minimizer(calibration)[3]
	df2ratiostar=1+0.5*exp(df2ratiostar)/(1+exp(df2ratiostar))
	println(" \n Star parameters:")
	showcompact([bbetastar, defcost1star, df2ratiostar])
	
	basecompuparams=ComputationParams(
		# Output Grid Lenght
		25,		# ynum::Int64
		# Debt grid parameters
		0.0,	# debtmin::Float64
		1.0,	# debtmax::Float64
		41,		# debtnum::Int64
		# Reserves grid parameters
		0.0,	# resmin::Float64
		0.90, 	# resmax::Float64
		31,		# resnum::Int64
		# Temporary (smoothing shock parameters)
		0.01, 	# msigma::Float64
		2.0,	# msdwidth::Float64 
		13,		# mnum::Int64
		-100.0,	# thrmin::Float64
		)

	stareconparams=EconParams(
		# Consumer
		bbetastar,	# bbeta::Float64
		2,		# ggamma::Int64;  # HAS TO BE EQUAL TO 2. This cannot change. Will destroy threshold solution.
		# Risk free rate
		0.01, 	# rfree::Float64 # 4% yearly
		# Bond Maturity and coupon
		0.05, 	# llambda::Float64 # 5 year avg maturity
		0.0185, 	# coupon:: Float64 
		# Expected Output grid parameters
		0.8, 	# logOutputRho::Float64
		0.0716, 	# logOutputSigma::Float64
		# Default output cost parameters
		defcost1star,# defcost1::Float64
		-defcost1star*df2ratiostar, 	# defcost2::Float64
		0.1, # reentry::Float64
		# Sudden Stop Probability
		32.0, 	# panicfrequency::Float64 -- One every 16 quarters
		6.0   # panicduration::Float64 -- 8 quarters
		)
	
	starmodel=ReservesModel(basecompuparams,stareconparams);
	
	modelinitialize!(starmodel)

	basesolverparams=SolverParams(
		0.25, 	# updatespeed::Float64 
		0, 		# startiternum::Int64
		200,	# interprint::Int64 
		1200,	# itermax::Int64
		4000, 	# intermediatesave::Int64
		false,	# policiesout::Bool
		1e-5, 	# valtol::Float64 
		)

	solvereservesmodel!(starmodel, basesolverparams)	
	
	# 3 Results
	# 3.1. Simulate model
	starsimul=ModelSimulation(100000)
	simulatemodel!(starsimul,starmodel,true)
	# 3.2. Obtain moments
	starmoments=ModelMoments()
	getmoments!(starmoments, starsimul, starmodel.grids, 1000) # burnin 1000 observations
	# 3.3 Print Moments
	println("\n    debt	|	reserves	|	spravg		|	sprvar		|	defstat		|	defchoice	|")
	println("---------------------------------------------------------------------------------------------")
	showcompact(starmoments.debtmean)
	print("	|	")
	showcompact(starmoments.reservesmean)
	print("	|	")
	showcompact(starmoments.spreadmean)
	print("	|	")
	showcompact(starmoments.spreadsigma)
	print("	|	")
	showcompact(starmoments.defaultstatemean)
	print("	|	")
	showcompact(starmoments.defaultchoicemean)
	println("	|	")
	println("=============================================================================================")
	# 3.4. Save solved
	#@save "firstsolved.jld" basemodel basesimul basemoments
	
	
end

