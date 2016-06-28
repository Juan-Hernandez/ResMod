module ReservesTypes

using Base.LinAlg.BLAS
using JLD

export ComputationParams, EconParams, ModelGrids, ReservesModel, ModelSimulation, ModelMoments
export modelinitialize!, valuepayexpectation!, solvedefaultvalue!, solvereservesmodel!
export updatevaluepaydrm!, mthresholds!, mergethresholds!, defaultthresholds!
export suddenstopthresholds!, integratethresholds!, priceexpectation!, getpolicies!
export simulatemodel!, moutputregimeexpectation!, getmoments!, reservesfigures


immutable ComputationParams
	# Output grid lenght
	ynum::Int64
	# Debt grid parameters
	debtmin::Float64
	debtmax::Float64
	debtnum::Int64
	# Reserves grid parameters
	resmin::Float64
	resmax::Float64
	resnum::Int64
	# Temporary (smoothing shock parameters)
	msigma::Float64
	mnum::Int64
	thrmin::Float64
	# Tolerances
	valtol::Float64 
	# Price updating step
	updatespeed::Float64
end

immutable EconParams
	# Consumer
	bbeta::Float64
	ggamma::Int64  # HAS TO BE EQUAL TO 2. This cannot change. Will destroy threshold solution.
	# Risk free rate
	rfree::Float64 
	# Bond Maturity and coupon
	llambda::Float64 # inverse of quarterly avg maturity
	coupon:: Float64 
	# Total Output correlation and variance
	logoutputrho::Float64
	logoutputsigma::Float64
	# Default output cost parameters
	defcost1::Float64
	defcost2::Float64
	reentry::Float64
	# Sudden Stop transition Probability
	norm2ss::Float64
	ss2ss::Float64
end

include("rouwenhorst.jl")
include("makegrids.jl")

immutable ModelGrids
	# Output
	y::Array{Float64,1}
	ytrans::Array{Float64,2}
	yergodic::Array{Float64,1}
	ydefault::Array{Float64,1}
	# m-shock
	mmidpoints::Array{Float64,1}
	mextremes::Array{Float64,1}
	mmass::Array{Float64,1}
	# debt
	debt::Array{Float64,1}
	debtmaxss::Array{Int64,1}
	# reserves
	reserves::Array{Float64,1}	
end

immutable ModelPolicies
	reservesindefault::Array{Int64,3}
	debt::Array{Int64,5}
	reserves::Array{Int64,5}
	default::BitArray{5}
	function ModelPolicies(compuparams::ComputationParams)
		reservesindefault=Array{Int64}(compuparams.resnum, compuparams.ynum, 2)
		debt=Array{Int64}(compuparams.debtnum, compuparams.resnum, compuparams.mnum, compuparams.ynum, 2)
		reserves=Array{Int64}(compuparams.debtnum, compuparams.resnum, compuparams.mnum, compuparams.ynum, 2)
		default=BitArray(compuparams.debtnum, compuparams.resnum, compuparams.mnum, compuparams.ynum, 2)
		new(reservesindefault, debt, reserves, default)
	end
end

immutable ReservesModel
	compuparams::ComputationParams
	econparams::EconParams
	grids::ModelGrids
	valuepay::Array{Float64,5}
	valuedefault::Array{Float64,3}
	bondprice::Array{Float64,4}
	policies::ModelPolicies
	cashinhandpay::Array{Float64,3}
	function ReservesModel(compuparams::ComputationParams, econparams::EconParams)
		grids::ModelGrids=makegrids(compuparams, econparams)
		valuepay=Array{Float64}(compuparams.debtnum, compuparams.resnum, compuparams.mnum, compuparams.ynum, 2)
		valuedefault=Array{Float64}(compuparams.resnum, compuparams.ynum, 2)
		bondprice=Array{Float64}(compuparams.debtnum, compuparams.resnum, compuparams.ynum, 2)
		policies=ModelPolicies(compuparams)
		cashinhandpay=broadcast( +, -(econparams.llambda+(1-econparams.llambda)*econparams.coupon)*grids.debt, 
								reshape( grids.reserves, 1, compuparams.resnum),
								reshape( grids.y, 1, 1, compuparams.ynum) )
		new(compuparams, econparams, grids, valuepay, valuedefault, bondprice, policies, cashinhandpay)
	end
end

immutable ModelSimulation
	periods::Int64
	debtind::Array{Int64,1}
	reservesind::Array{Int64,1}
	mind::Array{Int64,1}
	yind::Array{Int64,1}
    regime::Array{Int64,1}
	defaultstate::Array{Bool,1}
	randomshocks::Array{Float64,2}	
	output::Array{Float64,1}
	bondprice::Array{Float64,1}
	bondspread::Array{Float64,1}
	function ModelSimulation(per::Int64)
		periods=per
		debtind=Array{Int64}(periods+1)
		reservesind=Array{Int64}(periods+1)
		mind=Array{Int64}(periods)
		yind=Array{Int64}(periods)
	    regime=Array{Int64}(periods)
		defaultstate=Array{Bool}(periods+1)
		randomshocks=Array{Float64}(periods, 4)
		output=Array{Float64}(periods)
		bondprice=Array{Float64}(periods)	
		bondspread=Array{Float64}(periods)
		new(periods, debtind, reservesind, mind, yind, regime, defaultstate, randomshocks, output, bondprice, bondspread)	
	end
end


type ModelMoments
	# Mean indices
	debtmeanind::Int64
	reservesmeanind::Int64
	mmeanind::Int64
	ymeanind::Int64
	# Mean Values
	debtmean::Float64
	reservesmean::Float64
	outputmean::Float64
	spreadmean::Float64
	defaultstatemean::Float64
	defaultchoicemean::Float64
	# Variances
	outputsigma::Float64
	spreadsigma::Float64
	function ModelMoments()
		new()
	end
end

# Functions 
include("modelinitialize!.jl")
# Functions to solve model (maybe no need to upload every time)
include("valuepayexpectation!.jl")
include("solvedefaultvalue!.jl")
include("solvereservesmodel!.jl")
include("updatevaluepaydrm!.jl")	
include("mthresholds!.jl")
include("mergethresholds!.jl")
include("defaultthresholds!.jl")
include("suddenstopthresholds!.jl")
include("integratethresholds!.jl")
include("priceexpectation!.jl")
include("getpolicies!.jl")
# Simulation and Figures
include("simulatemodel!.jl")
include("getmoments!.jl")
include("moutputregimeexpectation!.jl")
end # Module end




