function findnewpriceCARA!(newbondprice::Array{Float64,4}, # Output
						bondcashflow::Array{Float64,5}, grids::ModelGrids, # Inputs
						grossrate::Float64, ggammalender::Float64, # Inputs
						expgcashflow::Array{Float64,5}, tempdryw::Array{Float64,4}, tempdry::Array{Float64,3}) # preallocated temps

# 1. Calculate stochastic discount factor at bond holdings
	# 1.1. Marginal utility at bond holdings cashflow: exp(-γ*b*cashflow)
	broadcast!( *, expgcashflow, -ggammalender, grids.debt, bondcashflow)
	expgcashflow=exp.(expgcashflow)
	# 1.2 Expected marginal utility of risk free bond: (1+r)*E[exp(-γ*b*cashflow)]
	mexpectation!(tempdryw, expgcashflow, grids.mmass)
	ywexpectation!(newbondprice, tempdryw, grids.ytrans, grids.regimetrans, grossrate, tempdry)

# 2. Calculate marginal utility of holding sovereign bond 
	# 2.1 Bond cahsflow times stochastic dicount factor: exp(-γ*b*cashflow)*cashflow 
	expgcashflow=expgcashflow.*bondcashflow
	# 2.2 Expected maginal utility of bond holding E[exp(-γ*b*cashflow)*cashflow]
	mexpectation!(tempdryw, expgcashflow, grids.mmass)
	# Keep the expected marginal utility of risk free bond that is stored in newbondprice
	setindex!(expgcashflow, newbondprice, : ,: ,1 , :, :)
	ywexpectation!(newbondprice, tempdryw, grids.ytrans, grids.regimetrans, 1.0, tempdry)

# 3. Calculate bondprice: q = E[exp(-γ*b*cashflow)*cashflow] / E[exp(-γ*b*cashflow)] / (1+r)
	newbondprice.=newbondprice./expgcashflow[:,:,1,:,:]
	return nothing
end
