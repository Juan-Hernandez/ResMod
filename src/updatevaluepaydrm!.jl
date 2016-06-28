function updatevaluepaydrm!( newvaluepaydrm::AbstractArray{Float64,3}, newbondprice::AbstractArray{Float64,2}, 
							 debtpolicy::AbstractArray{Int64,3}, reservespolicy::AbstractArray{Int64,3}, defaultpolicy::AbstractArray{Bool,3}, #Outputs
							 expvalue::Array{Float64,2}, cashinhandpay::Array{Float64,2}, bondprice::Array{Float64,2},
							 valuedefault::Array{Float64,1}, defaultreservesoptim::Array{Int64,1}, regime::Int64, 
							 econparams::EconParams, compuparams::ComputationParams, grids::ModelGrids, policiesout::Bool)
"""	This function returns updated array of value (current debt,reserves and m), policies: next debt, 
	reserves and default (curr debt, res, m) for a fixed output and regime (output preallocation is done here),
	taking as input expected continuation value, cash in hand if pay, current bondprice matrix (future debt and reserves)
	default value and reserves choice, and regime. 
	New bond price will be an intermediate step for updating bond price. Average over mshok of future price times repay prob.
	the other exogenous expectations are done outside """

	# 3.0 Loops and preliminaries.
	# No need to preallocate outside since function will be parallelized
	consexm=Array{Float64}(compuparams.debtnum, compuparams.resnum)
	# This is the main thresholds vector policies and size
	thresholds=Array{Float64}(compuparams.debtnum*compuparams.resnum)
	threspolicy=Array{Int64}(compuparams.debtnum*compuparams.resnum, 2)
	thresnum::Int64=0
	thresdefault=falses(compuparams.debtnum*compuparams.resnum)
	# New threshold vector, policies and size to be merged (smaller for fixed reserves)
	newthresholds=Array{Float64}(compuparams.debtnum)
	newthrespolicy=Array{Int64}(compuparams.debtnum)
	newthresnum::Int64=0
	# Minimum future debt given future reserves
	debtminfres::Int64=0
	# incidence of sudden stop on thresholds
	relevantss::Bool=false

	# Interim allocations
	interimthresholds=Array{Float64}(compuparams.debtnum)
	interimthresdebt=Array{Int64}(compuparams.debtnum)

	interimbigthresholds=Array{Float64}(compuparams.debtnum*(compuparams.resnum+1))
	bigthresnum::Int64=0
	interimthrespolicy=Array{Int64}(compuparams.debtnum*compuparams.resnum, 2)
	interimthresnum::Int64=0
	mergedsortidx=Array{Int64}(compuparams.debtnum*(compuparams.resnum+1) ) 
	interimnewthresholds=Array{Float64}(compuparams.debtnum*(compuparams.resnum+1) )
	# small grid for values
	smallvalgrid=Array{Float64}(compuparams.mnum)
	smallpricegrid=Array{Float64}(compuparams.mnum)
	smallmassgrid=Array{Float64}(compuparams.mnum)
	# small grid for policies
	smallpolicygrid=Array{Int64}(compuparams.mnum, 2)
	smalldefaultgrid=falses(compuparams.mnum)

	cons0::Float64=0.0
	cons1::Float64=0.0
	vdiff::Float64=0.0
	tempthres::Float64=0.0
	debt1::Int64=0
	debt2::Int64=0
	res1::Int64=0
	res2::Int64=0
	tempindex::Int64=0

	for ires=1:compuparams.resnum
		# Consumption excluding M given future reserves.
		for idebt=1:compuparams.debtnum
			thresnum=0
			# Revenue from issuance			
			broadcast!(*, consexm, bondprice, grids.debt-(1-econparams.llambda)*grids.debt[idebt] )
			# Consumption excluding M given future reserves, future debt.
			broadcast!(+, consexm, cashinhandpay[idebt,ires], -grids.reserves'/(1+econparams.rfree), consexm)
			for idftres=1:compuparams.resnum
				debtminfres=findfirst(x->(x>grids.mextremes[1]), consexm[ :, idftres] )
				newthresnum=0
				if debtminfres!=0
					# 3.1 Find thresholds
					newthresnum=mthresholds!(newthresholds, newthrespolicy, # Outputs
									interimthresholds, interimthresdebt, # Interim allocations
									consexm[ :, idftres], expvalue[ :, idftres], compuparams,
									grids.mextremes, econparams.bbeta, debtminfres, 
									cons0, cons1, vdiff, tempthres )
					# 3.2 Merge with previous thresholds
					thresnum=mergethresholds!(thresholds, threspolicy, thresnum, # Outputs
										consexm, expvalue, econparams.bbeta, compuparams.valtol, grids.mextremes,
										newthresholds, newthrespolicy, newthresnum, idftres,
										interimbigthresholds, interimthrespolicy, bigthresnum,
										mergedsortidx, interimnewthresholds, interimthresnum,
										debt1, debt2, res1, res2, tempindex, cons0, cons1, vdiff, tempthres)
				end
			end
			# Here thresholds were merged for all future reserves.
			# 3.3 Enhance threshold with default decision
			(thresnum, debt1)=defaultthresholds!(thresholds, threspolicy, thresnum, thresdefault, # Outputs
								expvalue, valuedefault[ires], defaultreservesoptim[ires], consexm, compuparams.valtol,
								compuparams.thrmin, econparams.ggamma, econparams.bbeta, grids.mextremes[1],
								interimnewthresholds, interimthrespolicy, interimthresnum,
								vdiff, tempthres, cons0,
								debt1, res1)
			# 3.4 Check Sudden Stop impact
			relevantss=false
			if regime==2
				if (debt1!=0) && (debt1>grids.debtmaxss[idebt])
					# From defaultthresholds debt1 is the index of the first debt1 consistent with no default (lowest m)
					# if debt1=0 then default for all mshocks and no need to check SS
					(thresnum,relevantss)=suddenstopthresholds!(thresholds, threspolicy, thresnum, thresdefault, # Outputs	
                            	expvalue, valuedefault[ires], defaultreservesoptim[ires], consexm, grids.debtmaxss[idebt],
                            	compuparams.valtol, compuparams.thrmin, econparams.ggamma, econparams.bbeta, grids.mextremes,
								compuparams.resnum, interimnewthresholds, interimthrespolicy, interimthresnum,
								vdiff, tempthres, cons0, # cons0: token temporary Float64
                            	debtminfres, res1, relevantss)
				end
			end
			# 3.5 Integration
			integratethresholds!(smallvalgrid,	smallpricegrid, smallmassgrid,# Outputs
									thresholds, threspolicy, thresnum, thresdefault,
									grids.mextremes, grids.mmidpoints, bondprice, consexm, expvalue,
									valuedefault[ires], econparams, 
									vdiff, tempthres, cons0, cons1,
									res2, res1)
			newvaluepaydrm[idebt, ires, :]=smallvalgrid
			newbondprice[idebt, ires]=BLAS.dot(compuparams.mnum,grids.mmass,1,smallpricegrid,1)
			# 3.6 Get policies on grid.
			if policiesout
				getpolicies!( smallpolicygrid, smalldefaultgrid,
								thresholds, threspolicy, thresnum, thresdefault,
								grids.mmidpoints, compuparams.mnum)
				debtpolicy[idebt,ires,:]=smallpolicygrid[:,1]
				reservespolicy[idebt,ires,:]=smallpolicygrid[:,2]
				setindex!(defaultpolicy, smalldefaultgrid, idebt, ires, :)
			end
	 	end # Finish loop over current debt
	end # Finish loop over current reserves
end # Function end
