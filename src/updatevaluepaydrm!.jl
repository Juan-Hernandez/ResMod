function updatevaluepaydrm!( newvaluepaydrm::AbstractArray{Float64,3}, newbondprice::AbstractArray{Float64,2}, 
							 debtpolicy::AbstractArray{Int64,3}, reservespolicy::AbstractArray{Int64,3}, defaultpolicy::AbstractArray{Bool,3}, #Outputs
							 expvalue::Array{Float64,2}, cashinhandpay::Array{Float64,2}, bondprice::Array{Float64,2},
							 valuedefault::Array{Float64,1}, defaultreservesoptim::Array{Int64,1}, regime::Int64, valtol::Float64,
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
	# Minimum future debt given future reserves
	debtminfres::Int64=0
	# incidence of sudden stop on thresholds
	relevantss::Bool=false
	smallestMnodefdebt::Int64=0
	# Interim allocations
	valuesatm=Array{Float64}(compuparams.debtnum, compuparams.resnum)


	interimthresholds=Array{Float64}(compuparams.debtnum)
	interimthresdebt=Array{Int64}(compuparams.debtnum)


	interimthrespolicy=Array{Int64}(compuparams.debtnum*compuparams.resnum, 2)
	interimthresnum::Int64=0

	interimnewthresholds=Array{Float64}(compuparams.debtnum*(compuparams.resnum+1) )
	# small grid for values
	smallvalgrid=Array{Float64}(compuparams.mnum)
	smallpricegrid=Array{Float64}(compuparams.mnum)
	smallmassgrid=Array{Float64}(compuparams.mnum)
	# small grid for policies
	smallpolicygrid=Array{Int64}(compuparams.mnum, 2)
	smalldefaultgrid=falses(compuparams.mnum)
	



	for ires=1:compuparams.resnum
		# Consumption excluding M given future reserves.
		for idebt=1:compuparams.debtnum
			thresnum=0
			# Revenue from issuance			
			broadcast!(*, consexm, bondprice, grids.debt-(1-econparams.llambda)*grids.debt[idebt] )
			# Consumption excluding M given future reserves, future debt.
			broadcast!(+, consexm, cashinhandpay[idebt,ires], -grids.reserves'/(1+econparams.rfree), consexm)
			# 3.0  	# Check for positive consumption options, if not, default is sure and return only one threhold
		    if maximum(consexm)+grids.mextremes[end]<1e-13 # default is sure
        		thresnum=1
        		threspolicy[1,1]=1 								# Default debt
        		threspolicy[1,2]=defaultreservesoptim[ires] 	# Minimum reserves (will be changed to default reserves in defaulthresholds!.jl)
        		thresholds[1]=grids.mextremes[end] 				# Threshold is the biggest     
        		thresdefault[1]=true							# Always default
        		smallestMnodefdebt=0							# Default for all m-shocks
    		else
				# 3.1-2 GGQ algorithm.
				thresnum=GGQthresholds!(thresholds, threspolicy, 		# Outputs
								consexm, expvalue, grids.mextremes,		# Array inputs
								econparams.bbeta, econparams.ggamma, compuparams.debtnum, compuparams.resnum,	# Scalar inputs
	                            valuesatm)			# Temporary Array inputs 
				# Here thresholds were merged for all future reserves.
				# 3.3 Enhance threshold with default decision
				(thresnum, smallestMnodefdebt)=defaultthresholds!(thresholds, threspolicy, thresnum, thresdefault, # Outputs
									expvalue, valuedefault[ires], defaultreservesoptim[ires], consexm, valtol,
									compuparams.thrmin, econparams.ggamma, econparams.bbeta, grids.mextremes[1],
									interimnewthresholds, interimthrespolicy)
			end	
			# 3.4 Check Sudden Stop impact
			relevantss=false
			if regime==2
				if (smallestMnodefdebt!=0) && (smallestMnodefdebt>grids.debtmaxss[idebt])
					# From defaultthresholds smallestMnodefdebt is the index of the debt choice at the smallest M-shock 
					# consistent with no default. If smallestMnodefdebt=0 then default for all mshocks and no need to 
					# check SS. If it is positive but smaller than maxdebtss whenever repay is optimal it is still 
					# optimal without new lending, hence thresholds stay the same. 
					# Else, find new thresholds
					(thresnum,relevantss)=suddenstopthresholds!(thresholds, threspolicy, thresnum, thresdefault, # Outputs	
                            	expvalue, valuedefault[ires], defaultreservesoptim[ires], consexm, grids.debtmaxss[idebt],
                            	valtol, compuparams.thrmin, econparams.ggamma, econparams.bbeta, grids.mextremes,
								compuparams.resnum, interimnewthresholds, interimthrespolicy, interimthresnum)
				end
			end
			# 3.5 Integration
			integratethresholds!(smallvalgrid,	smallpricegrid, smallmassgrid,# Outputs
									thresholds, threspolicy, thresnum, thresdefault,
									grids.mextremes, grids.mmidpoints, bondprice, consexm, expvalue,
									valuedefault[ires], econparams)
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
