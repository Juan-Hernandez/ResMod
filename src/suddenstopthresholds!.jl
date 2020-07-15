function suddenstopthresholds!(thresholds::Array{Float64,1}, threspolicy::Array{Int64,2}, thresnum::Int64, thresdefault::BitArray{1}, # Outputs
								expvalue::Array{Float64,2}, valuedefault::Float64, defresoptim::Int64, consexm::Array{Float64,2}, debtmaxss::Int64,
								valtol::Float64, thrmin::Float64, ggamma::Int64, bbeta::Float64, mextremes::Array{Float64,1}, resnum::Int64, 
								interimnewthresholds::Array{Float64,1}, interimthrespolicy::Array{Int64,2}, interimthresnum::Int64)
	# Find thresholds for default if Sudden Stop (SS). THIS IS NOT THE SAME AS PUTTING AN UPPER BOUND ON DEBT!
	# This is important for solution. If the threshold for default if SS is mSS, then for all M's greater than mSS
	# the debt is unconstrained. This is, if m > mSS then V-(m) > X and unique Nash Equilibrium is lending which
	# leads to payoffs V+(m) and 0. Hence we must find the mSS threshold.
	# 
	
	
	# 0. Allocate temporary variables
	# 0.1 Temp Integers  
	toappend::Union{Nothing,Int64}=0
	# debtminfres::Int64=0
	# 0.2. Temp floats
	minthres::Float64=-thrmin   # Same as +MaxThreshold, always default
	tempthres::Float64=NaN
	vdiff::Float64=NaN
	# 0.3 Tempo booleans
	relevantss::Bool=true   # In principle Sudden Stop will be relevant
	# tempbitvec=BitArray(size(consexm,1))

	# 1. Compare all feasible debts against default. Get smallest threshold
	for idftres=1:resnum 	# Check all pairs of reserves and debt
		for idftdebt=1:debtmaxss # Will enter only if can consume if repay in SS


			
			if consexm[ idftdebt, idftres]>mextremes[1]     # Only enter for relevant future choice
				# Compare default with idftdebt, idftres
				vdiff=valuedefault-expvalue[idftdebt, idftres]
				if vdiff>-1e-6*(1.0-bbeta)*valtol # as difference in flow utility is negative, default is sure
					tempthres=-thrmin # Same as +MaxThreshold, always default
				else # Find threshold 
					# (this can be solved for all gamma)
					# tempthres=(vdiff*(1-ggamma)/(1.0-bbeta))^(1.0/(1-ggamma))-consexm[idftdebt, idftres]
					# Here we use gamma = 2
					tempthres=(bbeta-1.0)/vdiff-consexm[idftdebt, idftres]
					# Normalized utility
				end
				# Keep smallest m 
				minthres>tempthres && (minthres=tempthres)
				# Exit quickly
				minthres<mextremes[1] && break
			end
		end
	end
	# Note there is no need to know the debt
	# and reserves choice in this case, above
	# SSthresmin Vminus > X hence utility is V+
	if minthres<mextremes[1] # Repay for all M
		relevantss=false # Sudden Stop not relevant
	end
	# Enter only if Sudden Stop is Relevant
	if relevantss 
		# Again no default above threshold, append the rest. Default always below,
		@inbounds interimnewthresholds[1]=min(minthres, mextremes[end])
		@inbounds interimthrespolicy[1, 1]=1
		@inbounds interimthrespolicy[1, 2]=defresoptim
		interimthresnum=1
		@inbounds thresdefault[1]=true 	# Default on first interval. There can be no other default,
								# because default without SS implies default with SS
								# given that at this point continuation value of repayment is the same.
		# Append the other intervals, no default
		# toappend=findfirst(x->(x>interimnewthresholds[1]), thresholds[1:thresnum])
		toappend=findfirst(thresholds[1:thresnum].>interimnewthresholds[1])
		if toappend!=nothing
			@inbounds interimthresnum=thresnum-toappend+2
			@inbounds interimnewthresholds[2:interimthresnum]=thresholds[toappend:thresnum]
			@inbounds interimthrespolicy[2:interimthresnum, 1]=threspolicy[toappend:thresnum, 1]
			@inbounds interimthrespolicy[2:interimthresnum, 2]=threspolicy[toappend:thresnum, 2]
			@inbounds thresdefault[2:interimthresnum]=thresdefault[toappend:thresnum] # No default again    
		end
		thresnum=interimthresnum
		@inbounds thresholds[1:thresnum]=interimnewthresholds[1:thresnum]
		@inbounds threspolicy[1:thresnum, 1]=interimthrespolicy[1:thresnum, 1]
		@inbounds threspolicy[1:thresnum, 2]=interimthrespolicy[1:thresnum, 2]
	end
	return (thresnum, relevantss)
end # Function end