function GGQthresholds!(thresholds::Array{Float64,1}, threspolicy::Array{Int64,2}, 								# Outputs
							consexm::Array{Float64,2}, expvalue::Array{Float64,2}, mextremes::Array{Float64,1},	# Array inputs
							bbeta::Float64, ggamma::Int64,	debtnum::Int64, resnum::Int64,						# Scalar inputs
                            valuesatm::Array{Float64,2})					# Temporary Array inputs 
                            
	mnew::Float64=0.0
	cons1::Float64=0.0
	cons2::Float64=0.0
	val1::Float64=-Inf
	val2::Float64=-Inf
	exp1::Float64=-Inf
	exp2::Float64=0.0

	maxvaluesatmisneg::Bool=true
	consstar::Float64=0.0
	thresnum::Int64=0
	# Begin at highest m
	mstar::Float64=mextremes[end]	

	# First one outside: find optimal choice at highest m
	for ires=1:resnum
		for idebt=1:debtnum
			@inbounds cons2=consexm[idebt, ires]
			@inbounds exp2=expvalue[idebt, ires]
			if cons2+mstar<1e-12
				@inbounds setindex!(valuesatm, -Inf, idebt,ires)
			else
				# This depends on gamma==2
				if (abs(exp2-exp1)<1.1e-14 && cons2>cons1) || ((cons2-cons1)>(exp1-exp2)*(cons1+mstar)*(cons2+mstar)/(1.0-bbeta)) # Need to improve speed
					cons1=cons2
					exp1=exp2
					@inbounds threspolicy[1,1]=idebt
					@inbounds threspolicy[1,2]=ires
				end
				# Reuse loop to set values at LOWEST m (needed for second round)
				if cons2+mextremes[1]>1e-12
					# Next line is very expensive:
					# @inbounds setindex!(valuesatm, (cons2+mextremes[1])^(1-ggamma)*(1.0-bbeta)/(1.0-ggamma)+exp2, idebt, ires) 
					# Since this function fails for gamma not equal to 2. we use the easy setup
					#############
					# This line has to change to a Brent's solver for ggamma not equal to 2
					@inbounds setindex!(valuesatm, (bbeta-1.0)/(cons2+mextremes[1])+exp2, idebt, ires) 
					maxvaluesatmisneg::Bool=false
				else
					@inbounds setindex!(valuesatm, -Inf, idebt, ires)
				end
			end
		end
	end	
	
	# If all consumptions at LOWEST m are negative return just first threshold
	maxvaluesatmisneg && (return 1)	# Could comment but faster if out quickly	

	while mstar>mextremes[1]+1e-10
		# Set threshold and number (policies set before but not counted if thresnum not increased)
		thresnum=thresnum+1
		@inbounds setindex!(thresholds, mstar, thresnum)
		# Set values to find the next mstar 
		@inbounds cons1=getindex(consexm, threspolicy[thresnum, 1], threspolicy[thresnum, 2] )
		@inbounds val1=getindex(valuesatm, threspolicy[thresnum, 1], threspolicy[thresnum, 2] )
		@inbounds exp1=getindex(expvalue, threspolicy[thresnum, 1] , threspolicy[thresnum, 2] )
		
		consstar=0
		mstar=mextremes[1]
		# Check all choices that beat the current choice at LOWEST m.
		# Since current is optimal at current threshold, there should be a lower threshold among current and candidate choice.
		# Also new candidate should have more consumption and less continuation value.
		for ires=1:resnum
			for idebt=1:debtnum # Entering a lot here. Perhaps use monotonicity to speed
				@inbounds val2=valuesatm[idebt, ires]
				@inbounds exp2=expvalue[idebt, ires]
				if val2>val1+1e-14 && exp2<exp1-1e-14 # Here it is important to have strict inequality for the -Inf cases
					@inbounds cons2=consexm[idebt, ires] # No need to asing outside
					#############
					# This line has to change to a Brent's solver for gamma not equal to 2
					mnew=-0.5*(cons1+cons2)+0.5*sqrt((cons2-cons1)^2+4.0*(cons2-cons1)/(exp1-exp2)*(1.0-bbeta))			
					#############
					if mstar<mnew+1e-9 || (abs(mstar-mnew)<2e-9 && cons2>consstar+1e-14)		# Need to improve speed
						# Keep the highest mnew
						mstar=mnew
						consstar=cons2 
						@inbounds threspolicy[thresnum+1, 1]=idebt
						@inbounds threspolicy[thresnum+1, 2]=ires
					end
				end
			end
		end
	end
	# Revert the vectors
	@inbounds thresholds[1:thresnum]=flipdim(thresholds[1:thresnum], 1)
	@inbounds threspolicy[1:thresnum, :]=flipdim(threspolicy[1:thresnum, :], 1)
	return thresnum
end
