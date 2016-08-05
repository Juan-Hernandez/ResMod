function GGQthresholds!(thresholds::Array{Float64,1}, threspolicy::Array{Int64,2}, 								# Outputs
							consexm::Array{Float64,2}, expvalue::Array{Float64,2}, m1::Float64,					# Array inputs
							bbeta::Float64, ggamma::Int64,	debtnum::Int64, resnum::Int64,						# Scalar inputs
                            valuesatm::Array{Float64,2})					# Temporary Array inputs 
                            
	mnew::Float64=0.0
	cons1::Float64=0.0
	val1::Float64=-Inf
	exp1::Float64=0.0

	thresnum::Int64=0
	# first one outside
	mstar::Float64=-m1	
	for idebt=1:debtnum
		for ires=1:resnum
			if consexm[idebt,ires]+mstar<0
				@inbounds setindex!(valuesatm, -Inf, idebt,ires)
			else
				@inbounds setindex!(valuesatm, (consexm[idebt, ires]+mstar)^(1-ggamma)*(1-bbeta)/(1-ggamma)+expvalue[idebt, ires], idebt, ires)
				if valuesatm[idebt, ires]>val1
					@inbounds threspolicy[1, :]=[idebt, ires]
					@inbounds val1=valuesatm[idebt, ires]
				end
				@inbounds setindex!(valuesatm, (max(1e-10, consexm[idebt, ires]+m1[1]) )^(1-ggamma)*(1-bbeta)/(1-ggamma)+expvalue[idebt, ires], idebt, ires)
			end
		end
	end
	
	
	while mstar>m1+1e-10
		# Set threshold and number (policies set before but not counted if thresnum not increased)
		thresnum=thresnum+1
		@inbounds setindex!(thresholds, mstar, thresnum)
		# Set values to find the next mstar 
		@inbounds cons1=getindex(consexm, threspolicy[thresnum, 1], threspolicy[thresnum, 2] )
		@inbounds val1=getindex(valuesatm, threspolicy[thresnum, 1], threspolicy[thresnum, 2] )
		@inbounds exp1=getindex(expvalue, threspolicy[thresnum, 1] , threspolicy[thresnum, 2] )
		
		mstar=m1
		for idebt=1:debtnum
			for ires=1:resnum 
				@inbounds if valuesatm[idebt, ires]>val1 #&& consexm[idebt, ifres]>cons1
					#############
					# This line has to change to a Brent's solver for gamma not equal to 2
					@inbounds mnew=-0.5*(cons1+consexm[idebt, ires])+0.5*((cons1-consexm[idebt, ires])^2+4*(cons1-consexm[idebt, ires])/(expvalue[idebt, ires]-exp1)*(1-bbeta))^0.5;			
					#############
					if mstar<mnew
						mstar=mnew 
						@inbounds threspolicy[thresnum+1, :]=[idebt, ires]
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
