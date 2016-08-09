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

	consstar::Float64=0.0
	thresnum::Int64=0
	# first one outside
	mstar::Float64=mextremes[end]	
	for ires=1:resnum
		for idebt=1:debtnum
			@inbounds cons2=consexm[idebt,ires]
			@inbounds exp2=expvalue[idebt, ires]
			if cons2+mstar<1e-12
				@inbounds setindex!(valuesatm, -Inf, idebt,ires)
			else
				# This depends on gamma==2
				if (abs(exp2-exp1)<1.1e-14 && cons2>cons1) || ((cons2-cons1)>(exp1-exp2)*(cons1+mstar)*(cons2+mstar)/(1-bbeta))
					cons1=cons2
					exp1=exp2
					@inbounds threspolicy[1, :]=[idebt, ires]
				end
				if cons2+mextremes[1]>1e-12
					@inbounds setindex!(valuesatm, (cons2+mextremes[1])^(1-ggamma)*(1-bbeta)/(1-ggamma)+exp2, idebt, ires)
				else
					@inbounds setindex!(valuesatm, -Inf, idebt, ires)
				end
			end
		end
	end
	
	
	while mstar>mextremes[1]+1e-8
		# Set threshold and number (policies set before but not counted if thresnum not increased)
		thresnum=thresnum+1
		@inbounds setindex!(thresholds, mstar, thresnum)
		# Set values to find the next mstar 
		@inbounds cons1=getindex(consexm, threspolicy[thresnum, 1], threspolicy[thresnum, 2] )
		@inbounds val1=getindex(valuesatm, threspolicy[thresnum, 1], threspolicy[thresnum, 2] )
		@inbounds exp1=getindex(expvalue, threspolicy[thresnum, 1] , threspolicy[thresnum, 2] )
		
		consstar=0
		mstar=mextremes[1]
		for ires=1:resnum
			for idebt=1:debtnum
				@inbounds cons2=consexm[idebt, ires]
				@inbounds val2=valuesatm[idebt, ires]
				@inbounds exp2=expvalue[idebt, ires]
				if val2>=val1+1e-14 && exp2<exp1-1e-14 # Here consumption can only increase
					#############
					# This line has to change to a Brent's solver for gamma not equal to 2
					mnew=-0.5*(cons1+cons2)+0.5*((cons2-cons1)^2+4*(cons2-cons1)/(exp1-exp2)*(1-bbeta))^0.5;			
					#############
					if mstar<mnew+1e-8 && 
					mstar<mnew-1e-8 || (abs(mstar-mnew)<1.1e-8 && cons2>consstar+1e-14)
						mstar=mnew
						consstar=cons2 
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
