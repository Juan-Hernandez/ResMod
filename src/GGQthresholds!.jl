function GGQthresholds!(thresholds::Array{Float64,1}, threspolicy::Array{Int64,2}, 								# Outputs
							consexm::Array{Float64,2}, expvalue::Array{Float64,2}, mextremes::Array{Float64,1},	# Array inputs
							bbeta::Float64, ggamma::Int64, debtnum::Int64, resnum::Int64, debugbool::Bool,		# Scalar inputs
                            valuesatm::Array{Float64,2})					# Temporary Array inputs 
# 0 Initialize
    mnew::Float64=0.0
	cons1::Float64=0.0
	cons2::Float64=0.0
	val1::Float64=-Inf
	val2::Float64=-Inf
	exp1::Float64=-Inf
	expect2::Float64=0.0 # Do not overload exp2, this is a function in Base

	consdiff::Float64=0.0
	expdiff::Float64=0.0
	absexpdiff::Float64=0.0
	
	maxvaluesatmisneg::Bool=true
	consstar::Float64=0.0
	thresnum::Int64=0
	
# 1 Begin at highest m
	mstar::Float64=mextremes[end]	
	# First one outside: find optimal choice at highest m
	for ires=1:resnum
		for idebt=1:debtnum
			@inbounds cons2=consexm[idebt, ires]
			@inbounds expect2=expvalue[idebt, ires]
			if cons2+mstar<1e-12
				@inbounds setindex!(valuesatm, -Inf, idebt,ires)
			else
				# This depends on gamma==2
				consdiff=cons2-cons1
				expdiff=expect2-exp1
				if (expdiff>-1e-15 && consdiff>0.0) || (consdiff>-expdiff*(cons1+mstar)*(cons2+mstar)/(1.0-bbeta)) 
					cons1=cons2
					exp1=expect2
					@inbounds threspolicy[debtnum*resnum,1]=idebt
					@inbounds threspolicy[debtnum*resnum,2]=ires
				end
				# Reuse loop to set values at LOWEST m (needed for second round)
				if cons2+mextremes[1]>1e-12
					# Next line is very expensive:
					# # @inbounds 
					# setindex!(valuesatm, (cons2+mextremes[1])^(1-ggamma)*(1.0-bbeta)/(1.0-ggamma)+expect2, idebt, ires) 
					# Since this function fails for gamma not equal to 2. we use the easy setup
					#############
					# This line has to change to a Brent's solver for ggamma not equal to 2
					@inbounds setindex!(valuesatm, (bbeta-1.0)/(cons2+mextremes[1]+1e-12)+expect2, idebt, ires) 
					maxvaluesatmisneg=false
				else
					@inbounds setindex!(valuesatm, -Inf, idebt, ires)
				end
			end
		end
	end	

# 2 If all consumptions at LOWEST m are negative return just first threshold
	if maxvaluesatmisneg # Could comment but faster if out quickly	
		@inbounds setindex!(thresholds, mstar, 1)
		@inbounds setindex!(threspolicy, threspolicy[debtnum*resnum, 1:2], 1, 1:2)
		return 1
	end

# 3 Loop to find other thresholds
	while mstar>mextremes[1]+1e-12 # No threshold is less OR EQUAL to mextremes[1]
		# 3.1 Set threshold and number (policies set before but not counted if thresnum not increased)
		thresnum=thresnum+1
		@inbounds setindex!(thresholds, mstar, debtnum*resnum+1-thresnum)
		# 3.2 Get values to find the next mstar 
		@inbounds cons1=getindex(consexm, threspolicy[debtnum*resnum+1-thresnum, 1], threspolicy[debtnum*resnum+1-thresnum, 2] )
		@inbounds val1=getindex(valuesatm, threspolicy[debtnum*resnum+1-thresnum, 1], threspolicy[debtnum*resnum+1-thresnum, 2] )
		@inbounds exp1=getindex(expvalue, threspolicy[debtnum*resnum+1-thresnum, 1] , threspolicy[debtnum*resnum+1-thresnum, 2] )	
		consstar=0.0
		mstar=mextremes[1]+1e-12

		# 3.3 Check all choices that beat the current choice at LOWEST m.
		# Since current is optimal at current threshold, there should be a lower threshold among current and candidate choice.
		# Also new candidate should have more consumption and less continuation value.
		for ires=1:resnum
			for idebt=1:debtnum # Entering a lot here. Perhaps use monotonicity to speed
				# 3.3.0 Get value at mlow and continuation value for candidate
				@inbounds val2=valuesatm[idebt, ires]
				@inbounds expect2=expvalue[idebt, ires]
				expdiff=expect2-exp1
				# 3.3.1 Only enter if continuation value is lower and value at mlow is higher
				if val2>val1+1e-14 && expdiff<-1e-14 # Here it is important to have strict inequality for the -Inf cases
					@inbounds cons2=consexm[idebt, ires] # No need to assign outside
					consdiff=cons2-cons1
					# 3.3.2 Get breakeven threshold
					#############
					# WARNING #
					# This line has to change to a Brent's solver for gamma not equal to 2
					mnew=-0.5*(cons1+cons2)+0.5*sqrt(consdiff*consdiff-4.0*consdiff/expdiff*(1.0-bbeta))
					#############
					# 3.3.2.1 Check breakeven threshold is not bigger than current threshold
					if debugbool && mnew>thresholds[debtnum*resnum+1-thresnum] 
						display(threspolicy[(debtnum*resnum+1-thresnum):debtnum*resnum,:])
						println(thresholds[(debtnum*resnum+1-thresnum):debtnum*resnum])
						println([idebt,ires])
						error("mnew=", mnew, " . currthres=", thresholds[debtnum*resnum+1-thresnum], ". cons2=", cons2, " . consdiff=", consdiff, " . expdiff=", expdiff)
					end
					# 3.3.3 Keep highest breakeven threshold and policies, or similar with more consumption
					consdiff=mstar-mnew # use as tempholder for mstar-mnew
					if consdiff<0.0 || (consdiff<1e-9 && cons2>consstar+1e-14)	# Here consdiff is mstar-mnew.
						# Keep the highest mnew 
						mstar=mnew
						consstar=cons2 
						# Save at next thresnum policies, if mstar too close to border will be discarded
						@inbounds threspolicy[debtnum*resnum-thresnum, 1]=idebt
						@inbounds threspolicy[debtnum*resnum-thresnum, 2]=ires
					end
				end
			end
		end
	end
	# 4. Revert the vectors
	for idthres=1:thresnum
		@inbounds thresholds[idthres]=thresholds[debtnum*resnum-thresnum+idthres]
		@inbounds threspolicy[idthres, :]=threspolicy[debtnum*resnum-thresnum+idthres, :]
	end
	# 5. Return threshold number
	return thresnum
end
