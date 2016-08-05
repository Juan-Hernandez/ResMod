function GGQthresholds!(thresholds::Array{Float64,1}, threspolicy::Array{Int64,2}, 								# Outputs
							consexm::Array{Float64,2}, expvalue::Array{Float64,2}, mextremes::Array{Float64,1},	# Array inputs
							bbeta::Float64, ggamma::Int64, debtnum::Int64,						# Scalar inputs
                            valuesatm::Array{Float64,2}, feasible::Array{Int64,2})					# Temporary Array inputs 
                            
	mnew::Float64=0.0
	cons1::Float64=0.0
	val1::Float64=0.0

	thresnum::Int64=0
	# first one outside
	mstar::Float64=mextremes[end]	
	setindex!(feasible,1:length(feasible),:)
	feasible[consexm+mstar.<0]=0
	fill!(valuesatm, -Inf)	
	valuesatm[feasible.>0]= (consexm[feasible.>0]+mstar).^(1-ggamma)*(1-bbeta)/(1-ggamma)+expvalue[feasible.>0]
	Base.findmax!(sub(thresholds, 1), sub(threspolicy, 1,1), valuesatm)	
	linpolicy::Int64=feasible[ threspolicy[1] ]
	# For reducing the feasible set:
	valuesatm[feasible.>0]= ( max(1e-10, consexm[feasible.>0]+mextremes[1]) ).^(1-ggamma)*(1-bbeta)/(1-ggamma)+expvalue[feasible.>0]	
	while mstar>mextremes[1]+1e-10
		# Set threshold and policies
		thresnum=thresnum+1
		# Recover linear policy index
		setindex!(threspolicy, [mod1(linpolicy, debtnum), div(linpolicy-1,debtnum)+1], thresnum, : )
		setindex!(thresholds, mstar, thresnum)
		# Reduce feasible set
		feasible[ valuesatm.<valuesatm[linpolicy] ]=0
		feasible[ consexm.<consexm[linpolicy] ]=0
		# 2.
		# Find new mstar
		cons1=consexm[linpolicy]
		val1=expvalue[linpolicy]
		# intialize candidate
		mstar=mextremes[1]
		for polind in feasible[feasible.>0]
			#############
			# This line has to change to a Brent's solver for gamma not equal to 2
			mnew=-0.5*(cons1+consexm[polind])+0.5*((cons1-consexm[polind])^2+4*(cons1-consexm[polind])/(expvalue[polind]-val1)*(1-bbeta))^0.5;			
			#############
			if mstar<mnew
				mstar=mnew 
				linpolicy=polind
			end
		end
		feasible[consexm+mstar.<0]=0
	end
	# Revert the vectors
	thresholds[1:thresnum]=flipdim(thresholds[1:thresnum], 1)
	threspolicy[1:thresnum, :]=flipdim(threspolicy[1:thresnum, :], 1)

	return thresnum
end
