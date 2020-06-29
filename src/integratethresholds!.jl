function integratethresholds!( valmgrid::Array{Float64,1},	pricegrid::Array{Float64,1}, massgrid::Array{Float64,1},# Outputs
								thresholds::Array{Float64,1}, threspolicy::Array{Int64,2}, thresnum::Int64, thresdefault::BitArray{1},
								mextremes::Array{Float64,1}, mmidpoints::Array{Float64,1}, bondprice::Array{Float64,2}, consexm::Array{Float64,2}, 
								expvalue::Array{Float64,2}, valuedefault::Float64, econparams::EconParams, debugbool::Bool)
	# 1. Unpack
	llambda::Float64=econparams.llambda
	coupon::Float64=econparams.coupon
	ggamma::Int64=econparams.ggamma
	bbeta::Float64=econparams.bbeta
	mextnum::Int64=size(mextremes,1)
	# 2. Initialize
	# 2.1 Temporary integer varaibles
	idmtop::Int64=2 # First interval [mextremes[1],mextremes[2]
	idmlow::Int64=1
	# 2.2 Floating point temporary vars
	mdiff::Float64=(mextremes[mextnum]-mextremes[1])/(mextnum-1)
	mstar::Float64=mextremes[1] # Current point of integration, will grow until mextremes[end]
	currentprice::Float64=0.0
	currentvalue::Float64=0.0
	thisthresdebt::Int64=0
	thisthresres::Int64=0
	# 2.3 Fill output grids with zeros
	fill!(valmgrid, 0.0)
	fill!(pricegrid, 0.0)
	fill!(massgrid, 0.0)

# 3. Loop over thresholds
	for idthres=1:thresnum # Loop over M,
		# idmtop: Index of first mextreme above threshold to integrate
		idmtop=min( mextnum-floor( Int64, (mextremes[mextnum]-thresholds[idthres])/mdiff ), mextnum)
	# 3.1 Deal with default (always only on first threshold)	
		if (idthres==1) && (thresdefault[idthres]) # Default only for low M
		# 3.1.1 Fill current interval: [idmtop-1] < threshold , [idmtop]
			@inbounds valmgrid[idmtop-1]+=(thresholds[idthres]-mextremes[idmtop-1])/mdiff*valuedefault
			@inbounds mstar=thresholds[idthres]
			debugbool && massgrid[idmtop-1]+=(thresholds[idthres]-mextremes[idmtop-1])/mdiff
		# 3.1.2 Fill all intervals before with default utility (no need to add anything to price)
			if idmtop>2 
				@inbounds valmgrid[1:(idmtop-2)].+=valuedefault
				debugbool && massgrid[1:(idmtop-2)].+=1.0			# 
			end
		else 
	# 3.2 Deal with non-default thresholds
			# idmlow: Index of first mextreme above current integration marker (mstar) 
			mstar<=mextremes[1] ? idmlow=2 : idmlow=ceil( Int64, (mstar-mextremes[1])/mdiff)+1
			# Check above needed for mstar=mextremes[1] at first pass (if no default)		 
			thisthresdebt=threspolicy[idthres, 1]
			thisthresres=threspolicy[idthres, 2]	    
			currentprice=llambda+coupon+(1.0-llambda)*bondprice[ thisthresdebt, thisthresres ]
			""" Recall price is an expectation over exogenous variables (and actions contingent in both exo and endo states)
				of (1-default)*(bond.service + remaining.bondprice). We began with price on current exo states and future endogenous,
				seek to obtain price on yesterday exogenous and current endogenous, which requires current default and future 
				debt and reserves decisions. Variable currentprice is (bond.service+remaining.bondprice) for yesterdays exogenous and current 
				endogenous states. """
		# 3.2.1 Fill all intervals from mstar < [idmlow] until [idmtop-1]		
			for idj=idmlow:(idmtop-1)
				currentvalue=consexm[ thisthresdebt, thisthresres ] + 0.5*max(mextremes[idj-1],mstar)+0.5*mextremes[idj] # consumption
				# currentvalue=currentvalue^(1-ggamma)/(1.0-ggamma)*(1.0-bbeta) # flow utility
				# using gamma = 2
				currentvalue=(bbeta-1.0)/currentvalue # flow utility
				@inbounds currentvalue+=expvalue[ thisthresdebt, thisthresres ] # plus continuation value
				@inbounds valmgrid[idj-1]+=min((mextremes[idj]-mstar)/mdiff, 1.0)*currentvalue
				@inbounds pricegrid[idj-1]+=min((mextremes[idj]-mstar)/mdiff, 1.0)*currentprice
				debugbool && massgrid[idj-1]+=min((mextremes[idj]-mstar)/mdiff, 1.0)
			end
		# 3.2.2 Fill current interval [idmtop-1] < threshold	
			currentvalue=consexm[ thisthresdebt, thisthresres ] + 0.5*(max(mextremes[idmtop-1], mstar)+thresholds[idthres])
			# currentvalue=currentvalue^(1-ggamma)/(1-ggamma)*(1.0-bbeta) # flow utility
			# using gamma = 2
			currentvalue=(bbeta-1.0)/currentvalue # flow utility
			   @inbounds currentvalue+=expvalue[ thisthresdebt, thisthresres ]
			   @inbounds valmgrid[idmtop-1]+=(thresholds[idthres]-max(mextremes[idmtop-1], mstar))/mdiff*currentvalue
			   @inbounds pricegrid[idmtop-1]+=(thresholds[idthres]-max(mextremes[idmtop-1], mstar))/mdiff*currentprice
			   debugbool && massgrid[idmtop-1]+=(thresholds[idthres]-max(mextremes[idmtop-1], mstar))/mdiff
			   @inbounds mstar=thresholds[idthres]
		end
	# 3.3 Check no excess mass integrated in any interval
		if debugbool && maximum(massgrid)>(1.0+1e-12)
			error("Relative mass on interval exceeds 1. Parameters: $econparams")
		end
	end
# 4. Check full mas was integrated
	if debugbool && minimum(massgrid)<(1.0-1e-12)
		error("Relative mass on interval below 1 after integration loop")
	end

end #function end
