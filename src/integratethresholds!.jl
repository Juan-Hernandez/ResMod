function integratethresholds!( valmgrid::Array{Float64,1},	pricegrid::Array{Float64,1}, massgrid::Array{Float64,1},# Outputs
								thresholds::Array{Float64,1}, threspolicy::Array{Int64,2}, thresnum::Int64, thresdefault::BitArray{1},
								mextremes::Array{Float64,1}, mmidpoints::Array{Float64,1}, bondprice::Array{Float64,2}, consexm::Array{Float64,2}, 
								expvalue::Array{Float64,2}, valuedefault::Float64, econparams::EconParams)
	# 1. Unpack
	llambda::Float64=econparams.llambda
	coupon::Float64=econparams.coupon
	ggamma::Int64=econparams.ggamma
	bbeta::Float64=econparams.bbeta

	# 2. Initialize
	# 2.1 Temporary integer varaibles
	idmtop::Int64=2 # First interval [mextremes[1],mextremes[2]
	idmlow::Int64=1
	# 2.2 Floating point temporary vars
	mdiff::Float64=mextremes[2]-mextremes[1]
	mstar::Float64=mextremes[1] # Current point of integration, will grow until mextremes[end]
	currentprice::Float64=0.0
	currentvalue::Float64=0.0
	thisthresdebt::Int64=0
	thisthresres::Int64=0
	# 2.3 Fill output grids with zeros
	fill!(valmgrid, 0.0)
	fill!(pricegrid, 0.0)
	fill!(massgrid, 0.0)

	# 3. Loop
	for idthres=1:thresnum # Loop over M,
	    idmtop=findfirst(mextremes.>=thresholds[idthres])
	    idmlow=findfirst(mextremes.>mstar)
	    if (idthres==1) && (thresdefault[idthres]) # Default only for low M
	        @inbounds valmgrid[idmtop-1]+=(thresholds[idthres]-mextremes[idmtop-1])/mdiff*valuedefault
	        @inbounds mstar=thresholds[idthres]
	        @inbounds massgrid[idmtop-1]+=(thresholds[idthres]-mextremes[idmtop-1])/mdiff
	        if idmtop>2 # Over all intervals before, default utility
	            @inbounds valmgrid[1:(idmtop-2)]+=valuedefault
	            @inbounds massgrid[1:(idmtop-2)]+=1
	        end
	    else # No default fill with care
			@inbounds thisthresdebt=threspolicy[idthres, 1]
			@inbounds thisthresres=threspolicy[idthres, 2]	    
	        @inbounds currentprice=llambda+coupon+(1.0-llambda)*bondprice[ thisthresdebt, thisthresres ]
	        # Recall price is an expectation over exogenous variables (and actions contingent in both exo and endo states)
	        # of (1-default)*(bond.service + remaining.bondprice). We begun with price on current exo states and future endogenous,
	        # seek to obtain price on yesterday exogenous and current endogenous, which requires current default and future 
	        # debt and reserves decisions. Variable currentprice is (bond.service+remaining.bondprice) for yesterdays exogenous and current 
	        # endogenous states.
	        for idj=idmlow:(idmtop-1)
	            @inbounds currentvalue=consexm[ thisthresdebt, thisthresres ] + 0.5*max(mextremes[idj-1],mstar)+0.5*mextremes[idj] # consumption
	            currentvalue=currentvalue^(1-ggamma)/(1.0-ggamma)*(1.0-bbeta) # flow utility
	           	@inbounds currentvalue+=expvalue[ thisthresdebt, thisthresres ] # plus continuation value
	            @inbounds valmgrid[idj-1]+=min((mextremes[idj]-mstar)/mdiff, 1)*currentvalue
	            @inbounds pricegrid[idj-1]+=min((mextremes[idj]-mstar)/mdiff, 1)*currentprice
	            @inbounds massgrid[idj-1]+=min((mextremes[idj]-mstar)/mdiff, 1)
	        end
	        @inbounds currentvalue=consexm[ thisthresdebt, thisthresres ] + 0.5*(max(mextremes[idmtop-1], mstar)+thresholds[idthres])
	        currentvalue=currentvalue^(1-ggamma)/(1.0-ggamma)*(1.0-bbeta) # flow utility
	        @inbounds currentvalue+=expvalue[ thisthresdebt, thisthresres ]
	        @inbounds valmgrid[idmtop-1]+=(thresholds[idthres]-max(mextremes[idmtop-1], mstar))/mdiff*currentvalue
	        @inbounds pricegrid[idmtop-1]+=(thresholds[idthres]-max(mextremes[idmtop-1], mstar))/mdiff*currentprice
	        @inbounds massgrid[idmtop-1]+=(thresholds[idthres]-max(mextremes[idmtop-1], mstar))/mdiff
	        @inbounds mstar=thresholds[idthres]
	    end
	    # if maximum(pricegrid)>(Qrfree*(1+rfree)+0.01*ValTol)
	    #     error('Integrated value for q too big')
	    # end
	    if maximum(massgrid)>(1+1e-12)
	        error("Relative mass on interval exceeds 1")
	    end
	end

	if minimum(massgrid)<(1-1e-12)
	    error("Relative mass on interval below 1 after integration loop")
	end

end #function end