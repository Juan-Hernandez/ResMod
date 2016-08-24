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
	# 2.3 Fill output grids with zeros
	fill!(valmgrid, 0.0)
	fill!(pricegrid, 0.0)
	fill!(massgrid, 0.0)

	# 3. Loop
	for idthres=1:thresnum # Loop over M,
	    idmtop=findfirst( x->(x>=thresholds[idthres]) , mextremes)
	    idmlow=findfirst( x->(x>mstar), mextremes)
	    if (idthres==1) && (thresdefault[idthres]) # Default only for low M
	        valmgrid[idmtop-1]+=(thresholds[idthres]-mextremes[idmtop-1])/mdiff*valuedefault
	        mstar=thresholds[idthres]
	        massgrid[idmtop-1]+=(thresholds[idthres]-mextremes[idmtop-1])/mdiff
	        if idmtop>2 # Over all intervals before, default utility
	            valmgrid[1:(idmtop-2)]+=valuedefault
	            massgrid[1:(idmtop-2)]+=1
	        end
	    else # No default fill with care
	        currentprice=llambda+(1-llambda)*coupon+(1-llambda)*bondprice[ threspolicy[idthres, 1], threspolicy[idthres, 2] ]
	        # Recall price is an expectation over exogenous variables (and actions contingent in both exo and endo states)
	        # of (1-default)*(bond.service + remaining.bondprice). We begun with price on current exo states and future endogenous,
	        # seek to obtain price on yesterday exogenous and current endogenous, which requires current default and future 
	        # debt and reserves decisions. Variable currentprice is (bond.service+remaining.bondprice) for yesterdays exogenous and current 
	        # endogenous states.
	        for idj=idmlow:(idmtop-1)
	            currentvalue=consexm[ threspolicy[idthres, 1], threspolicy[idthres, 2] ] + 0.5*max(mextremes[idj-1],mstar)+0.5*mextremes[idj] # consumption
	            currentvalue=currentvalue^(1-ggamma)/(1-ggamma)*(1-bbeta) # flow utility
	           	currentvalue+=expvalue[ threspolicy[idthres, 1], threspolicy[idthres, 2] ] # plus continuation value
	            valmgrid[idj-1]+=min((mextremes[idj]-mstar)/mdiff, 1)*currentvalue
	            pricegrid[idj-1]+=min((mextremes[idj]-mstar)/mdiff, 1)*currentprice
	            massgrid[idj-1]+=min((mextremes[idj]-mstar)/mdiff, 1)
	        end
	        currentvalue=consexm[ threspolicy[idthres, 1], threspolicy[idthres, 2] ] + 0.5*(max(mextremes[idmtop-1], mstar)+thresholds[idthres])
	        currentvalue=currentvalue^(1-ggamma)/(1-ggamma)*(1-bbeta) # flow utility
	        currentvalue+=expvalue[ threspolicy[idthres, 1], threspolicy[idthres, 2] ]
	        valmgrid[idmtop-1]+=(thresholds[idthres]-max(mextremes[idmtop-1], mstar))/mdiff*currentvalue
	        pricegrid[idmtop-1]+=(thresholds[idthres]-max(mextremes[idmtop-1], mstar))/mdiff*currentprice
	        massgrid[idmtop-1]+=(thresholds[idthres]-max(mextremes[idmtop-1], mstar))/mdiff
	        mstar=thresholds[idthres]
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