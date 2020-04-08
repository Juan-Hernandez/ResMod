function modelinitialize!(model::ReservesModel)
	# unpack:
	rfree=model.econparams.rfree
	llambda=model.econparams.llambda
	coupon=model.econparams.coupon
	yrho=model.econparams.logoutputrho*model.econparams.logoutputsigma^2/(model.econparams.logoutputsigma^2-model.compuparams.msigma^2)
	""" Prices depend on y, a', b': cannot add E[y'|y]+a' because of non
	linearity. It is not the same a high "a" but a low "y" than a low "a" but a
	high "y" because of persistance. """
	# To initialize price at risk free
	fill!(model.bondprice,(llambda+coupon)/(rfree+llambda))
	""" Initialize Value if default which depends only on Reserves, output and regime. 
	Assume constant consumption equivalent. This will be average of 
	current income and average default incoms plus return on reserves """
	# First reserves, then output, then mshock (recall mshock to ZERO when default), finally regime
	broadcast!( +, model.valuedefault, model.grids.reserves*rfree, 
				reshape( rfree/(1.0+rfree-yrho).*model.grids.ydefault, 1, model.compuparams.ynum), 
				(1.0-yrho)/(1.0+rfree-yrho)*dot(model.grids.yergodic,model.grids.ydefault),
				reshape([0.0,-0.001],1,1,2) )  					
	#	ValDefReg=ValDefReg.^(1-ggamma)./((1-ggamma)*(1-bbeta)); # Not Normalized
	model.valuedefault[:]=model.valuedefault[:].^(1-model.econparams.ggamma)./(1-model.econparams.ggamma) # Normalized Utility		
	""" Initialize value if not default which depends also on debt and Mshock.
	for debt we assume a consumption equivalent cost of coupon: rfree+mean(spread)=coupon """
	broadcast!( +, model.valuepay, -model.grids.debt*coupon, model.grids.reserves'*rfree,
    			reshape( model.grids.mmidpoints, 1, 1, model.compuparams.mnum ),
				reshape( rfree/(1+rfree-yrho).*model.grids.y, 1, 1, 1, model.compuparams.ynum )
				.+(1-yrho)/(1+rfree-yrho)*dot(model.grids.yergodic,model.grids.y),
    			reshape( [0.0,-0.005], 1, 1, 1, 1, 2) )
	model.valuepay[model.valuepay.<0].=0.0
	# ValPayAll=ValPayAll.^(1-ggamma)./((1-ggamma)*(1-bbeta)); # Not Normalized
	model.valuepay[:]=model.valuepay.^(1-model.econparams.ggamma)./(1-model.econparams.ggamma) #Normalized Utility	
	nothing
end
