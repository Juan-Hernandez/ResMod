function valuepayexpectation!(expectedvaluepay::Array{Float64,4},model::ReservesModel,interimvaluenom::Array{Float64,4})
	# Clean futuraveluepay (allocated outside)
	fill!(interimvaluenom,0.0)  #Remember this has size 3 on last dimension
	# Expectation with respect to m-shock. Reduce dimension
	innersize::Int64=model.compuparams.debtnum*model.compuparams.resnum
	outersize::Int64=model.compuparams.ynum*2
	for iout=1:outersize
		for im=1:model.compuparams.mnum
			for iin=1:innersize
				interimvaluenom[(iout-1)*innersize+iin]+=model.grids.mmass[im]*model.valuepay[(iout-1)*model.compuparams.mnum*innersize+(im-1)*innersize+iin]
			end
		end
	end
	# Expectation w.r.t. Regime and output
	fill!(expectedvaluepay,0.0)
	for iy=1:model.compuparams.ynum
		broadcast!( *, sub(interimvaluenom, :, :, :, 3), interimvaluenom[:,:,iy,1], reshape(model.grids.ytrans[:,iy], 1, 1, model.compuparams.ynum) )
		axpy!( model.econparams.bbeta*(1-model.econparams.norm2ss), interimvaluenom[:,:,:,3], sub(expectedvaluepay, :, :, :, 1) )
		axpy!( model.econparams.bbeta*(1-model.econparams.ss2ss), interimvaluenom[:,:,:,3], sub(expectedvaluepay, :, :, :, 2) )	
		broadcast!( *, sub(interimvaluenom, :, :, :, 3), interimvaluenom[:,:,iy,2], reshape(model.grids.ytrans[:,iy], 1, 1, model.compuparams.ynum) )
		axpy!( model.econparams.bbeta*model.econparams.norm2ss, interimvaluenom[:,:,:,3], sub(expectedvaluepay, :, :, :, 1) )
		axpy!( model.econparams.bbeta*model.econparams.ss2ss, interimvaluenom[:,:,:,3], sub(expectedvaluepay, :, :, :, 2) )
	end
	nothing
end