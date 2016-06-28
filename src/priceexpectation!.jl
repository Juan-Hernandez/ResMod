function priceexpectation!(newbondprice::Union{Array{Float64,4},SharedArray{Float64,4}}, interimprice::Array{Float64,4},
							econparams::EconParams, ytrans::Array{Float64,2}, ynum::Int64)
	# Clean futuraveluepay (allocated outside)
	fill!(interimprice,0.0)  #Remember this has size 3 on last dimension
	# Expectation w.r.t. Regime and output
	for iy=1:ynum
		broadcast!( *, sub(interimprice, :, :, :, 3), newbondprice[:,:,iy,1], reshape(ytrans[:,iy], 1, 1, ynum) )
		axpy!( (1-econparams.norm2ss)/(1+econparams.rfree), interimprice[:,:,:,3], sub(interimprice, :, :, :, 1) )
		axpy!( (1-econparams.ss2ss)/(1+econparams.rfree), interimprice[:,:,:,3], sub(interimprice, :, :, :, 2) )	
		broadcast!( *, sub(interimprice, :, :, :, 3), newbondprice[:,:,iy,2], reshape(ytrans[:,iy], 1, 1, ynum) )
		axpy!( econparams.norm2ss/(1+econparams.rfree), interimprice[:,:,:,3], sub(interimprice, :, :, :, 1) )
		axpy!( econparams.ss2ss/(1+econparams.rfree), interimprice[:,:,:,3], sub(interimprice, :, :, :, 2) )
	end
	setindex!(newbondprice, interimprice[:,:,:,1:2], :)
	nothing
end