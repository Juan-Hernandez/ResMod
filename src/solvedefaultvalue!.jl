function solvedefaultvalue!(model::ReservesModel, expectedvaluepay::Array{Float64,4}, defaultflowutility::Array{Float64,3}, 
                            newvaluedefault::Array{Float64,3}, interimdefaultvalue::Array{Float64,3}, 
                            reservesmaxtemp::Array{Float64,2}, reservesidtemp::Array{Int64,2}, resrestemp::Array{Float64,2})
    defaultgap::Float64=10*model.compuparams.valtol
    while defaultgap>(1-model.econparams.bbeta)*model.compuparams.valtol*0.0001
        # Expectations over future output and regime for default given current output, regime and FUTURE reserves.
        # The constant includes reentry prob
        fill!(interimdefaultvalue,0.0) 
        gemm!('N', 'T', (1-model.econparams.reentry)*model.econparams.bbeta*(1-model.econparams.norm2ss), 
        	sub(model.valuedefault, :, :, 1), model.grids.ytrans, 1.0, sub(interimdefaultvalue, :, :, 1) )
        gemm!('N', 'T', (1-model.econparams.reentry)*model.econparams.bbeta*(1-model.econparams.ss2ss), 
        	sub(model.valuedefault, :, :, 1), model.grids.ytrans, 1.0, sub(interimdefaultvalue, :, :, 2) )
        gemm!('N', 'T', (1-model.econparams.reentry)*model.econparams.bbeta*model.econparams.norm2ss, 
        	sub(model.valuedefault, :, :, 2), model.grids.ytrans, 1.0, sub(interimdefaultvalue, :, :, 1) )
        gemm!('N', 'T', (1-model.econparams.reentry)*model.econparams.bbeta*model.econparams.ss2ss, 
        	sub(model.valuedefault, :, :, 2), model.grids.ytrans, 1.0, sub(interimdefaultvalue, :, :, 2) )
		# update to continuation value, including reentry
        
		axpy!(model.econparams.reentry, slice(expectedvaluepay,1,:,:,1), sub( interimdefaultvalue, :, :, 1) )
        axpy!(model.econparams.reentry, slice(expectedvaluepay,1,:,:,2), sub( interimdefaultvalue, :, :, 2) )
		for ireg=1:2
			for iy=1:model.compuparams.ynum
				broadcast!(+, resrestemp, sub(defaultflowutility, :, :, iy), sub(interimdefaultvalue, :, iy,ireg) )
                #maximum!(sub(reservesmaxtemp , 1,:) , resrestemp)
                Base.findmax!(reservesmaxtemp , reservesidtemp, resrestemp)
                # Need to turn reservesindex from linear to cartesian, just care about first component
                broadcast!(mod1, reservesidtemp, reservesidtemp, model.compuparams.resnum)  
                setindex!(newvaluedefault, reservesmaxtemp, :, iy, ireg)
                setindex!(model.policies.reservesindefault, reservesidtemp, :, iy, ireg)                
            end
        end
                 
        axpy!(-1.0, newvaluedefault, model.valuedefault)
        maxabs!(sub(reservesmaxtemp, 1), model.valuedefault)
        defaultgap=reservesmaxtemp[1]        
        setindex!(model.valuedefault, newvaluedefault, :)                
    end
	return defaultgap
end

