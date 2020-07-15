function solvedefaultvalue!(model::ReservesModel, expectedvaluepay::Array{Float64,4}, defaultflowutility::Array{Float64,3}, 
                            newvaluedefault::Array{Float64,3}, interimdefaultvalue::Array{Float64,3}, 
                            reservesmaxtemp::Array{Float64,2}, reservesidtemp::Array{CartesianIndex{2},2}, resrestemp::Array{Float64,2},
                            tempry::Array{Float64,2}, valtol::Float64)
    
    # No while loop, solved in general iteration

    # defaultgap::Float64=10*valtol
    # while defaultgap>(1-model.econparams.bbeta)*valtol*0.001
        # Expectations over future output and regime for default given current output, regime and FUTURE reserves.
        # Reshaping to use ywexpectation
        fill!(interimdefaultvalue,0.0) 
        # The constant includes reentry prob
        ywexpectation!(reshape(interimdefaultvalue,(1,size(interimdefaultvalue)...)), # output matrix
                        reshape(model.valuedefault,(1,size(model.valuedefault)...)), # input matrix
                        model.grids.ytrans, model.grids.regimetrans, # transitions
                        (1.0-model.econparams.reentry)*model.econparams.bbeta, # scalar
                        reshape(tempry,(1,size(tempry)...))) # tempry

        # gemm!('N', 'T', (1-model.econparams.reentry)*model.econparams.bbeta*(1-model.econparams.norm2ss), 
        # 	view(model.valuedefault, :, :, 1), model.grids.ytrans, 1.0, view(interimdefaultvalue, :, :, 1) )
        # gemm!('N', 'T', (1-model.econparams.reentry)*model.econparams.bbeta*(1-model.econparams.ss2ss), 
        # 	view(model.valuedefault, :, :, 1), model.grids.ytrans, 1.0, view(interimdefaultvalue, :, :, 2) )
        # gemm!('N', 'T', (1-model.econparams.reentry)*model.econparams.bbeta*model.econparams.norm2ss, 
        # 	view(model.valuedefault, :, :, 2), model.grids.ytrans, 1.0, view(interimdefaultvalue, :, :, 1) )
        # gemm!('N', 'T', (1-model.econparams.reentry)*model.econparams.bbeta*model.econparams.ss2ss, 
        # 	view(model.valuedefault, :, :, 2), model.grids.ytrans, 1.0, view(interimdefaultvalue, :, :, 2) )

        regimenum=model.econparams.regimenum
        for iregime=1:regimenum
            # update to continuation value, including reentry
            axpy!(model.econparams.reentry, view(expectedvaluepay,1,:,:,iregime), view( interimdefaultvalue, :, :, iregime) )
			# Find optimal reserves for each y
            for iy=1:model.compuparams.ynum
				broadcast!(+, resrestemp, view(defaultflowutility, :, :, iy), view(interimdefaultvalue, :, iy,iregime) )
                Base.findmax!(reservesmaxtemp , reservesidtemp, resrestemp)
                # Need to extract optim reserves index from Cartesian on 2x2 matrix 
                setindex!(newvaluedefault, reservesmaxtemp, :, iy, iregime)
                for ires = 1:model.compuparams.resnum 
                    setindex!(model.policies.reservesindefault, reservesidtemp[ires][1], ires, iy, iregime)               
                end

            end
        end
        # Temporarily store gap of valuedefault in model.valuedefault         
        axpy!(-1.0, newvaluedefault, model.valuedefault)
        # Slow update speed for default not necesary. Not depending on price but on Value
        # when repaying, wich is updated slowly.
        defaultgap::Float64=maximum( abs.(model.valuedefault))   
        setindex!(model.valuedefault, newvaluedefault, :)                
    # end
	return defaultgap
end

