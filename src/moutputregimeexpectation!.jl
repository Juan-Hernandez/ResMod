function moutputregimeexpectation!(expected::Array{Float64,4}, current::Union{Array{Float64,5}, BitArray{5}}, 
									mmass::Array{Float64,1}, ytrans::Array{Float64,2}, regimetrans::Array{Float64,2},
									interimnom::Array{Float64,4})
	# Clean futuraveluepay (allocated outside)
	fill!(interimnom,0.0)  #Remember this has size 3 on last dimension
	# Expectation with respect to m-shock. Reduce dimension
	(debtnum,resnum,mnum,ynum,regimenum)=size(current)
	innersize::Int64=debtnum*resnum
	outersize::Int64=ynum*regimenum
	for iout=1:outersize
		for im=1:mnum
			for iin=1:innersize
				interimnom[(iout-1)*innersize+iin]+=mmass[im]*current[(iout-1)*mnum*innersize+(im-1)*innersize+iin]
			end
		end
	end
	# Expectation w.r.t. Regime and output
	fill!(expected,0.0)
	for ifutreg=1:regimenum
		for ifuty=1:ynum
			broadcast!( *, sub(interimnom, :, :, :, (regimenum+1)), interimnom[:,:,ifuty,ifutreg], reshape(ytrans[:,ifuty], 1, 1, ynum))
			for ireg=1:regimenum
				axpy!( regimetrans[ireg, ifutreg], interimnom[:,:,:,3], sub(expected, :, :, :, ireg) )
			end
		end
	end
	nothing
end