function mexpectation!(expected::Array{Float64,4}, # output
						current::Union{Array{Float64,5}, BitArray{5}}, mmass::Array{Float64,1}) # inputs
	# Clean expected (allocated outside)
	fill!(expected,0.0)  # Must have corresponding reduced size
	# Expectation with respect to m-shock. Reduce dimension
	(debtnum,resnum,mnum,ynum,regimenum)=size(current)
	innersize::Int64=debtnum*resnum
	outersize::Int64=ynum*regimenum
	if length(expected)*mnum!=length(current)
		error("Size mismatch")
		return
	end
	for iout=1:outersize
		for im=1:mnum
			for iin=1:innersize
				@inbounds expected[(iout-1)*innersize+iin]+=mmass[im]*current[(iout-1)*mnum*innersize+(im-1)*innersize+iin]
			end
		end
	end
	nothing
end