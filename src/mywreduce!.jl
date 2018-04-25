function mywreduce!(expected::AbstractArray{Float64,2}, # output
						current::Array{Float64,5}, mmass::Array{Float64,1}, ytrans::Array{Float64,1}, regimetrans::Array{Float64,1}) # inputs
	# Clean expected (allocated outside)
	fill!(expected,0.0)  # Must have corresponding reduced size
	# Expectation with respect to m-shock. Reduce dimension
	(debtnum,resnum,mnum,ynum,regimenum)=size(current)
	innersize::Int64=debtnum*resnum
	if length(expected)*mnum*ynum*regimenum!=length(current)
		error("Size mismatch")
		return
	end
	for iregime=1:regimenum
		for iy=1:ynum
			for imshock=1:mnum
				for ires=1:resnum
					for idebt=1:debtnum
						@inbounds expected[idebt, ires] += mmass[imshock]*ytrans[iy]*regimetrans[iregime]*current[idebt, ires, imshock, iy, iregime]
					end
				end
			end
		end
	end
	nothing
end