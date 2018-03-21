function getpolicies!(policygrid::Array{Int64,2}, defaultgrid::BitArray{1},
						thresholds::Array{Float64,1}, threspolicy::Array{Int64,2}, thresnum::Int64, thresdefault::BitArray{1},
						mmidpoints::Array{Float64,1}, mnum::Int64)
	# For speed need to write in similar fashion as integratethresholds
	# Not critical since run only once
	idxthres::Int64=0
	for idm=1:mnum # Loop over M,
	    idxthres=findfirst(thresholds[1:thresnum].>=mmidpoints[idm])
	    @inbounds policygrid[idm, 1:2]=threspolicy[idxthres, 1:2]
	    defaultgrid[idm]=thresdefault[idxthres]
	end
end