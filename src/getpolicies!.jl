function getpolicies!(policygrid::Array{Int64,2}, defaultgrid::BitArray{1},
						thresholds::Array{Float64,1}, threspolicy::Array{Int64,2}, thresnum::Int64, thresdefault::BitArray{1},
						mmidpoints::Array{Float64,1}, mnum::Int64)
	idxthres::Int64=0
	for idm=1:mnum # Loop over M,
	    idxthres=findfirst( x->(x>=mmidpoints[idm]) , thresholds[1:thresnum])
	    policygrid[idm, :]=threspolicy[idxthres, :]
	    defaultgrid[idm]=thresdefault[idxthres]
	end
end