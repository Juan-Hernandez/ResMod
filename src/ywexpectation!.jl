function ywexpectation!(expected::Array{Float64,4}, current::DenseArray{Float64,4}, 
							ytrans::Array{Float64,2}, regimetrans::Array{Float64,2}, scalar::Float64,
							tempdry::Array{Float64,3})
	# Clean expected (allocated outside)
	fill!(expected,0.0)  # Same size as current
	(debtnum,resnum,ynum,regimenum)=size(current)
	# Expectation w.r.t. Regime and output
	if size(expected)!=size(current)
		error("Size mismatch")
		return
	end
	for ifutreg=1:regimenum
		for ifuty=1:ynum
			broadcast!( *, tempdry, current[:,:,ifuty,ifutreg], reshape(ytrans[:,ifuty], 1, 1, ynum))
			@inbounds for ireg=1:regimenum
				axpy!( scalar*regimetrans[ireg, ifutreg], tempdry, sub(expected, :, :, :, ireg) )
			end
		end
	end
	nothing
end