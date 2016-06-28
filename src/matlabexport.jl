module matlabexport

using ReservesTypes
using JLD
using MAT 


@load "firstsolved.jld"

	file=matopen("juliaexport.mat","w")
	try
		write(file, "valuedefault", basemodel.valuedefault)
		write(file, "valuepay", basemodel.valuepay)
		write(file, "bondprice", basemodel.bondprice)
		write(file, "policies_debt", basemodel.policies.debt)
		write(file, "policies_reserves", basemodel.policies.reserves)
		write(file, "policies_reservesindefault", basemodel.policies.reservesindefault)
		write(file, "policies_default", bitunpack(basemodel.policies.default) )
		write(file, "grids", basemodel.grids)
		parnames=fieldnames(basemodel.econparams)
		strnames=map(string, parnames) 
		econdict=Dict{ASCIIString,Float64}(zip(strnames, map(getfield, repeated(basemodel.econparams,nfields(basemodel.econparams)), parnames)))
		write(file, "econparams", econdict)
		parnames=fieldnames(basemodel.compuparams)
		strnames=map(string, parnames) 
		compudict=Dict{ASCIIString,Float64}(zip(strnames, map(getfield, repeated(basemodel.compuparams,nfields(basemodel.compuparams)), parnames)))
		write(file, "compuparams", compudict)
	finally
		close(file)
	end

end # Module end