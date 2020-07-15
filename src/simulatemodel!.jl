function simulatemodel!(simulated::ModelSimulation, model::ReservesModel, fixseed::Bool=false)
	# 0. preallocation and initialization
	mtransbools=falses(model.compuparams.mnum)
	ytransbools=falses(model.compuparams.ynum)
	regtransbools=falses(model.econparams.regimenum)

	# Initialize first period endogenous states:
	simulated.debtind[1]=cld(model.compuparams.debtnum, 2)
	simulated.reservesind[1]=cld(model.compuparams.resnum, 2)
	simulated.defaultstate[1]=false
	# Simulate all exogenous shocks: 1. y, 2. m, 3. regime, 4. reentry
	if fixseed
		fixrng=MersenneTwister(7010)
		rand!(fixrng, simulated.randomshocks)
	else
		rand!(simulated.randomshocks)
	end
	cumytrans=Array{Float64}(undef, model.compuparams.ynum, model.compuparams.ynum)
	cumsum!(cumytrans, model.grids.ytrans, dims=2)
	cummmass=Array{Float64}(undef, model.compuparams.mnum)
	cumsum!(cummmass, model.grids.mmass)
	cumregtrans=Array{Float64}(undef, model.econparams.regimenum, model.econparams.regimenum)
	cumsum!(cumregtrans, model.grids.regimetrans, dims=2)

	# Longindex to avoid several calls
	longstateindex::Int64=0

	# 1. Simulation
	for idper=1:simulated.periods
	   # Enter with previously chosen debt, reserves, and default state  
		if idper==1 # Since Markov exogenous procesess assume initial value, t=0 output was median index, t=0 regime was first ss
			ytransbools[:]=simulated.randomshocks[idper,1].>cumytrans[cld(model.compuparams.ynum, 2), :]
			simulated.yind[idper]=findfirst(ytransbools)
			regtransbools[:]= simulated.randomshocks[idper, 3].<cumregtrans[1,:]
			simulated.regime[idper]=findfirst(regtransbools)
		else
			ytransbools[:]=simulated.randomshocks[idper,1].<cumytrans[simulated.yind[idper-1], :]
			simulated.yind[idper]=findfirst(ytransbools)
			regtransbools[:]= simulated.randomshocks[idper, 3].<cumregtrans[simulated.regime[idper-1], :] 
			simulated.regime[idper]=findfirst(regtransbools)
		end
		mtransbools[:]=simulated.randomshocks[idper,2].<cummmass
		simulated.mind[idper]=findfirst(mtransbools)
		# Reentry
		if simulated.defaultstate[idper] 
			simulated.defaultstate[idper]=(simulated.randomshocks[idper,4]>model.econparams.reentry) # Chance of reentry
		end
		# Default in current period
		longstateindex=LinearIndices(model.policies.debt)[ simulated.debtind[idper], simulated.reservesind[idper], 
									simulated.mind[idper], simulated.yind[idper], simulated.regime[idper] ]
		if !simulated.defaultstate[idper]
			simulated.defaultstate[idper]=model.policies.default[longstateindex]
		end
		# if not reentry or default in current period, still in default state
		if simulated.defaultstate[idper] 
			simulated.output[idper]=model.grids.ydefault[simulated.yind[idper]]+model.econparams.govtspend # mshock to zero
			simulated.debtind[idper+1]=1 # Zero debt in default
			simulated.reservesind[idper+1]=model.policies.reservesindefault[ simulated.reservesind[idper], simulated.yind[idper], simulated.regime[idper] ]
			simulated.bondprice[idper]=0
			simulated.defaultstate[idper+1]=true # Still in default state tomorrow ()
			simulated.bondspread[idper]=Inf
		else 
			# below, the country has acces to markets
			simulated.output[idper]=model.grids.y[simulated.yind[idper]]+model.grids.mmidpoints[simulated.mind[idper]]+model.econparams.govtspend		
			simulated.debtind[idper+1]=model.policies.debt[longstateindex]
			simulated.reservesind[idper+1]=model.policies.reserves[longstateindex]
			simulated.defaultstate[idper+1]=false 
			# Bond price today depends on debt and reserves tomorrow
			simulated.bondprice[idper]=model.bondprice[ simulated.debtind[idper+1], simulated.reservesind[idper+1], simulated.yind[idper], simulated.regime[idper] ]
			simulated.bondspread[idper]=((model.econparams.coupon+model.econparams.llambda)/simulated.bondprice[idper]-
											model.econparams.llambda+1)^4-(1+model.econparams.rfree)^4
		end
	end
	nothing
end # Function end

