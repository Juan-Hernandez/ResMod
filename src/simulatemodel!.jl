function simulatemodel!(simulated::ModelSimulation, model::ReservesModel, fixseed::Bool=false)
	# 0. preallocation and initialization
	mtransbools=falses(model.compuparams.mnum)
	ytransbools=falses(model.compuparams.ynum)
	regtransbools=falses(model.econparams.regimenum)
	# Longindex to avoid several calls
	longstateindex::Int64=0

	# Initialize first period endogenous states:
	nextdebtind::Int64=cld(model.compuparams.debtnum, 2)
	nextreservesind::Int64=cld(model.compuparams.resnum, 2)
	nextdefaultstate::Bool=false
	
	# 1. Create cumulative transition probabilities
	cumytrans=Array{Float64}(undef, model.compuparams.ynum, model.compuparams.ynum)
	cumsum!(cumytrans, model.grids.ytrans, dims=2)
	cummmass=Array{Float64}(undef, model.compuparams.mnum)
	cumsum!(cummmass, model.grids.mmass)
	cumregtrans=Array{Float64}(undef, model.econparams.regimenum, model.econparams.regimenum)
	cumsum!(cumregtrans, model.grids.regimetrans, dims=2)

	# 2. Simulate all exogenous shocks: 1. y, 2. m, 3. regime, 4. reentry
	if fixseed
		fixrng=MersenneTwister(7010)
		rand!(fixrng, simulated.randomshocks)
	else
		rand!(simulated.randomshocks)
	end

	# 3-5. Big loop simulation
	for idper=1:simulated.periods
		# 3.0 Enter with previously chosen debt, reserves, and default state
		simulated.debtind[idper]=nextdebtind
		simulated.reservesind[idper]=nextreservesind
		simulated.defaultstate[idper]=nextdefaultstate  
		# 3. Find exogneous realizations
		# 3.1-2 Output and regime Markov transitions
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
		# 3.3 mshock iid realization
		mtransbools[:]=simulated.randomshocks[idper,2].<cummmass
		simulated.mind[idper]=findfirst(mtransbools)
		# 3.4 Reentry after default iid transition
		if simulated.defaultstate[idper] 
			simulated.defaultstate[idper]=(simulated.randomshocks[idper,4]>model.econparams.reentry) # Chance of reentry
		end

		# 4. Policies simulation
		# 4.0 Create current state index
		longstateindex=LinearIndices(model.policies.debt)[ simulated.debtind[idper], simulated.reservesind[idper], 
									simulated.mind[idper], simulated.yind[idper], simulated.regime[idper] ]
		
		# 4.1 Decide if default in current period ( if not already in default)
		if !simulated.defaultstate[idper]
			simulated.defaultstate[idper]=model.policies.default[longstateindex]
		end
		# 4.2 Policies in default state 
		# If not reentry or default in current period, still in default state
		if simulated.defaultstate[idper] 
			simulated.output[idper]=model.grids.ydefault[simulated.yind[idper]]+model.econparams.govtspend # mshock to zero
			nextdebtind=1 # Zero debt in default
			nextreservesind=model.policies.reservesindefault[ simulated.reservesind[idper], simulated.yind[idper], simulated.regime[idper] ]
			simulated.bondprice[idper]=0
			nextdefaultstate=true # Still in default state tomorrow ()
			simulated.bondspread[idper]=Inf
		else 
		# 4.3. Policies in repayment state: country has access to markets
			simulated.output[idper]=model.grids.y[simulated.yind[idper]]+model.grids.mmidpoints[simulated.mind[idper]]+model.econparams.govtspend		
			nextdebtind=model.policies.debt[longstateindex]
			nextreservesind=model.policies.reserves[longstateindex]
			nextdefaultstate=false 
			# Bond price today depends on debt and reserves tomorrow
			simulated.bondprice[idper]=model.bondprice[ nextdebtind, nextreservesind, simulated.yind[idper], simulated.regime[idper] ]
			simulated.bondspread[idper]=((model.econparams.coupon+model.econparams.llambda)/simulated.bondprice[idper]-
											model.econparams.llambda+1)^4-(1+model.econparams.rfree)^4
		end
	end
	nothing
end # Function end

