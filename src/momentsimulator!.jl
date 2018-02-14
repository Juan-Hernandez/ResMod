function momentsimulator!(compuparams::ComputationParams, econparams::EconParams, solverparams::SolverParams, outfilename::AbstractString)
	
	innermodel=ReservesModel(compuparams,econparams)
	modelinitialize!(innermodel)
	
	# 2.3. # Call solver routine

	(resiternum,valuegap,pricegap,defaultgap)=solvereservesmodel!(innermodel, solverparams)	
	
	# 3. Simulate model with seed=true
	
	innersimul=ModelSimulation(100000)
	simulatemodel!(innersimul,innermodel,true)
	
	# 4. Obtain moments
	innermoments=ModelMoments()
	flag=getmoments!(innermoments, innersimul, innermodel.grids, 1000) # burnin 1000 observations
	
	# 5. Return relevant moments
	calout=open(outfilename,"a")
		showcompact(calout, econparams.bbeta)
		print(calout, " |  ")
		showcompact(calout, econparams.defcost1)
		print(calout, " |  ")
		showcompact(calout, econparams.defcost2)
		print(calout, " |  ")
		showcompact(calout, econparams.panicfrequency)
		print(calout, " |  ")
		showcompact(calout, innermoments.debtmean)
		print(calout, "  | ")
		showcompact(calout, innermoments.reservesmean)
		print(calout, "  | ")
		showcompact(calout, innermoments.spreadmean)
		print(calout, " | ")
		showcompact(calout, innermoments.spreadsigma)
		print(calout, " | ")
		showcompact(calout, innermoments.defaultstatemean)
		print(calout, "  | ")
		showcompact(calout, innermoments.defaultchoicemean)
		print(calout, "  | ")
		showcompact(calout, innermoments.spreadXgdp)
		print(calout, "  | ")
		showcompact(calout, innermoments.spreadXgrowth)
		print(calout, "  | ")
		showcompact(calout, innermoments.deltaspreadXgdp)
		print(calout, "  | ")
		showcompact(calout, innermoments.deltaspreadXgrowth)
		print(calout, "  | ")
		showcompact(calout, maximum([valuegap,pricegap,defaultgap]) )
		println(calout, "  |")
	close(calout)
	return nothing
end	


