function momentsimulator!(compuparams::ComputationParams, econparams::EconParams, solverparams::SolverParams, outfilename::AbstractString)
	
innermodel=ReservesModel(compuparams,econparams)
modelinitialize!(innermodel)

# 2.3. # Call solver routine

(resiternum,valuegap,pricegap,defaultgap)=solvereservesmodel!(innermodel, solverparams, false)	
# 2.3.1 Solve again at lower update speed (if converged before, will exit quickly)
solverparams.updatespeed=solverparams.updatespeed*0.1
solverparams.startiternum=floor(Int64, resiternum)
solverparams.itermax=solverparams.itermax+800
solverparams.intermediatesave=solverparams.intermediatesave+800
solverparams.policiesout=false
(resiternum, valuegap, pricegap, defaultgap)=solvereservesmodel!(basemodel, solverparams, true)

# 3. Simulate model with seed=true

innersimul=ModelSimulation(101000)
simulatemodel!(innersimul,innermodel,true)

# 4. Obtain moments
innermoments=ModelMoments()
flag=getmoments!(innermoments, innersimul, innermodel.grids, 1000) # burnin 1000 observations

# 5. Return relevant moments
calout=open(outfilename,"a")
	iodef=IOContext(calout, :compact => true)
	show(iodef, econparams.bbeta)
	print(calout, " |  ")
	show(iodef, econparams.defcost1)
	print(calout, " |  ")
	show(iodef, econparams.defcost2)
	print(calout, " |  ")
	show(iodef, econparams.safeduration)
	print(calout, " |  ")
	show(iodef, econparams.riskduration)
	print(calout, " |  ")
	show(iodef, econparams.panicduration)
	print(calout, " |  ")
	show(iodef, econparams.panicfrequency)
	print(calout, " |  ")
	show(iodef, innermoments.debtmean)
	print(calout, "  | ")
	show(iodef, innermoments.reservesmean)
	print(calout, "  | ")
	show(iodef, innermoments.spreadmean)
	print(calout, " | ")
	show(iodef, innermoments.spreadsigma)
	print(calout, " | ")
	show(iodef, innermoments.defaultstatemean)
	print(calout, "  | ")
	show(iodef, innermoments.defaultchoicemean)
	print(calout, "  | ")
	show(iodef, innermoments.spreadXgdp)
	print(calout, "  | ")
	show(iodef, innermoments.spreadXgrowth)
	print(calout, "  | ")
	show(iodef, innermoments.deltaspreadXgdp)
	print(calout, "  | ")
	show(iodef, innermoments.deltaspreadXgrowth)
	print(calout, "  | ")
	show(iodef, max(valuegap,pricegap,defaultgap) )
	println(calout, "  |")
close(calout)
return nothing

end # Function end	


