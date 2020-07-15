function momentsimulator!(compuparams::ComputationParams, econparams::EconParams, solverparams::SolverParams, outfilename::AbstractString)
	
innermodel=ReservesModel(compuparams,econparams)
modelinitialize!(innermodel)

# 2.3. # Call solver routine

(resiternum,valuegap,pricegap,defaultgap)=solvereservesmodel!(innermodel, solverparams, false)	

# 3. Simulate model with seed=true

innersimul=ModelSimulation(101000)
simulatemodel!(innersimul,innermodel,true)

# 4. Obtain moments
innermoments=ModelMoments()
flag=getmoments!(innermoments, innersimul, innermodel.grids, 1000) # burnin 1000 observations

# 5. Return relevant moments
calout=open(outfilename,"a")
	show(IOContext(calout, :compact => true), econparams.bbeta)
	print(calout, " |  ")
	show(IOContext(calout, :compact => true), econparams.riskduration)
	print(calout, " |  ")
	show(IOContext(calout, :compact => true), econparams.defcost1)
	print(calout, " |  ")
	show(IOContext(calout, :compact => true), econparams.defcost2)
	print(calout, " |  ")
	show(IOContext(calout, :compact => true), econparams.panicfrequency)
	print(calout, " |  ")
	show(IOContext(calout, :compact => true), econparams.panicduration)
	print(calout, " |  ")
	show(IOContext(calout, :compact => true), econparams.llambda)
	print(calout, " |  ")
	show(IOContext(calout, :compact => true), innermoments.debtmean)
	print(calout, "  | ")
	show(IOContext(calout, :compact => true), innermoments.reservesmean)
	print(calout, "  | ")
	show(IOContext(calout, :compact => true), innermoments.spreadmean)
	print(calout, " | ")
	show(IOContext(calout, :compact => true), innermoments.spreadsigma)
	print(calout, " | ")
	show(IOContext(calout, :compact => true), innermoments.defaultstatemean)
	print(calout, "  | ")
	show(IOContext(calout, :compact => true), innermoments.defaultchoicemean)
	print(calout, "  | ")
	show(IOContext(calout, :compact => true), innermoments.spreadXgdp)
	print(calout, "  | ")
	show(IOContext(calout, :compact => true), innermoments.spreadXgrowth)
	print(calout, "  | ")
	show(IOContext(calout, :compact => true), innermoments.deltaspreadXgdp)
	print(calout, "  | ")
	show(IOContext(calout, :compact => true), innermoments.deltaspreadXgrowth)
	print(calout, "  | ")
	show(IOContext(calout, :compact => true), max(valuegap,pricegap,defaultgap) )
	println(calout, "  |")
close(calout)
return nothing

end # Function end	


