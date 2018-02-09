function solvereservesmodel!(model::ReservesModel, solverparams=SolverParams)

	resiternum::Int64=solverparams.startiternum
	iterprint::Int64=solverparams.iterprint 
	itermax::Int64=solverparams.itermax 
	intermediatesave::Int64=solverparams.intermediatesave
	policiesout::Bool=solverparams.policiesout
	updatespeed::Float64=solverparams.updatespeed

	printbool::Bool=iterprint!=0

	# Unpack counters
	debtnum::Int64=model.compuparams.debtnum
	resnum::Int64=model.compuparams.resnum
	ynum::Int64=model.compuparams.ynum
	regimenum::Int64=size(model.grids.regimetrans,1)
	exonum::Int64=regimenum*ynum

	# Holder for new value functions
	newbondprice=SharedArray(Float64, debtnum, resnum, ynum, regimenum)
	newvaluedefault=Array{Float64}(resnum, ynum, regimenum)
	newvaluepay=SharedArray(Float64, debtnum, resnum, model.compuparams.mnum, ynum, regimenum)
	# Also need to allocate policies in shared arrays
	debtpolicy=SharedArray(Int64, debtnum, resnum, model.compuparams.mnum, ynum, regimenum)
	reservespolicy=SharedArray(Int64, debtnum, resnum, model.compuparams.mnum, ynum, regimenum)
	defaultpolicy=SharedArray(Bool, debtnum, resnum, model.compuparams.mnum, ynum, regimenum)

	# initial values for error
	valuegap::Float64=10.0*solverparams.valtol
	pricegap::Float64=10.0*solverparams.valtol
	defaultgap::Float64=10.0*solverparams.valtol
	# Holder for exogenous expectation
	expectedvaluepay=Array{Float64}(debtnum,resnum,ynum,regimenum)
	# Preallocate intermediate stages

	# Flow utility matrix after default and chossing FUTURE reserves.
	# (futurereserves, currentreserves, currentoutput)
	defaultflowutility=Array{Float64}(resnum,resnum,ynum)
	broadcast!(+, defaultflowutility, -model.grids.reserves./(1+model.econparams.rfree), model.grids.reserves', reshape(model.grids.ydefault,1,1,ynum))
	defaultflowutility[defaultflowutility.<0.0]=0.0
	defaultflowutility=defaultflowutility.^(1-model.econparams.ggamma)./(1-model.econparams.ggamma).*(1-model.econparams.bbeta)
	# Preallocation of temporary arrays
	tempdry=Array{Float64}(debtnum,resnum,ynum)
	tempdryw=Array{Float64}(debtnum,resnum,ynum,regimenum)
	
	tempry=Array{Float64}(resnum,ynum) 
	temprr=Array{Float64}(resnum, resnum)
	tempryw=Array{Float64}(resnum,ynum,regimenum) 
	
    reservesmaxtemp=Array{Float64}(1, resnum)
    reservesidtemp=Array{Int64}(1, resnum)
    

    printbool && println("        valuegap         |        pricegap         |   iternum  | time (sec) |")

    # Main loop
    tic()
	while resiternum<itermax && ( max(valuegap, pricegap, defaultgap)>=solverparams.valtol || !policiesout )
		# 0. Output and update controls
		# Increase iteration counter
		resiternum+=1
		# Set intermediate save
		if mod1(resiternum,intermediatesave)==intermediatesave
			jldopen("debugoldmodel.jld", "w") do file
				write(file, "oldmodel", model)
				write(file, "valuegap", valuegap)
				write(file, "pricegap", 100*pricegap)  
				write(file, "iternum", resiternum-1)
			end
		end
		# Intermediate print and timing
		printbool && ( (mod1(resiternum,iterprint)==1) && tic() )
		# Policies out in last iteration
		(max(valuegap, pricegap, defaultgap)>=solverparams.valtol) && (policiesout=true)
		(resiternum==itermax-1) && (policiesout=true)

		# 1. Expectation on exogenous varaibles
		mexpectation!(tempdryw, model.valuepay, model.grids.mmass)
		ywexpectation!(expectedvaluepay, tempdryw, # here tempdry has the expectation over mshock 
							model.grids.ytrans, model.grids.regimetrans,model.econparams.bbeta,tempdry)
		
		# 2. update value of default
		defaultgap=solvedefaultvalue!(model, expectedvaluepay, defaultflowutility, newvaluedefault, 
								tempryw, reservesmaxtemp, reservesidtemp, temprr, tempry, solverparams.valtol)
		# 3. Update value function: parallel for each set of exogenous states. 
			# Here i need to use comprehensions to pass values (no need to predefine oustide the loop)
			# pmap is enough since function return is inplaced
			# Careful, use sub() for those arrays you want to be referenced. [:,:,:] is good for those that can be copied.
			# Since we are pmapping, 		
""" For prototyping we use just one evaluation, pmap below"""		
		# updatevaluepaydrm!( sub(newvaluepay, :, :, :, 3), sub(newbondprice, :, :, 3), 
		# 		sub(debtpolicy, :, :, :, 3), sub(reservespolicy, :, :, :, 3), sub(defaultpolicy, :, :, :, 3), # Outputs
		# 		expectedvaluepay[ :, :, 3, 1], model.cashinhandpay[ :, :, 3], model.bondprice[:,:,3,1],
		# 		newvaluedefault[ :, 3], model.policies.reservesindefault[:, 3], 1, model.econparams, model.compuparams, model.grids, true )		
		# pmap requires shared arrays for inplace! outputs
		pmap(updatevaluepaydrm!, [ sub(newvaluepay, :, :, :, iyr) for iyr=1:exonum], [sub(newbondprice, :, :, iyr) for iyr=1:exonum], # Outputs
			[sub(debtpolicy, :, :, :, iyr) for iyr=1:exonum], [sub(reservespolicy, :, :, :, iyr) for iyr=1:exonum], # Outputs
			[sub(defaultpolicy, :, :, :, iyr) for iyr=1:exonum],  # Outputs
			[expectedvaluepay[ :, :, iyr] for iyr=1:exonum], [model.cashinhandpay[ :, :, mod1(iyr,ynum)] for iyr=1:exonum],
			[model.bondprice[ :, :, iyr] for iyr=1:exonum], [newvaluedefault[ :, iyr] for iyr=1:exonum], 
			[model.policies.reservesindefault[:, iyr] for iyr=1:exonum], [cld(iyr, ynum) for iyr=1:exonum], 
			repeated( solverparams.valtol, exonum), repeated( model.econparams, exonum), repeated( model.compuparams, exonum), 
			repeated( model.grids, exonum), repeated(policiesout, exonum) )
		# 4. Find new price: take expectation over regime and output
		ywexpectation!(tempdryw, newbondprice, 
							model.grids.ytrans, model.grids.regimetrans, 1.0/(1.0+model.econparams.rfree),
							tempdry)
		setindex!(newbondprice, tempdryw, :)
		# 5. Find gaps and update 
		axpy!(-1.0, newvaluepay, model.valuepay)
		maxabs!(sub(reservesmaxtemp,1:1), model.valuepay)
		valuegap=reservesmaxtemp[1]/(1-model.econparams.bbeta) # Because for higher beta changes in V are more meaningful
		setindex!(model.valuepay, newvaluepay, :)
		# Cannot do the same, old price cannot overwritten because of updatespeed. Recall tempdryw also holds newbondprice
		axpy!(-1.0, model.bondprice, tempdryw )
		maxabs!(sub(reservesmaxtemp,1:1), tempdryw)
		pricegap=reservesmaxtemp[1]
		# Update Control
		scal!(debtnum*resnum*ynum*regimenum, 1-updatespeed, model.bondprice, 1 )
		axpy!(updatespeed, newbondprice, model.bondprice )
		# update policies
	    if policiesout
			setindex!(model.policies.debt, debtpolicy, :)
			setindex!(model.policies.reserves, reservespolicy, :)
			setindex!(model.policies.default, defaultpolicy, :)
		end		
    	# 5.1 Print intermediate output
    	if printbool && mod1(resiternum,iterprint)==iterprint
	    	print("       ")
		    showcompact(valuegap)
	    	print("        |       ")
	    	showcompact(100*pricegap) 
	    	print("       |     ")
	    	showcompact(resiternum)
	    	print("    | ")
	    	showcompact(toq())
	    	println("  |")
	    	
	    end
	    # 5.2 Save intermediate step
		if mod1(resiternum,intermediatesave)==intermediatesave
			jldopen("debugmodels.jld", "w") do file
				write(file, "basemodel", model)
				write(file, "valuegap", valuegap)
				write(file, "pricegap", 100*pricegap)  
				write(file, "iternum", resiternum)
			end
		end
		# # Just to exit quickly
		# valuegap=0
		# pricegap=0
	end

	if printbool && resiternum%iterprint!=0 # Print last iteration
    	print("       ")
	    showcompact(valuegap)
    	print("        |       ")
    	showcompact(100*pricegap) 
    	print("       |     ")
    	showcompact(resiternum)
    	print("    |   ")
    	showcompact(toq())
    	println("  |")
	end
	return (resiternum,valuegap,pricegap,defaultgap)
end # Function End

solvereservesmodel!(model::ReservesModel)=solvereservesmodel!(model, SolverParams( ))
#solvereservesmodel!(model::ReservesModel, x::Any)=solvereservesmodel!(model, SolverParams(x...))


