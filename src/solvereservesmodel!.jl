function solvereservesmodel!(model::ReservesModel, policiesout::Bool=false, 
								resiternum::Int64=0, itermax::Int64=6001, iterprint::Int64=10, intermediatesave::Int64=1000)
	# Unpack counters
	debtnum::Int64=model.compuparams.debtnum
	resnum::Int64=model.compuparams.resnum
	ynum::Int64=model.compuparams.ynum

	# Holder for new value functions
	newbondprice=SharedArray(Float64, debtnum, resnum, ynum, 2)
	newvaluedefault=Array{Float64}(resnum, ynum, 2)
	newvaluepay=SharedArray(Float64, debtnum, resnum, model.compuparams.mnum, ynum, 2)
	# Also need to allocate policies in shared arrays
	debtpolicy=SharedArray(Int64, debtnum, resnum, model.compuparams.mnum, ynum, 2)
	reservespolicy=SharedArray(Int64, debtnum, resnum, model.compuparams.mnum, ynum, 2)
	defaultpolicy=SharedArray(Bool, debtnum, resnum, model.compuparams.mnum, ynum, 2)

	# initial values for error
	valuegap::Float64=10.0*model.compuparams.valtol
	pricegap::Float64=10.0*model.compuparams.valtol
	defaultgap::Float64=10.0*model.compuparams.valtol
	# Holder for exogenous expectation
	expectedvaluepay=Array{Float64}(debtnum,resnum,ynum,2)
	# Preallocate intermediate stages


	# Flow utility matrix after default and chossing FUTURE reserves.
	# (futurereserves, currentreserves, currentoutput)
	defaultflowutility=Array{Float64}(resnum,resnum,ynum)
	broadcast!(+, defaultflowutility, -model.grids.reserves./(1+model.econparams.rfree), model.grids.reserves', reshape(model.grids.ydefault,1,1,ynum))
	defaultflowutility[defaultflowutility.<0.0]=0.0
	defaultflowutility=defaultflowutility.^(1-model.econparams.ggamma)./(1-model.econparams.ggamma).*(1-model.econparams.bbeta)
	# Preallocation for default solver
	interimvaluenom=Array{Float64}(debtnum,resnum,ynum,3)
	interimdefaultvalue=Array{Float64}(resnum,ynum,2) 
    reservesmaxtemp=Array{Float64}(1, resnum)
    reservesidtemp=Array{Int64}(1, resnum)
    resrestemp=Array{Float64}(resnum, resnum)

    println("	valuegap	|	pricegap	|	iternum")

    # Main loop

	while resiternum<itermax && ( valuegap>=model.compuparams.valtol || pricegap>=model.compuparams.valtol || !policiesout )
		resiternum+=1
		if mod1(resiternum,intermediatesave)==intermediatesave
			jldopen("debugoldmodel.jld", "w") do file
				write(file, "oldmodel", model)
				write(file, "valuegap", valuegap)
				write(file, "pricegap", 100*pricegap)  
				write(file, "iternum", resiternum-1)
			end
		end
		(mod1(resiternum,iterprint)==1) && tic()
		# 0. Just to make sure in the last iteration policies are recovered
		(valuegap<model.compuparams.valtol) && (pricegap<model.compuparams.valtol) && (policiesout=true)
		(resiternum==itermax-1) && (policiesout=true)
		# 1. Expectation on exogenous varaibles
		valuepayexpectation!(expectedvaluepay, model, interimvaluenom)
		# 2. Solve value of default
		defaultgap=solvedefaultvalue!(model, expectedvaluepay, defaultflowutility, 
								newvaluedefault, 
								interimdefaultvalue, reservesmaxtemp, reservesidtemp, resrestemp)
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
		pmap(updatevaluepaydrm!, [ sub(newvaluepay, :, :, :, iyr) for iyr=1:2*ynum], [sub(newbondprice, :, :, iyr) for iyr=1:2*ynum], # Outputs
				[sub(debtpolicy, :, :, :, iyr) for iyr=1:2*ynum], [sub(reservespolicy, :, :, :, iyr) for iyr=1:2*ynum], # Outputs
				[sub(defaultpolicy, :, :, :, iyr) for iyr=1:2*ynum],  # Outputs
				[expectedvaluepay[ :, :, iyr] for iyr=1:2*ynum], [model.cashinhandpay[ :, :, mod1(iyr,ynum)] for iyr=1:2*ynum],
				[model.bondprice[ :, :, iyr] for iyr=1:2*ynum], [newvaluedefault[ :, iyr] for iyr=1:2*ynum], 
				[model.policies.reservesindefault[:, iyr] for iyr=1:2*ynum], [cld(iyr, ynum) for iyr=1:2*ynum], 
				repeated( model.econparams, 2*ynum), repeated( model.compuparams, 2*ynum), repeated( model.grids, 2*ynum), repeated(policiesout, 2*ynum) )
		# 4. Find new price: take expectation over regime and output
		priceexpectation!(newbondprice, interimvaluenom,
							model.econparams, model.grids.ytrans, model.compuparams.ynum)
		# 5. Find gaps and update 
		axpy!(-1.0, newvaluepay, model.valuepay)
		maxabs!(sub(reservesmaxtemp,1:1), model.valuepay)
		valuegap=reservesmaxtemp[1]/(1-model.econparams.bbeta) # Because for higher beta changes in V are more meaningful
		setindex!(model.valuepay, newvaluepay, :)
		# Cannot do the same, old price cannot overwritten because of updatespeed. Recall interimvaluenom[:,:,:,1:2] also holds newbondprice
		axpy!(-1.0, model.bondprice, sub(interimvaluenom, : ,: ,: , 1:2) )
		maxabs!(sub(reservesmaxtemp,1:1), interimvaluenom[:,:,:,1:2])
		pricegap=reservesmaxtemp[1]
		# Update Control
		scal!(debtnum*resnum*ynum*2, 1-model.compuparams.updatespeed*138/(resiternum+138), model.bondprice, 1 )
		axpy!(model.compuparams.updatespeed*138/(resiternum+138), newbondprice, model.bondprice )
		# update policies
	    if policiesout
			setindex!(model.policies.debt, debtpolicy, :)
			setindex!(model.policies.reserves, reservespolicy, :)
			setindex!(model.policies.default, defaultpolicy, :)
		end		
    	# 5.1 Print intermediate output
    	if mod1(resiternum,iterprint)==iterprint
	    	print("	")
		    showcompact(valuegap)
	    	print("	|	")
	    	showcompact(100*pricegap) 
	    	print("	| 	")
	    	showcompact(resiternum)
	    	print("	| ")
	    	toc()
	    	println()
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

	if resiternum%iterprint!=0 # Print last iteration
    	print("	")
	    showcompact(valuegap)
    	print("	|	")
    	showcompact(100*pricegap) 
    	print("	|	")
    	showcompact(resiternum)
    	print("	|	")
    	toc()
    	println("|")
	end
	policiesout && max(valuegap,pricegap)
end # Function End

