function solvereservesmodel!(model::ReservesModel, solverparams=SolverParams())

	resiternum::Int64=solverparams.startiternum
	iterprint::Int64=solverparams.iterprint 
	itermax::Int64=solverparams.itermax 
	intermediatesave::Int64=solverparams.intermediatesave
	policiesout::Bool=solverparams.policiesout
	updatespeed::Float64=solverparams.updatespeed
	solveroutvec::Array{Float64}(4)
	printbool::Bool=iterprint!=0

	# Unpack counters
	debtnum::Int64=model.compuparams.debtnum
	resnum::Int64=model.compuparams.resnum
	ynum::Int64=model.compuparams.ynum
	mnum::Int64=model.compuparams.mnum
	regimenum::Int64=size(model.grids.regimetrans,1)
	exonum::Int64=regimenum*ynum

	# Holder for new value functions
	newbondprice=Array{Float64}(debtnum, resnum, ynum, regimenum)
	bondcashflow=SharedArray{Float64}(debtnum, resnum, mnum, ynum, regimenum)
	newvaluedefault=Array{Float64}(resnum, ynum, regimenum)
	newvaluepay=SharedArray{Float64}(debtnum, resnum, mnum, ynum, regimenum)
	# Also need to allocate policies in shared arrays
	debtpolicy=SharedArray{Int64}(debtnum, resnum, mnum, ynum, regimenum)
	reservespolicy=SharedArray{Int64}(debtnum, resnum, mnum, ynum, regimenum)
	defaultpolicy=SharedArray{Bool}(debtnum, resnum, mnum, ynum, regimenum)

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
	broadcast!(+, defaultflowutility, -model.grids.reserves./(1.0+model.econparams.rfree), model.grids.reserves', reshape(model.grids.ydefault,1,1,ynum))
	defaultflowutility[defaultflowutility.<0.0]=0.0
	defaultflowutility=defaultflowutility.^(1-model.econparams.ggamma)./(1.0-model.econparams.ggamma).*(1.0-model.econparams.bbeta)
	# Preallocation of temporary arrays
	tempdry=Array{Float64}(debtnum,resnum,ynum)
	tempdryw=Array{Float64}(debtnum,resnum,ynum,regimenum)
	temp2dryw=Array{Float64}(debtnum,resnum,ynum,regimenum)
	tempdrmyw=Array{Float64}(debtnum,resnum,mnum,ynum,regimenum)
	
	tempry=Array{Float64}(resnum,ynum) 
	temprr=Array{Float64}(resnum, resnum)
	tempryw=Array{Float64}(resnum,ynum,regimenum) 
	
    reservesmaxtemp=Array{Float64}(1, resnum)
    reservesidtemp=Array{Int64}(1, resnum)
    

    printbool && println("        valuegap         |        pricegap         |   iternum  | time (sec) |")

    # Main loop
    tic()
	while resiternum<itermax && (max(valuegap, pricegap, defaultgap)>=solverparams.valtol || !policiesout)
		# 0. Output and update controls
		# Increase iteration counter
		resiternum+=1
		# Set intermediate save
		if mod1(resiternum,intermediatesave)==intermediatesave
			jldopen("debugoldmodel.jld", "w") do file
				write(file, "oldmodel", model)
				write(file, "valuegap", valuegap)
				write(file, "pricegap", 100.0*pricegap)  
				write(file, "iternum", resiternum-1)
			end
		end
		# Intermediate print and timing
		printbool && ( (mod1(resiternum,iterprint)==1) && tic() )
		# Policies out in last iteration
		(max(valuegap, pricegap, defaultgap)<=solverparams.valtol) && (policiesout=true)
		(resiternum==itermax) && (policiesout=true)

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
			# Careful, use view() for those arrays you want to be referenced. [:,:,:] is good for those that can be copied.
			# Since we are pmapping, 		
""" For prototyping we use just one evaluation, pmap below"""		
		# updatevaluepaydrm!( view(newvaluepay, :, :, :, 3), view(bondcashflow, :, :, :, 3), 
		# 		view(debtpolicy, :, :, :, 3), view(reservespolicy, :, :, :, 3), view(defaultpolicy, :, :, :, 3), # Outputs
		# 		expectedvaluepay[ :, :, 3, 1], model.cashinhandpay[ :, :, 3], model.bondprice[:,:,3,1],
		# 		newvaluedefault[ :, 3], model.policies.reservesindefault[:, 3], 1, model.econparams, model.compuparams, model.grids, true )		
		# pmap requires shared arrays for inplace! outputs
		pmap(updatevaluepaydrm!, [ view(newvaluepay, :, :, :, iy, ir) for iy=1:ynum, ir=1:regimenum], # Output: new value of paying
			[view(bondcashflow, :, :, :, iy, ir) for iy=1:ynum, ir=1:regimenum], # Output: bond cashflow after choices
			[view(debtpolicy, :, :, :, iy, ir) for iy=1:ynum, ir=1:regimenum], # Output: new debt choice
			[view(reservespolicy, :, :, :, iy, ir) for iy=1:ynum, ir=1:regimenum], # Output: reseves choice
			[view(defaultpolicy, :, :, :, iy, ir) for iy=1:ynum, ir=1:regimenum],  # Output: new default choice 
			[expectedvaluepay[ :, :, iy, ir] for iy=1:ynum, ir=1:regimenum], # Input: current continuation value of paying 
			[model.cashinhandpay[ :, :, mod1(iyr,ynum)] for iyr=1:exonum], # Input: cash in hand after servicing debt
			[model.bondprice[ :, :, iy, ir] for iy=1:ynum, ir=1:regimenum], # Input: current bond price
			[newvaluedefault[ :, iy, ir] for iy=1:ynum, ir=1:regimenum], # Input: value of default from 2.
			[model.policies.reservesindefault[:, iy, ir] for iy=1:ynum, ir=1:regimenum], # Input: reserves choice in default from 2.
			[cld(iyr, ynum) for iyr=1:exonum], # Input: regimenum
			Iterators.repeated( solverparams.valtol, exonum), Iterators.repeated( model.econparams, exonum), Iterators.repeated( model.compuparams, exonum), 
			Iterators.repeated( model.grids, exonum), Iterators.repeated(policiesout, exonum) )
		# 4. Find new price: take expectation over regime and output
		if model.econparams.ggammalender!=0
			# Create wtilde^(-gamma) from bondcashflow ad others
			tempdrmyw.= ( (1.0+model.econparams.rfree).*model.econparams.wealthmean .- model.grids.debt.*
							( (1.0+model.econparams.rfree).*reshape(model.bondprice, debtnum, resnum, 1, ynum, regimenum).-bondcashflow ) 
							).^(-model.econparams.ggammalender)
			# Create marginal utility of bond: cashflow*wtilde^(-gamma). 
			bondcashflow.= tempdrmyw .* bondcashflow 
			# Create E[wtilde^(-gamma)]
			mexpectation!(tempdryw, tempdrmyw, model.grids.mmass)				
			ywexpectation!(temp2dryw, tempdryw, 
							model.grids.ytrans, model.grids.regimetrans, 1.0, 
							tempdry)			
		end
		mexpectation!(tempdryw, sdata(bondcashflow), model.grids.mmass)
		ywexpectation!(newbondprice, tempdryw, 
							model.grids.ytrans, model.grids.regimetrans, 1.0/(1.0+model.econparams.rfree),
							tempdry)
		if model.econparams.ggammalender!=0
			newbondprice.= newbondprice./temp2dryw
		end
		# 5. Find gaps and update 
		Base.LinAlg.BLAS.axpy!(-1.0, newvaluepay, model.valuepay)
		maximum!(view(reservesmaxtemp,1:1), abs.(model.valuepay))
		valuegap=reservesmaxtemp[1]/(1-model.econparams.bbeta) # Because for higher beta changes in V are more meaningful
		setindex!(model.valuepay, newvaluepay, :)
		# Cannot do the same, old price cannot be overwritten because of updatespeed. Recall tempdryw holds new bond price
		setindex!(tempdryw, newbondprice, :)
		Base.LinAlg.BLAS.axpy!(-1.0, model.bondprice, tempdryw )
		maximum!(view(reservesmaxtemp,1:1), abs.(tempdryw))
		pricegap=reservesmaxtemp[1]
		# Update Control
		scal!(debtnum*resnum*ynum*regimenum, 1.0-updatespeed, model.bondprice, 1 )
		Base.LinAlg.BLAS.axpy!(updatespeed, newbondprice, model.bondprice )
		# update policies
	    if policiesout
			setindex!(model.policies.debt, debtpolicy, :)
			setindex!(model.policies.reserves, reservespolicy, :)
			setindex!(model.policies.default, BitArray(defaultpolicy), :)
		end		
    	# 5.1 Print intermediate output
    	if printbool && mod1(resiternum,iterprint)==iterprint
	    	print("       ")
		    showcompact(valuegap)
	    	print("        |       ")
	    	showcompact(100.0*pricegap) 
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
				write(file, "pricegap", 100.0*pricegap)
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
    	showcompact(100.0*pricegap) 
    	print("       |     ")
    	showcompact(resiternum)
    	print("    |   ")
    	showcompact(toq())
    	println("  |")
	end
	solveroutvec = [convert(Float64, resiternum),valuegap, pricegap, defaultgap]
	return solveroutvec
end # Function End

#solvereservesmodel!(model::ReservesModel)=solvereservesmodel!(model, SolverParams( ))
#solvereservesmodel!(model::ReservesModel, x::Any)=solvereservesmodel!(model, SolverParams(x...))


