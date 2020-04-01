function findnewprice!(newbondprice::Array{Float64,4}, # Output
						bondcashflow::Array{Float64,5}, oldbondprice::Array{Float64,4}, grids::ModelGrids, # Inputs
						grossrate::Float64, ggammalender::Float64, wealthmean::Float64, # Inputs
						expwtilde::Array{Float64,2}) # preallocated temps

# 0. Unpack 
(debtnum,resnum,mnum,ynum,regimenum)=size(bondcashflow)
fill!(newbondprice, 0.0)
fill!(expwtilde, 0.0)
wtilde::Float64=0.0
exotrans::Float64=0.0
# 1. Main Loop
for ioldregime=1:regimenum
	for ioldy=1:ynum
		# Here fixing the previous exogenous state. 
		for iregime=1:regimenum
			for iy=1:ynum
				for imshock=1:mnum
				# Probability of transitioning from previuos old exogenous state to new exogenous state	
				exotrans=grids.mmass[imshock]*grids.ytrans[ioldy, iy]*grids.regimetrans[ioldregime, iregime]
					for ires=1:resnum
						for idebt=1:debtnum
							# 2. Create wtilde^(-gamma) from bondcashflow and others. 
							# @inbounds wtildegamma=( grossrate*wealthmean 
							# 	 					-grids.debt[idebt]*( grossrate*oldbondprice[idebt, ires, ioldy, ioldregime]
							# 	 										-bondcashflow[idebt, ires, imshock, iy, iregime]) 
							# 						)^(-ggammalender)
							#
							# 3. Create expwtilde=E[wtilde^(-gamma)| ioldy, ioldregime]
							# @inbounds expwtilde[idebt, ires]+=exotrans*wtildegamma
							# 4. Create marginal utility of bond: bondcashflow*wtilde^(-gamma), divided by grossrate.  
							# 5. Create E[bondcashflow*wtilde^(-gamma)| ioldy, ioldregime]/(1+r)			
							# @inbounds newbondprice[idebt, ires, ioldy, ioldregime]+=exotrans*wtildegamma/grossrate
							# 															*bondcashflow[idebt, ires, imshock, iy, iregime]
							###################################
							###
							###		USE ggammalender=2
							###
							###################################
							# 2. Create wtilde from bondcashflow and others. 
							@inbounds wtilde=grossrate*wealthmean-grids.debt[idebt]*( grossrate*oldbondprice[idebt, ires, ioldy, ioldregime]
								 														-bondcashflow[idebt, ires, imshock, iy, iregime]) 
							# 3. Create expwtilde=E[wtilde^(-gamma)| ioldy, ioldregime]
							@inbounds expwtilde[idebt, ires]+=exotrans/wtilde/wtilde
							# 4. Create marginal utility of bond: bondcashflow*wtilde^(-gamma), divided by grossrate.  
							# 5. Create E[bondcashflow*wtilde^(-gamma)| ioldy, ioldregime]/(1+r)			
							@inbounds newbondprice[idebt, ires, ioldy, ioldregime]+=exotrans/wtilde/wtilde/grossrate*bondcashflow[idebt, ires, imshock, iy, iregime]
						end
					end
				end
			end
		end
		# 6. New price: E[bondcashflow*wtilde^(-gamma)| ioldy, ioldregime] / ((1+r)expwtilde)
		for ifutres=1:resnum
			for ifutdebt=1:debtnum
				@inbounds newbondprice[ifutdebt, ifutres, ioldy, ioldregime]= newbondprice[ifutdebt, ifutres, ioldy, ioldregime]/expwtilde[ifutdebt, ifutres]
				# Clear expwtilde for next loop -- better than fill! because loop is already here
				@inbounds expwtilde[ifutdebt, ifutres]=0.0
			end
		end
	end
end
nothing
end # Function end