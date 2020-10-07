function getmoments!(moments::ModelMoments, oldsimul::ModelSimulation, grids::ModelGrids, burnin::Int64)
	# Local


	# 1. Create new simulation structure eliminating burnin
	simulated=ModelSimulation(oldsimul.periods-burnin)
	simulated.debtind[:]=oldsimul.debtind[(burnin+1):end]
	simulated.reservesind[:]=oldsimul.reservesind[(burnin+1):end]
	simulated.mind[:]=oldsimul.mind[(burnin+1):end]
	simulated.yind[:]=oldsimul.yind[(burnin+1):end]
    simulated.regime[:]=oldsimul.regime[(burnin+1):end]
	simulated.randomshocks[:]=oldsimul.randomshocks[(burnin+1):end, :]
	simulated.output[:]=oldsimul.output[(burnin+1):end]
	simulated.bondprice[:]=oldsimul.bondprice[(burnin+1):end]
	simulated.bondspread[:]=oldsimul.bondspread[(burnin+1):end]
	simulated.defaultstate[:]=oldsimul.defaultstate[(burnin+1):end]
	deltaspread=diff(oldsimul.bondspread[burnin:end])
	growth=diff(log.(oldsimul.output[burnin:end]))
	# 2. Crete repayment moment indicator
	nodefindicator=falses(simulated.periods)
	for index=1:simulated.periods
		# Check for any defaultstate in the last 24 quarters (6 years). If true, then exclude !(). 
		nodefindicator[index]=!( maximum(oldsimul.defaultstate[(burnin-23+index):(burnin+index)]) )
	end

	# 3. Mean Indices
	moments.debtmeanind=round(Int64, mean(simulated.debtind[nodefindicator]) )
	moments.reservesmeanind=round(Int64, mean(simulated.reservesind[nodefindicator]) )
	moments.mmeanind=round(Int64, mean(simulated.mind) )
	moments.ymeanind=round(Int64, mean(simulated.yind) )

	# 4. Mean values to average yearly gdp (4.0)
	moments.debtmean=0.25*mean(grids.debt[ simulated.debtind[nodefindicator] ] ) 
	moments.reservesmean=0.25*mean(grids.reserves[ simulated.reservesind[nodefindicator] ] )
	# Output just to check
	moments.outputmean=mean(simulated.output)
	# Spread Mean
	moments.spreadmean=mean(simulated.bondspread[nodefindicator] ) 
	# Defauls state and choice mean
	moments.defaultstatemean=mean(simulated.defaultstate)
	moments.defaultchoicemean=sum(abs, diff(simulated.defaultstate) )/simulated.periods/2
	
	# 5. Moments of ratios to GDP
	# 5.1 Means
	moments.debt2gdpmean=0.25*mean(grids.debt[ simulated.debtind[nodefindicator] ]
									./simulated.output[nodefindicator] )
	moments.reserves2gdpmean=0.25*mean(grids.reserves[ simulated.reservesind[nodefindicator] ]
										./simulated.output[nodefindicator])
	# 5.2 Variances 
	moments.outputsigma=std(simulated.output[nodefindicator])
	moments.spreadsigma=std(simulated.bondspread[nodefindicator])
	# 5.3 Correlations
	moments.spreadXgdp=cor( simulated.bondspread[nodefindicator], 
								simulated.output[nodefindicator] )
	moments.spreadXgrowth=cor( simulated.bondspread[nodefindicator],
								growth[nodefindicator] )  
	moments.deltaspreadXgdp=cor( deltaspread[nodefindicator] , 
								simulated.output[nodefindicator] )
	moments.deltaspreadXgrowth=cor( deltaspread[nodefindicator],
								growth[nodefindicator] )
	
end
