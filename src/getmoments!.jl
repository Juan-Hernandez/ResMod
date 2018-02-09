function getmoments!(moments::ModelMoments, oldsimul::ModelSimulation, grids::ModelGrids, burnin::Int64)
	# Start index
	simulated=ModelSimulation(oldsimul.periods-burnin)
	simulated.debtind[:]=oldsimul.debtind[(burnin+1):end]
	simulated.reservesind[:]=oldsimul.reservesind[(burnin+1):end]
	simulated.mind[:]=oldsimul.mind[(burnin+1):end]
	simulated.yind[:]=oldsimul.yind[(burnin+1):end]
    simulated.regime[:]=oldsimul.regime[(burnin+1):end]
	simulated.defaultstate[:]=oldsimul.defaultstate[(burnin+1):end]
	simulated.randomshocks[:]=oldsimul.randomshocks[(burnin+1):oldsimul.periods, 1:4]
	simulated.output[:]=oldsimul.output[(burnin+1):end]
	simulated.bondprice[:]=oldsimul.bondprice[(burnin+1):end]
	simulated.bondspread[:]=oldsimul.bondspread[(burnin+1):end]

	# Mean Indices
	moments.debtmeanind=round(Int64, mean(simulated.debtind[!simulated.defaultstate]) )
	moments.reservesmeanind=round(Int64, mean(simulated.reservesind[!simulated.defaultstate]) )
	moments.mmeanind=round(Int64, mean(simulated.mind) )
	moments.ymeanind=round(Int64, mean(simulated.yind) )
	# Mean Values to potential yearly gdp
	moments.debtmean=0.25*mean(grids.debt[ simulated.debtind[!simulated.defaultstate] ] ) 
	moments.reservesmean=0.25*mean(grids.reserves[ 	simulated.reservesind[!simulated.defaultstate] ] )
	# Output just to check
	moments.outputmean=mean(simulated.output)
	# Spread Mean
	moments.spreadmean=mean(simulated.bondspread[!simulated.defaultstate[1:(end-1)] ] ) 
	
	moments.defaultstatemean=mean(simulated.defaultstate)
	moments.defaultchoicemean=sumabs( diff(simulated.defaultstate) )/simulated.periods/2
	# mean ratios

	moments.debt2gdpmean=0.25*mean(grids.debt[ simulated.debtind[1:(end-1)][!simulated.defaultstate[1:(end-1)]] ]
									./simulated.output[!simulated.defaultstate[1:(end-1)]])
	moments.reserves2gdpmean=0.25*mean(grids.reserves[ simulated.reservesind[1:(end-1)][!simulated.defaultstate[1:(end-1)]] ]
										./simulated.output[!simulated.defaultstate[1:(end-1)]])
	# Variances 
	moments.outputsigma=std(simulated.output[!simulated.defaultstate[1:(end-1)] ])
	moments.spreadsigma=std(simulated.bondspread[!simulated.defaultstate[1:(end-1)] ])
	# Correlations
	moments.spreadXgdp=cor( simulated.bondspread[!simulated.defaultstate[1:(end-1)] ], 
								simulated.output[!simulated.defaultstate[1:(end-1)] ] )
	moments.spreadXgrowth=cor( simulated.bondspread[!simulated.defaultstate[1:(end-1)] ][2:end],
								diff( log(simulated.output[!simulated.defaultstate[1:(end-1)] ]) ) )
	moments.deltaspreadXgdp=cor( diff( simulated.bondspread[!simulated.defaultstate[1:(end-1)] ] ), 
								simulated.output[!simulated.defaultstate[1:(end-1)] ][2:end] )
	moments.deltaspreadXgrowth=cor( diff(simulated.bondspread[!simulated.defaultstate[1:(end-1)] ]),
								diff( log(simulated.output[!simulated.defaultstate[1:(end-1)] ]) ) )
end
