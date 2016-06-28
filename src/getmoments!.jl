function getmoments!(moments::ModelMoments, simulated::ModelSimulation, grids::ModelGrids, burnin::Int64)
	# Start index

	burnin+=1
	# Mean Indices
	moments.debtmeanind=round(Int64, mean(simulated.debtind[!simulated.defaultstate][burnin:end]) )
	moments.reservesmeanind=round(Int64, mean(simulated.reservesind[!simulated.defaultstate][burnin:end]) )
	moments.mmeanind=round(Int64, mean(simulated.mind[burnin:end]))
	moments.ymeanind=round(Int64, mean(simulated.yind[burnin:end]) )
	# Mean Values to potential yearly gdp
	moments.debtmean=0.25*mean(grids.debt[ simulated.debtind[!simulated.defaultstate][burnin:end] ] ) 
	moments.reservesmean=0.25*mean(grids.reserves[ 	simulated.reservesind[!simulated.defaultstate][burnin:end] ] )
	# Output just to check
	moments.outputmean=mean(simulated.output[burnin:end])
	# Spread Mean
	moments.spreadmean=mean(simulated.bondspread[!simulated.defaultstate[1:(end-1)]][burnin:end])
	
	moments.defaultstatemean=mean(simulated.defaultstate[burnin:end])
	moments.defaultchoicemean=sumabs( diff(simulated.defaultstate[burnin:end]) )/(simulated.periods-burnin)/2

	# Variances 
	moments.outputsigma=std(simulated.output[!simulated.defaultstate[1:(end-1)] ][burnin:end])
	moments.spreadsigma=std(simulated.bondspread[!simulated.defaultstate[1:(end-1)] ][burnin:end])
	
end
