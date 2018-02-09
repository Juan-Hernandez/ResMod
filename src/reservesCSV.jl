# reservesCSV export script
# May want to make this a function

using ReservesTypes
using JLD
using DataFrames

workdir=pwd()
cd(homedir()"\\dropbox\\U-penn\\research\\ReservesProject\\Julia\\Results")
@load "solved100.jld"
cd(workdir)

simulationframe=DataFrame(Debt=0.25*basemodel.grids.debt[ basesimul.debtind[1:end-1] ],
							Reserves=0.25*basemodel.grids.reserves[ basesimul.reservesind[1:end-1] ],
							Mshock=basemodel.grids.mmidpoints[basesimul.mind],
							Output=basemodel.grids.y[basesimul.yind]+basemodel.grids.mmidpoints[basesimul.mind],
							Regime=basesimul.regime,
							DefaultState=basesimul.defaultstate[1:end-1],
							Bondprice=basesimul.bondprice,
							Bondspread=basesimul.bondspread )

writetable("newsimulatedreserves.csv", simulationframe)