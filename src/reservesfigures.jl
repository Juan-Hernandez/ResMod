function reservesfigures(model::ReservesModel, simulated::ModelSimulation, moments::ModelMoments)
	# 0. Preliminaries
	mmass=model.grids.mmass


	# 1. Extra matrices for plot	
	bondspread=Array{Float64}( model.compuparams.debtnum, model.compuparams.resnum, model.compuparams.ynum, 2)
	broadcast!(/, bondspread, model.econparams.coupon*(1-model.econparams.llambda)+model.econparams.llambda , model.bondprice)
	broadcast!(+, bondspread, bondspread, 1 , -model.econparams.llambda)
	broadcast!(^, bondspread, bondspread, 4)
	broadcast!(+, bondspread, bondspread, -(1+model.econparams.rfree)^4)
	broadcast!(*, bondspread, bondspread, 10000)
	
	# 2. Point Values for figures

	# In general output mean is close to YN+1/2
	outputmeanind=findfirst(x->(x>moments.outputmean),model.grids.y)
	
	outputsigmalowind=findfirst(x->(x>moments.outputmean-1.6*moments.outputsigma),model.grids.y)
	ylow16stdind=findlast( x->(x<moments.outputmean-1.6*moments.outputsigma)  ,model.grids.y)
	yfirstpayind=findfirst(!model.policies.default[moments.debtmeanind, moments.reservesmeanind, moments.mmeanind, :, 1] )
	ylowstdind=max(ylow16stdind, yfirstpayind)

	# Next debt and reserves levels at  debt, reserves; mean m,
	# mean and 1.6sd lower income, both regimes

	debtnextmean=model.policies.debt[moments.debtmeanind, moments.reservesmeanind, moments.mmeanind, moments.ymeanind, 1]
	debtnextylowstd=model.policies.debt[moments.debtmeanind, moments.reservesmeanind, moments.mmeanind,ylowstdind, 1]
	debtnextmeanss=model.policies.debt[moments.debtmeanind, moments.reservesmeanind, moments.mmeanind, moments.ymeanind, 2]
	debtnextylowstdss=model.policies.debt[moments.debtmeanind, moments.reservesmeanind, moments.mmeanind, ylowstdind, 2]

	reservesnextmean=model.policies.reserves[moments.debtmeanind, moments.reservesmeanind, moments.mmeanind, moments.ymeanind, 1]
	reservesnextylowstd=model.policies.reserves[moments.debtmeanind, moments.reservesmeanind, moments.mmeanind, ylowstdind, 1]
	reservesnextmeanss=model.policies.reserves[moments.debtmeanind, moments.reservesmeanind, moments.mmeanind, moments.ymeanind, 2]
	reservesnextylowstdss=model.policies.reserves[moments.debtmeanind, moments.reservesmeanind, moments.mmeanind, ylowstdind, 2]

	# 3. Average mshock out on policies
	debtpolicysmooth=Array{Float64}(size(model.bondprice))
	reservespolicysmooth=Array{Float64}(size(model.bondprice))
	defaultpolicysmooth=Array{Float64}(size(model.bondprice))
	mexpectation!(debtpolicysmooth, model.grids.debt[model.policies.debt], model.grids.mmass)
	mexpectation!(reservespolicysmooth, model.grids.reserves[model.policies.reserves], model.grids.mmass)
	mexpectation!(defaultpolicysmooth, model.policies.default, model.grids.mmass)

	# THIS IS probability of default next period for all possible endogenous state choices. 
	# Expectation is taken over exogenous states
	defaultnextprob=Array{Float64}(model.compuparams.debtnum, model.compuparams.resnum, model.compuparams.ynum, 2)
	tempdry=Array{Float64}(model.compuparams.debtnum, model.compuparams.resnum, model.compuparams.ynum)
	ywexpectation!(defaultnextprob, defaultpolicysmooth, 
							model.grids.ytrans, model.grids.regimetrans, 1.0,
							tempdry)
	
	# 4. Figures
	""" Make data frames: better for plotting with Gadfly """
	# Spreads too big overflow dataframe or plot
	bondspread[bondspread.>1.0e10]=Inf
	# Recall price and spread are contingent on current income and regime but future debt and reserves
	
	# 4.1 Spread vs next period debt
	# Make dataframe for figure
	figdata=DataFrame(Debt=repmat(0.25*model.grids.debt, 4), 
						Spread=[bondspread[:, reservesnextmean, moments.ymeanind, 1]; bondspread[:, reservesnextmeanss, moments.ymeanind, 2];
								bondspread[:, reservesnextylowstd, ylowstdind, 1]; bondspread[:, reservesnextylowstdss, ylowstdind, 2] ],
						ExoStates=repeat(["Yavg";"YavgPanic";"Ylow";"YlowPanic"], inner=[model.compuparams.debtnum]) )
	# Just bind x and y aesthetic for points
	figpointsx=0.25*model.grids.debt[ [debtnextmean debtnextmeanss debtnextylowstd  debtnextylowstdss] ]
	figpointsy=[bondspread[debtnextmean, reservesnextmean, moments.ymeanind, 1];
				bondspread[debtnextmeanss, reservesnextmeanss, moments.ymeanind, 2];
				bondspread[debtnextylowstd, reservesnextylowstdss, ylowstdind, 1]
				bondspread[debtnextylowstdss, reservesnextylowstdss, ylowstdind, 2] ]
	# Make layers
	lineslayer=layer(figdata, x=:Debt, y=:Spread, color=:ExoStates, Geom.line)
	pointslayer=layer(x=figpointsx, y=figpointsy, Geom.point)
	# Make plot
	figure1=plot(lineslayer, pointslayer, Coord.cartesian(xmax=0.25*model.compuparams.debtmax+0.01, ymax=1000),
					Guide.title("Spread given next period Debt"))
	# Save plot
	currdir=pwd()
	cd("..\\..\\figures")
	draw(SVG("Spread-debt.svg", 9inch, 6inch), figure1)	

	# 3.2 Spread vs Next Period Reserves 
	# Make dataframe for figure
	figdata=DataFrame(Reserves=repmat(0.25*model.grids.reserves, 4), 
						Spread=vec([bondspread[debtnextmean, :, moments.ymeanind, 1]' bondspread[debtnextmeanss, :, moments.ymeanind, 2]'
								bondspread[debtnextylowstd, :, ylowstdind, 1]' bondspread[debtnextylowstdss, :, ylowstdind, 2]' ]),
						DefProb=vec([defaultnextprob[debtnextmean, :, moments.ymeanind, 1]' defaultnextprob[debtnextmeanss, :, moments.ymeanind, 2]'
								defaultnextprob[debtnextylowstd, :, ylowstdind, 1]' defaultnextprob[debtnextylowstdss, :, ylowstdind, 2]' ]),
						NextDebt=0.25*vcat(slice(debtpolicysmooth, moments.debtmeanind, :, moments.ymeanind, 1), 
						 					slice(debtpolicysmooth, moments.debtmeanind, :, moments.ymeanind, 2),
		 									slice(debtpolicysmooth, moments.debtmeanind, :, ylowstdind, 1), 
						 					slice(debtpolicysmooth, moments.debtmeanind, :, ylowstdind, 2) ),
						ExoStates=repeat(["Yavg";"YavgPanic";"Ylow";"YlowPanic"], inner=[model.compuparams.resnum]) )
	# Just bind x and y aesthetic for points
	figpointsx=0.25*model.grids.reserves[ [reservesnextmean reservesnextmeanss reservesnextylowstd reservesnextylowstdss] ]
	# Same figpointsy as in figure 1
	# Make layers
	lineslayer=layer(figdata, x=:Reserves, y=:Spread, color=:ExoStates, Geom.line)
	pointslayer=layer(x=figpointsx, y=figpointsy, Geom.point)
	# Make plot
	figure2=plot(lineslayer, pointslayer, Coord.cartesian(xmax=0.25*model.compuparams.resmax+0.01, ymax=1000),
					Guide.title("Spread given next period Reserves"))
	# Save plot
	draw(SVG("Spread-Reserves.svg", 9inch, 6inch), figure2)		

	# 3.3 Defalult prob next period given reserves
	# Make dataframe for figure

	# Just bind x and y aesthetic for points
	# same figpointsx as before 
	figpointsy=[defaultnextprob[debtnextmean, reservesnextmean, moments.ymeanind, 1];
				defaultnextprob[debtnextmeanss, reservesnextmeanss, moments.ymeanind, 2];
				defaultnextprob[debtnextylowstd, reservesnextylowstd, ylowstdind, 1]
				defaultnextprob[debtnextylowstdss, reservesnextylowstdss, ylowstdind, 2] ]
	# Make layers
	lineslayer=layer(figdata, x=:Reserves, y=:DefProb, color=:ExoStates, Geom.line)
	pointslayer=layer(x=figpointsx, y=figpointsy, Geom.point)
	# Make plot
	figure3=plot(lineslayer, pointslayer, Coord.cartesian(xmax=0.25*model.compuparams.resmax+0.01, ymax=1),
					Guide.title("Default probability given next period Reserves"))
	# Save plot
	draw(SVG("DefProb-Reserves.svg", 9inch, 6inch), figure3)	

	# 3.4 Debt given Reserves THIS IS CRAZY
	# Plot next period debt given current reserves, zero temporary shock around mean debt mean level
	# Same dataframe
	# Just bind x and y aesthetic for points
	# same figpointsx as before 
	figpointsy=0.25*model.grids.debt[ [debtnextmean, debtnextmeanss, debtnextylowstd, debtnextylowstdss] ]
	# Make layers
	lineslayer=layer(figdata, x=:Reserves, y=:NextDebt, color=:ExoStates, Geom.line)
	pointslayer=layer(x=figpointsx, y=figpointsy, Geom.point)
	# Make plot
	figure4=plot(lineslayer, pointslayer, Guide.title("Next Debt given current Reserves"),
					Coord.cartesian(xmax=0.25*model.compuparams.resmax+0.01, ymax=0.25*model.compuparams.debtmax+0.01) )
	# Save plot
	draw(SVG("NextDebt-Reserves.svg", 9inch, 6inch), figure4)	

	# 3.5 Debt given current output
	defoutput1=findfirst(!model.policies.default[moments.debtmeanind, moments.reservesmeanind, moments.mmeanind, :, 1])
	defoutput2=findfirst(!model.policies.default[moments.debtmeanind, moments.reservesmeanind, moments.mmeanind, :, 2])
	


	# Make dataframe for figure
	figdata=DataFrame(Output=repmat(model.grids.y, 2), 
						NextDebt=0.25*vec(slice(debtpolicysmooth,moments.debtmeanind, moments.reservesmeanind, :, 1:2) ),
						NextReserves=0.25*vec(slice(reservespolicysmooth, moments.debtmeanind, moments.reservesmeanind, :, 1:2 ) ),
						Sunspot=repeat(["Normal";"Panic"], inner=[model.compuparams.ynum]),
						Default=vec(slice(defaultpolicysmooth, moments.debtmeanind, moments.reservesmeanind, :, 1:2) ) )
	figdata2=by(figdata,[:Sunspot, :Default], df->string(df[1,:Sunspot],'-',mean(df[1,:Default]) ) )
	figdata=join(figdata, figdata2, on=[:Sunspot, :Default], kind=:inner)
	# Just bind x and y aesthetic for points
	figpointsx=model.grids.y[ [moments.ymeanind moments.ymeanind ylowstdind ylowstdind] ]
	# Same figpointsy as in figure 4
	# Make layers
	lineslayer=layer(figdata, x=:Output, y=:NextDebt, color=:x1, Geom.line)
	pointslayer=layer(x=figpointsx, y=figpointsy, Geom.point)
	# Make plot
	figure5=plot(lineslayer, pointslayer, Coord.cartesian(xmin=model.grids.y[1], xmax=model.grids.y[end]), 
					Guide.title("Next Debt given current Output"))
	# Save plot
	draw(SVG("NextDebt-Output.svg", 9inch, 6inch), figure5)		

	# 3.6 Debt given current output
	# Same dataframe
	# Just bind x and y aesthetic for points
	# same figpointsx as before 
	figpointsy=0.25*model.grids.reserves[ [reservesnextmean reservesnextmeanss reservesnextylowstd reservesnextylowstdss] ]
	# Make layers
	lineslayer=layer(figdata, x=:Output, y=:NextReserves, color=:x1, Geom.line)
	pointslayer=layer(x=figpointsx, y=figpointsy, Geom.point)
	figure6=plot(lineslayer, pointslayer, Coord.cartesian(xmin=model.grids.y[1], xmax=model.grids.y[end]),
					Guide.title("Next Reserves given current Output"))
	draw(SVG("NextReserves-Output.svg", 9inch, 6inch), figure6)		
# %% Reserves-Output
# D1=find(Policy.Default(MmeanId,:,debtMeanId,reservesMeanId,1)<1,1);
# D2=find(Policy.Default(MmeanId,:,debtMeanId,reservesMeanId,2)<1,1);
# figOutreserves=figure('Position',[50,50,450,350],'Name',...
#     'reserves-Output','Color','w');
# plot(Ygrid(1:D1),[model.grids.reserves(Policy.Reserves(MmeanId,[1:(D1-1),(D1-1)],debtMeanId,reservesMeanId,1))],'LineWidth',2);
#     hold all
#     %plot(Ygrid(D1:YN),model.grids.reserves(Policy.Reserves(MmeanId,D1:YN,debtMeanId,reservesMeanId,1)),'LineWidth',2);
#     % M smoothed figure
#     plot(Ygrid(D1:YN),reservesMax/(reservesN-1)*(Mmass'*Policy.Reserves(:,D1:YN,debtMeanId,reservesMeanId,1)-1),'LineWidth',2);
#     line([Ygrid(D1),Ygrid(D1)],[0,0.7],'Color','r','LineStyle','--')
#     % plot(Ygrid(1:D2),[model.grids.reserves(Policy.Reserves(MmeanId,[1:(D2-1),(D2-1)],debtMeanId,reservesMeanId,2))],'LineWidth',2);
#     % plot(Ygrid(D2:YN),model.grids.reserves(Policy.Reserves(MmeanId,D2:YN,debtMeanId,reservesMeanId,2)),'LineWidth',2);
#     line([Ygrid(D2),Ygrid(D2)],[0,0.7],'Color','r','LineStyle','--')
#     legend('Default','No default','Threshold','Location','East')
#     title('Next Period Reserves Given Current Output')
#     xlabel('Output');
#     ylabel('Next Period Reserves');
#     axis([Ygrid(1),Ygrid(YN),0,0.7]);
# %    plot(model.grids.reserves(2),bondspread(11,D1,2),'bo','MarkerFaceColor','b');
# %    plot(model.grids.reserves(1),bondspread(8,D2,1),'go','MarkerFaceColor','g');
# %cd 'C:\Users\user\Dropbox\U-Penn\JMP\ReservesProject\Latex';
#     set(figOutreserves,'PaperPositionMode','auto')
#     saveas(figOutreserves,'Reserves-Output.png')
# %cd(CurrDir);

















	cd(currdir)
end







# %%



# %% JUST TO PLAY 
# % debt - Output
# D1=find(Policy.Default(MmeanId,:,debtMeanId,4,1)<1,1);
# D2=find(Policy.Default(MmeanId,:,debtMeanId,4,2)<1,1);
# figOutdebt=figure('Position',[50,50,450,350],'Name',...
#     'debt-Output','Color','w');
# plot(Ygrid(1:D1),[model.grids.debt(Policy.debt(MmeanId,1:(D1-1),debtMeanId,4,1));0],'LineWidth',2);
#     hold all
# plot(Ygrid(D1:YN),model.grids.debt(Policy.debt(MmeanId,D1:YN,debtMeanId,4,1)),'LineWidth',2);
#     line([Ygrid(D1),Ygrid(D1)],[0,model.grids.debt(debtMaxId)+0.02],'Color','r','LineStyle','--')
# plot(Ygrid(1:D2),[model.grids.debt(Policy.debt(MmeanId,1:(D2-1),debtMeanId,4,2));0],'LineWidth',2);
# plot(Ygrid(D2:YN),model.grids.debt(Policy.debt(MmeanId,D2:YN,debtMeanId,4,2)),'LineWidth',2);
#     line([Ygrid(D2),Ygrid(D2)],[0,model.grids.debt(debtMaxId)+0.02],'Color','y','LineStyle','--')

#     legend('Default','No default','Threshold','Location','East')
#     title('Next Period debt Given Current Output')
#     xlabel('Output');
#     ylabel('Next Period debt');
#     axis([0.775,1.225,0,model.grids.debt(debtMaxId)+0.02]);
# %    plot(model.grids.reserves(2),bondspread(moments.ymeanind,D1,2),'bo','MarkerFaceColor','b');
# %    plot(model.grids.reserves(1),bondspread(ylowstdind,D2,1),'go','MarkerFaceColor','g');




# end