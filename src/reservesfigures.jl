function reservesfigures(model::ReservesModel, simulated::ModelSimulation, moments::ModelMoments)
	# 0. Preliminaries
	mmass=model.grids.mmass;
	debtpoints=Array{Float64}(undef, 4);
	reservespoints=Array{Float64}(undef, 4);
	outputpoints=Array{Float64}(undef, 4);
	debtnum=model.compuparams.debtnum
	resnum=model.compuparams.resnum
	ynum=model.compuparams.ynum

	# 1. Extra matrices for plot
	bondspread=Array{Float64}(undef,  debtnum, resnum, ynum, 2);
	broadcast!(/, bondspread, model.econparams.coupon*(1-model.econparams.llambda)+model.econparams.llambda , model.bondprice);
	broadcast!(+, bondspread, bondspread, 1 , -model.econparams.llambda);
	broadcast!(^, bondspread, bondspread, 4);
	broadcast!(+, bondspread, bondspread, -(1+model.econparams.rfree)^4);
	broadcast!(*, bondspread, bondspread, 10000);
	
	# 2. Average mshock out on policies
	debtpolicysmooth=Array{Float64}(undef, size(model.bondprice));
	reservespolicysmooth=Array{Float64}(undef, size(model.bondprice));
	defaultpolicysmooth=Array{Float64}(undef, size(model.bondprice));
	mexpectation!(debtpolicysmooth, model.grids.debt[model.policies.debt], model.grids.mmass)
	mexpectation!(reservespolicysmooth, model.grids.reserves[model.policies.reserves], model.grids.mmass)
	mexpectation!(defaultpolicysmooth, model.policies.default, model.grids.mmass)

	# THIS IS probability of default next period for all possible endogenous state choices. ;
	# Expectation is taken over exogenous states;
	defaultnextprob=Array{Float64}(undef, debtnum, resnum, ynum, 2);
	tempdry=Array{Float64}(undef, debtnum, resnum, ynum);
	ywexpectation!(defaultnextprob, defaultpolicysmooth, 
							model.grids.ytrans, model.grids.regimetrans, 1.0,
							tempdry)

	# 3. Point Values for figures

	# In general output mean is close to YN+1/2
	outputmeanind::Int64 = findfirst(x->(x>moments.outputmean), model.grids.y)
	outputsigmalowind::Int64 = findfirst(x->(x>moments.outputmean-1.6*moments.outputsigma), model.grids.y)
	ylow16stdind::Int64 = findlast( x->(x<moments.outputmean-1.6*moments.outputsigma) ,model.grids.y)
	yfirstpayind::Int64 = findfirst(.!model.policies.default[moments.debtmeanind, moments.reservesmeanind, moments.mmeanind, :, 1] )
	ylowstdind::Int64 = max(ylow16stdind, yfirstpayind)

	# Next debt and reserves levels at  debt, reserves; mean m,;
	# mean and 1.6sd lower income, both regimes;

	debtnextmean = debtpolicysmooth[ moments.debtmeanind, moments.reservesmeanind, moments.ymeanind, 1]
	debtnextylowstd = debtpolicysmooth[ moments.debtmeanind, moments.reservesmeanind, ylowstdind, 1]
	debtnextmeanss = debtpolicysmooth[ moments.debtmeanind, moments.reservesmeanind, moments.ymeanind, 2]
	debtnextylowstdss =debtpolicysmooth[ moments.debtmeanind, moments.reservesmeanind, ylowstdind, 2]

	iddebtnextmean = ceil(Int64, dot( view(model.policies.debt, moments.debtmeanind, moments.reservesmeanind, :, moments.ymeanind, 1), model.grids.mmass ) )
	iddebtnextylowstd = ceil(Int64, dot( view(model.policies.debt, moments.debtmeanind, moments.reservesmeanind, :, ylowstdind, 1), model.grids.mmass ) )
	iddebtnextmeanss = ceil(Int64, dot( view(model.policies.debt, moments.debtmeanind, moments.reservesmeanind, :, moments.ymeanind, 2), model.grids.mmass ) )
	iddebtnextylowstdss = ceil(Int64, dot( view(model.policies.debt, moments.debtmeanind, moments.reservesmeanind, :, ylowstdind, 2), model.grids.mmass ) )

	reservesnextmean = reservespolicysmooth[ moments.debtmeanind, moments.reservesmeanind, moments.ymeanind, 1]
	reservesnextylowstd = reservespolicysmooth[ moments.debtmeanind, moments.reservesmeanind, ylowstdind, 1]
	reservesnextmeanss = reservespolicysmooth[ moments.debtmeanind, moments.reservesmeanind, moments.ymeanind, 2]
	reservesnextylowstdss =reservespolicysmooth[ moments.debtmeanind, moments.reservesmeanind, ylowstdind, 2]

	idreservesnextmean::Int64 = ceil(Int64, dot( view(model.policies.reserves, moments.debtmeanind, moments.reservesmeanind, :, moments.ymeanind, 1), model.grids.mmass ) )
	idreservesnextylowstd::Int64 = ceil(Int64, dot( view(model.policies.reserves, moments.debtmeanind, moments.reservesmeanind, :, ylowstdind, 1), model.grids.mmass ) )
	idreservesnextmeanss::Int64 = ceil(Int64, dot( view(model.policies.reserves, moments.debtmeanind, moments.reservesmeanind, :, moments.ymeanind, 2), model.grids.mmass ) )
	idreservesnextylowstdss::Int64 = ceil(Int64, dot( view(model.policies.reserves, moments.debtmeanind, moments.reservesmeanind, :, ylowstdind, 2), model.grids.mmass ) )


	# 4. Figures
	""" Make data frames: better for plotting with Gadfly """
	# Spreads too big overflow dataframe or plot
	bondspread[bondspread.>1.0e10].=Inf;
	# Recall price and spread are contingent on current income and regime but future debt and reserves
	
	latexstyle=style(major_label_font="CMU Serif",minor_label_font="CMU Serif", major_label_font_size=12pt, minor_label_font_size=10pt)

	# 4.1 Spread vs next period debt
	# Make dataframe for figure
	debtframe=DataFrame(Debt=repeat(0.25*model.grids.debt, 4), 
						Spread=[ bondspread[:, idreservesnextmean, moments.ymeanind, 1]; bondspread[:, idreservesnextmeanss, moments.ymeanind, 2];
									bondspread[:, idreservesnextylowstd, ylowstdind, 1]; bondspread[:, idreservesnextylowstdss, ylowstdind, 2] ],
						ExoStates=repeat(["Yavg";"Yavg-Panic";"Ylow";"Ylow-Panic"], inner=[debtnum]) )
	# Just bind x and y aesthetic for points
	debtpoints=0.25*[debtnextmean debtnextmeanss debtnextylowstd  debtnextylowstdss]
	spreadpoints=[bondspread[iddebtnextmean, idreservesnextmean, moments.ymeanind, 1];
				bondspread[iddebtnextmeanss, idreservesnextmeanss, moments.ymeanind, 2];
				bondspread[iddebtnextylowstd, idreservesnextylowstdss, ylowstdind, 1]
				bondspread[iddebtnextylowstdss, idreservesnextylowstdss, ylowstdind, 2] ]
	# Make layers
	lineslayer=layer(debtframe, x=:Debt, y=:Spread, color=:ExoStates, Geom.line)
	pointslayer=layer(x=debtpoints, y=spreadpoints, Geom.point)
	# Make plot
	figure1=plot(lineslayer, pointslayer, Coord.cartesian(xmax=0.25*model.compuparams.debtmax+0.01, ymax=500),
					Guide.title("Spread given next period Debt"), latexstyle, Guide.colorkey(title="Legend"))
	# Save plot
	currdir=pwd()
	cd("..\\..\\figures\\drafts")
	draw(SVG("Spread-debt.svg", 4.5inch, 3.5inch), figure1)
	draw(PDF("Spread-debt.pdf", 4.5inch, 3.5inch), figure1)
	draw(PNG("Spread-debt.png", 4.5inch, 3.5inch), figure1)

	# 4.2 Spread vs Next Period Reserves 
	# Make dataframe for figure
	reservesframe=DataFrame(Reserves=repeat(0.25*model.grids.reserves, 4), 
						Spread=vec([bondspread[iddebtnextmean, :, moments.ymeanind, 1]' bondspread[iddebtnextmeanss, :, moments.ymeanind, 2]'
								bondspread[iddebtnextylowstd, :, ylowstdind, 1]' bondspread[iddebtnextylowstdss, :, ylowstdind, 2]' ]),
						DefProb=100*vec([defaultnextprob[iddebtnextmean, :, moments.ymeanind, 1]' defaultnextprob[iddebtnextmeanss, :, moments.ymeanind, 2]'
								defaultnextprob[iddebtnextylowstd, :, ylowstdind, 1]' defaultnextprob[iddebtnextylowstdss, :, ylowstdind, 2]' ]),
						NextDebt=0.25*vcat(view(debtpolicysmooth, moments.debtmeanind, :, moments.ymeanind, 1), 
						 					view(debtpolicysmooth, moments.debtmeanind, :, moments.ymeanind, 2),
		 									view(debtpolicysmooth, moments.debtmeanind, :, ylowstdind, 1), 
						 					view(debtpolicysmooth, moments.debtmeanind, :, ylowstdind, 2) ),
						ExoStates=repeat(["Yavg";"YavgPanic";"Ylow";"YlowPanic"], inner=[resnum]) )

	# Just bind x and y aesthetic for points
	reservespoints=0.25*[reservesnextmean reservesnextmeanss reservesnextylowstd reservesnextylowstdss]
	# Make layers
	lineslayer=layer(reservesframe, x=:Reserves, y=:Spread, color=:ExoStates, Geom.smooth)
	pointslayer=layer(x=reservespoints, y=spreadpoints, Geom.point)
	# Make plot
	figure2=plot(lineslayer, pointslayer, Coord.cartesian(xmax=0.25*model.compuparams.resmax+0.01, ymax=500),
					Guide.title("Spread given next period Reserves"), latexstyle, Guide.colorkey(title="Legend"))
	# Save plot
	draw(SVG("Spread-Reserves.svg", 4.5inch, 3.5inch), figure2)		
	draw(PDF("Spread-Reserves.pdf", 4.5inch, 3.5inch), figure2)		
	draw(PNG("Spread-Reserves.png", 4.5inch, 3.5inch), figure2)		

	# 4.3 Defalult prob next period given reserves
	# Make dataframe for figure

	# Just bind x and y aesthetic for points
	defaultpoints=100*[defaultnextprob[iddebtnextmean, idreservesnextmean, moments.ymeanind, 1];
				defaultnextprob[iddebtnextmeanss, idreservesnextmeanss, moments.ymeanind, 2];
				defaultnextprob[iddebtnextylowstd, idreservesnextylowstd, ylowstdind, 1]
				defaultnextprob[iddebtnextylowstdss, idreservesnextylowstdss, ylowstdind, 2] ]
	# Make layers
	lineslayer=layer(reservesframe, x=:Reserves, y=:DefProb, color=:ExoStates, Geom.smooth)
	pointslayer=layer(x=reservespoints, y=defaultpoints, Geom.point)
	legendguide=Guide.manual_color_key("Legend", [ "\$\\bar y \$, \$\\omega=1\$", "\$\\bar y \$, \$\\omega=0\$", 
													"\$\\bar y-1.6\$sigma \$, \$\\omega=1\$", "\$\\bar y-1.6\$sigma \$, \$\\omega=0\$"],
										[ RGB{N0f8}(0.996,0.263,0.396) , RGB{N0f8}(0.925,0.635,0.361),  
											RGB{N0f8}(0.137,0.431,0.678), RGB{N0f8}(0.482,0.969,0.882) ] )
	# Make plot
	figure3=plot(lineslayer, pointslayer, Coord.cartesian(xmax=0.15, ymax=5.0), legendguide, Guide.colorkey(title="Legend"),
					Guide.title("Default probability given next period Reserves"))
	# Save plot
	draw(SVG("DefProb-Reserves.svg", 4.5inch, 3.5inch), figure3)	
	draw(PDF("DefProb-Reserves.pdf", 4.5inch, 3.5inch), figure3)	
	draw(PNG("DefProb-Reserves.png", 4.5inch, 3.5inch), figure3)	

	# 3.4 Debt given Reserves THIS IS CRAZY
	# Plot next period debt given current reserves, zero temporary shock around mean debt mean level
	# Same dataframe
	# Make layers
	lineslayer=layer(reservesframe, x=:Reserves, y=:NextDebt, color=:ExoStates, Geom.smooth)
	pointslayer=layer(x=reservespoints, y=debtpoints, Geom.point)
	# Make plot
	figure4=plot(lineslayer, pointslayer, Guide.title("Next Debt given current Reserves"),
					Coord.cartesian(xmax=0.25*model.compuparams.resmax+0.01, ymax=0.25*model.compuparams.debtmax+0.01) )
	# Save plot
	draw(SVG("NextDebt-Reserves.svg", 4.5inch, 3.5inch), figure4)	
	draw(PDF("NextDebt-Reserves.pdf", 4.5inch, 3.5inch), figure4)	
	draw(PNG("NextDebt-Reserves.png", 4.5inch, 3.5inch), figure4)	

	# 3.5 Debt given current output

	# Find default thresholds
	defoutput1=findfirst(.!model.policies.default[moments.debtmeanind, moments.reservesmeanind, moments.mmeanind, :, 1])
	defoutput2=findfirst(.!model.policies.default[moments.debtmeanind, moments.reservesmeanind, moments.mmeanind, :, 2])
	
	# Make dataframe for figure
	outputframe=DataFrame(Output=repeat(model.grids.y, 2), 
						NextDebt=0.25*vec(view(debtpolicysmooth,moments.debtmeanind, moments.reservesmeanind, :, 1:2) ),
						NextReserves=0.25*vec(view(reservespolicysmooth, moments.debtmeanind, moments.reservesmeanind, :, 1:2 ) ),
						Sunspot=repeat(["Normal";"Panic"], inner=[ynum]),
						Default=vec(view(defaultpolicysmooth, moments.debtmeanind, moments.reservesmeanind, :, 1:2) ) )
	tempframe=by(outputframe,[:Sunspot, :Default], df->string(df[1,:Sunspot],'-', mean(df[1,:Default])>0.01 ? "Default" : "Repay" ) )
	outputframe=join(outputframe, tempframe, on=[:Sunspot, :Default], kind=:inner)

	insertcols!(outputframe, 6, DefaultText=outputframe[!, :Sunspot] )
	outputframe[!, :DefaultText] .= "Repay" 
	outputframe.DefaultText[(outputframe[!, :Default].>0.1)] .= "Default"
	outputframe.NextDebt[(outputframe[!, :DefaultText].=="Default")].=0.0
	# Just bind x and y aesthetic for points
	outputpoints=model.grids.y[ [moments.ymeanind moments.ymeanind ylowstdind ylowstdind] ]
	# Same figpointsy as in figure 4
	# Make layers
	lineslayer=layer(outputframe, x=:Output, y=:NextDebt, color=:x1, Geom.smooth)
	pointslayer=layer(x=outputpoints, y=debtpoints, Geom.point)
	axisdebtoutput=Coord.cartesian(xmin=model.grids.y[minimum(simulated.yind)], xmax=model.grids.y[maximum(simulated.yind)-3])
	# Make plot
	figure5=plot(lineslayer, pointslayer, axisdebtoutput, Guide.colorkey(title="Legend"), 
					Guide.title("Next Debt given current Output"))
	# Save plot
	draw(SVG("NextDebt-Output.svg", 4.5inch, 3.5inch), figure5)		
	draw(PDF("NextDebt-Output.pdf", 4.5inch, 3.5inch), figure5)		
	draw(PNG("NextDebt-Output.png", 4.5inch, 3.5inch), figure5)		

	sort!(outputframe, (:Sunspot))
	lineslayer2=layer(outputframe[1:50,:], x=:Output, y=:NextDebt, color=:DefaultText, Geom.smooth)
	figure5bis=plot(lineslayer2, pointslayer, axisdebtoutput, 
					Guide.title("Next Reserves given current Output"), Guide.colorkey(title="Legend"))
	draw(SVG("NextDebt-Output2.svg", 4.5inch, 3.5inch), figure5bis)		
	draw(PDF("NextDebt-Output2.pdf", 4.5inch, 3.5inch), figure5bis)		
	draw(PNG("NextDebt-Output2.png", 4.5inch, 3.5inch), figure5bis)		



	# 3.6 Reserves given current output
	# Same dataframe
	# Make layers
	lineslayer=layer(outputframe, x=:Output, y=:NextReserves, color=:x1, Geom.smooth)
	pointslayer=layer(x=outputpoints, y=reservespoints, Geom.point)
	axisresoutput=Coord.cartesian(xmin=model.grids.y[minimum(simulated.yind)], xmax=model.grids.y[maximum(simulated.yind)-3], ymin=0.0, ymax=0.16)
	figure6=plot(lineslayer, pointslayer, axisresoutput,
					Guide.title("Next Reserves given current Output"))
	draw(SVG("NextReserves-Output.svg", 4.5inch, 3.5inch), figure6)		
	draw(PDF("NextReserves-Output.pdf", 4.5inch, 3.5inch), figure6)		
	draw(PNG("NextReserves-Output.png", 4.5inch, 3.5inch), figure6)		

	sort!(outputframe, (:Sunspot))
	lineslayer2=layer(outputframe[1:50,:], x=:Output, y=:NextReserves, color=:DefaultText, Geom.smooth)
	figure6bis=plot(lineslayer2, pointslayer, axisresoutput, 
					Guide.title("Next Reserves given current Output"), Guide.colorkey(title="Legend"))
	draw(SVG("NextReserves-Output2.svg", 4.5inch, 3.5inch), figure6bis)		
	draw(PDF("NextReserves-Output2.pdf", 4.5inch, 3.5inch), figure6bis)		
	draw(PNG("NextReserves-Output2.png", 4.5inch, 3.5inch), figure6bis)		
	
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
	return outputframe
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