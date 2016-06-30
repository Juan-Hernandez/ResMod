function makegrids(computationparameters::ComputationParams,econparams::EconParams)
	# Decompose autocorrelation and variance of output in y-m terms.
	
	ysigma::Float64=(econparams.logoutputsigma^2-computationparameters.msigma^2)^0.5
	yrho::Float64=econparams.logoutputrho*(econparams.logoutputsigma/ysigma)^2
	yerrsigma::Float64=ysigma*(1-yrho^2)^0.5
	# Output Grid using Rouwenhorst

	(lygrid::Array{Float64,1},ytrans::Array{Float64,2}) = rouwenhorst(computationparameters.ynum,0,yrho,yerrsigma);
	y::Array{Float64,1}=exp(lygrid)
	# ergodic dist of y
	yergodic::Array{Float64,1}=limitdist(ytrans)
	# Output grid when default
	ydefault::Array{Float64,1}=min(y,y.*(1-econparams.defcost1-econparams.defcost2*y))
	#Make sure YdefGrid is nondecreasing:
	for i=2:computationparameters.ynum
		ydefault[i]=maximum(ydefault[(i-1):i])
	end
	
	# Temporary shock grid for integration:
	# Define intervals and use midpoints
	mextremes::Array{Float64,1}=linspace(-1.6*computationparameters.msigma,1.6*computationparameters.msigma,computationparameters.mnum+1)
	mmidpoints::Array{Float64,1}=midpoints(mextremes)
	# Mass at each interval

	# Use erf function \Phi(z)=0.5+0.5*erf(z/sqrt(2))

	####
	mmass::Array{Float64,1}=diff(0.5+0.5*erf(mextremes/(sqrt(2)*computationparameters.msigma)))
	mmass=mmass/sum(mmass)

	# Debt Grid (vector): My debt is positive
	debt::Array{Float64,1}=linspace(computationparameters.debtmin, computationparameters.debtmax, computationparameters.debtnum)
	debtmaxss=Array(Int64, computationparameters.debtnum)
	for idebt=1:computationparameters.debtnum
		debtmaxss[idebt]=findlast( x->(x<=(1-econparams.llambda)*debt[idebt]), debt )
	end
	# Reserves Grid (vector):
	reserves::Array{Float64,1}=linspace(computationparameters.resmin, computationparameters.resmax, computationparameters.resnum)

	# Regime transition:
	regimetrans::Array{Float64,2}= [1.0-1.0/econparams.panicfrequency 1.0/econparams.panicfrequency; 1.0/econparams.panicduration 1.0-1.0/econparams.panicduration]
	# Remember, first debt, second reserves, third Mshock, fourth output, fifth sudden stop
	grids=ModelGrids( y, ytrans, yergodic, ydefault, mmidpoints, mextremes, mmass, debt, debtmaxss, reserves, regimetrans)
	return grids
end

function limitdist(originalmat::Matrix)
""" Obtain the stationary probability distribution
vvector p of an irreducible, recurrent Markov
chain by state reduction. transmat is the transition
probabilities matrix of a discrete-time Markov
chain or the generator matrix Q. """
	transmat=copy(originalmat)
	(rownum,colnum)=size(transmat)
	n=rownum::Int64
	p=zeros(n)::Array{Float64,1}

	while n>1
		n1=n-1::Int64
		s=sum(transmat[n,1:n1])::Float64
		transmat[1:n1,n]=transmat[1:n1,n]/s
		n2=n1;
		while n2>0
			transmat[1:n1,n2]=transmat[1:n1,n2]+transmat[1:n1,n]*transmat[n,n2]
			n2=n2-1
		end
		n=n-1
	end
	# backtracking
	p[1]=1
	j=2::Int64
	while j<=rownum
		j1=j-1::Int64
		p[j]=sum(p[1:j1].*(transmat[1:j1,j]))
		j=j+1
	end
	p=p/sum(p)
	return p 
end