function mthresholds!(newthresholds::Array{Float64,1}, newthrespolicy::Array{Int64,1}, # Outputs
						interimthresholds::Array{Float64,1}, interimthresdebt::Array{Int64,1},
						consumption::Array{Float64,1}, expvalue::Array{Float64,1}, compuparams::ComputationParams, 
						mextremes::Array{Float64,1}, bbeta::Float64, debtminfres::Int64,
                        c0::Float64, c1::Float64, vdiff::Float64, tempthres::Float64 )			
    # initialize debt policy index
	interimthresdebt[:]=1:(compuparams.debtnum)

	for ifutdebt=(debtminfres+1):(compuparams.debtnum)  
        # Gamma=2 so the value of M can be solved
        c0=consumption[ifutdebt-1]
        c1=consumption[ifutdebt]
        
        if c1>c0
            vdiff=max(expvalue[ifutdebt-1]-expvalue[ifutdebt],1e-12)
            interimthresholds[ifutdebt]= -0.5*(c1+c0) + 0.5*((c1-c0)^2+4*(c1-c0)/(vdiff/(1-bbeta)))^0.5
        else # always better to choose iFdebt-1
            interimthresholds[ifutdebt]=compuparams.thrmin
            interimthresdebt[ifutdebt]=0
        end
    end
     # In principle there exist some m such that debtminfres-1 is preferred to debtminfres. It
     # has to be high enough to turn consumption positive at debtminfres-1 and get relatively 
     # close enough to consumption at debtminfres. However this m is very big, bigger than 
     # mextremes(end) by definition, so just to keep it bigger than all other thresholds
    interimthresholds[debtminfres]=mextremes[compuparams.mnum+1]
    interimthresholds[debtminfres]=1e-9+maximum(interimthresholds[debtminfres:(compuparams.debtnum)])
    # Check back dominance of thresholds
	for ifutdebt=(debtminfres+1):(compuparams.debtnum)
        if interimthresdebt[ifutdebt]==0 # debt not relevant
            continue
        end
        # If continue, this debt is relevant and all relevant debt levels smaller than iFdebt are set.
        if interimthresdebt[ifutdebt-1]!=0
            if interimthresholds[ifutdebt]>interimthresholds[ifutdebt-1] # Never choose (ifutdebt-1)
                interimthresdebt[ifutdebt-1]=0
                interimthresholds[ifutdebt-1]=compuparams.thrmin
            else
                # Previous is relevant, current is relevant and does NOT dominate previous.
                # Threshold stays the same. No need to check if it dominates lower debt levels
                continue
            end
        end
        
        # Here ifutdebt-1 was dominated. Check if ifutdebt also dominates ifutdebt-2 and all before until debtminfres. 
        # Only check for relevants
        c1=consumption[ifutdebt]


        for irevdebt=1:(ifutdebt-1-debtminfres) # we need to do this backwards
            # if debtminfres> ifutdebt-2 never enter the loop
            # Find threshols between ifutdebt and ifutdebt-1-i
            if interimthresdebt[ifutdebt-1-irevdebt]==0
                continue
            end
            c0=consumption[ifutdebt-1-irevdebt]
            if c1>c0
                vdiff=max( expvalue[ifutdebt-1-irevdebt]-expvalue[ifutdebt], 1e-12)
                tempthres=-0.5*(c1+c0)+0.5*((c1-c0)^2+4*(c1-c0)/vdiff*(1-bbeta))^0.5
                if tempthres>interimthresholds[ifutdebt-1-irevdebt]
                    # In this case ifutdebt dominates
                    if (ifutdebt-1-irevdebt)==debtminfres #This can happen, TempTthres can be larger.
                        # if this is the smallest debt then
                        interimthresholds[ifutdebt]=tempthres
                        # just to keep it bigger than ohter thresholds
                        interimthresholds[debtminfres]=max(tempthres, mextremes[compuparams.mnum+1])+1e-9
                        # We are keeping debtminfres always relevant
                    else
                        # Debt level not chosen
                        #ifutdebt-1-irevdebt
                        interimthresdebt[ifutdebt-1-irevdebt]=0
                        interimthresholds[ifutdebt-1-irevdebt]=compuparams.thrmin
                    end
                else # ifutdebt does not dominate ifutdebt-1-i
                    interimthresholds[ifutdebt]=tempthres
                    break # end the loop
                end
            else # consumption at ifutdebt-1-irevdebt higher than c0
                # always better to choose ifutdebt-1-irevdebt
                interimthresholds[ifutdebt]=compuparams.thrmin
                interimthresdebt[ifutdebt]=0
                break
            end
        end # Reverse debt loop end
    end # Future debt loop end
    #####
    #####
    # This part used to be outside. Make sure threshold and debt vectors are inbounds (mextremes)
    newthresnum::Int64=0
    # want to begin with higher debt hence reverse index
    for irevdebt=1:(compuparams.debtnum-debtminfres+1)
        if interimthresdebt[compuparams.debtnum+1-irevdebt]!=0
            if  interimthresholds[compuparams.debtnum+1-irevdebt]>mextremes[1]
                newthresnum+=1
                newthrespolicy[newthresnum]=compuparams.debtnum+1-irevdebt
                newthresholds[newthresnum]=interimthresholds[compuparams.debtnum+1-irevdebt]
                if interimthresholds[compuparams.debtnum+1-irevdebt]>mextremes[compuparams.mnum+1]
                    newthresholds[newthresnum]=mextremes[compuparams.mnum+1]
                    break
                end
            end
        end
    end
    return newthresnum
end # Function end
