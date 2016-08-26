function defaultthresholds!(thresholds::Array{Float64,1}, threspolicy::Array{Int64,2}, thresnum::Int64, thresdefault::BitArray{1}, # Outputs
                            expvalue::Array{Float64,2}, valuedefault::Float64, defaultreserves::Int64, consexm::Array{Float64,2},
                            valtol::Float64, thrmin::Float64, ggamma::Int64, bbeta::Float64, mint1::Float64,
							interimnewthresholds::Array{Float64,1}, interimthrespolicy::Array{Int64,2})
                            
    # 0. Allocate temporary variables
    # 0.1 Temporary Integers                        
    interimthresnum::Int64=0
    firstdebtnodef::Int64=0
    toappend::Int64=0
    # 0.2 Temporary Floats
    tempthres::Float64=0.0
    lowthres::Float64=0.0
    vdiff::Float64=0.0

    # 1. Loop over current thresholds
    thresdefault[1:thresnum]=false # Initially no default.
    for idthres=1:thresnum # Run over current thresholds
        vdiff=valuedefault-expvalue[threspolicy[idthres, 1], threspolicy[idthres, 2]]
        if vdiff>-1e-6*(1-bbeta)*valtol # as flow utility is negative, default is sure
            tempthres=-thrmin # same as +ThrMax, always below, always default
        else # Find threshold (this can solve for all gamma)
            tempthres=(vdiff*(1-ggamma)/(1-bbeta))^(1/(1-ggamma))-consexm[threspolicy[idthres, 1], threspolicy[idthres, 2]]
            # Normalized Utility
        end
        if idthres==1
            lowthres=mint1
        else
            lowthres=thresholds[idthres-1]
        end
        if tempthres<lowthres
            # Note that value function has to be increasing on m. 
            # If no default over an entire interval, there wont be 
            # defaults for bigger M's. Hence we should append the rest of the vectors
            toappend=thresnum-idthres+1
            interimnewthresholds[(interimthresnum+1):(interimthresnum+toappend)]=thresholds[idthres:thresnum]
            interimthrespolicy[(interimthresnum+1):(interimthresnum+toappend), 1]=threspolicy[idthres:thresnum, 1]
            interimthrespolicy[(interimthresnum+1):(interimthresnum+toappend), 2]=threspolicy[idthres:thresnum, 2]
            thresdefault[(interimthresnum+1):(interimthresnum+toappend)]=false # No default again
            interimthresnum+=toappend
            firstdebtnodef=threspolicy[idthres, 1]
            break
        elseif tempthres>thresholds[idthres] # Default on whole interval
            # Again value function increasing on m. If default over entire interval, default 
            # will happen for sure for all smaller M's. Hence only 1 relevant threshold is left.
            interimthresnum=1
            interimnewthresholds[interimthresnum]=thresholds[idthres]
            interimthrespolicy[interimthresnum, 1]=1 # Index for zero debt
            # Index for optimal reserves given default
            interimthrespolicy[interimthresnum, 2]=defaultreserves
            thresdefault[interimthresnum]=true # Default on interval
        else # Theshold in the middle. Again no default above threshold, append the rest.
            # Default always below, exit loop
            interimthresnum=1
            interimnewthresholds[interimthresnum]=tempthres
            interimthrespolicy[interimthresnum, 1]=1 # Index for zero debt
            # Index for optimal reserves given default
            interimthrespolicy[interimthresnum, 2]=defaultreserves
            thresdefault[interimthresnum]=true # Default on first interval
            # Append the other intervals, no default
            toappend=thresnum-idthres+1
            interimnewthresholds[(interimthresnum+1):(interimthresnum+toappend)]=thresholds[idthres:thresnum]
            interimthrespolicy[(interimthresnum+1):(interimthresnum+toappend), 1]=threspolicy[idthres:thresnum, 1]
            interimthrespolicy[(interimthresnum+1):(interimthresnum+toappend), 2]=threspolicy[idthres:thresnum, 2]
            thresdefault[(interimthresnum+1):(interimthresnum+toappend)]=false # No default again
            interimthresnum+=toappend
            firstdebtnodef=threspolicy[idthres, 1]
            break
        end
    end
    # Now have thresholds and actions including default.
    thresnum=interimthresnum
    thresholds[1:thresnum]=interimnewthresholds[1:interimthresnum]
    threspolicy[1:thresnum, :]=interimthrespolicy[1:interimthresnum, :]
    return (thresnum, firstdebtnodef)
end # function end