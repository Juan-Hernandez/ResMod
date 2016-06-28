function mergethresholds!(thresholds::Array{Float64,1}, threspolicy::Array{Int64,2}, thresnum::Int64, # Outputs
							consexm::Array{Float64,2}, expvalue::Array{Float64,2}, bbeta::Float64, valtol::Float64, 
                            mextremes::Array{Float64,1},
							newthresholds::Array{Float64,1}, newthrespolicy::Array{Int64,1}, newthresnum::Int64, idftres::Int64, 
							interimbigthresholds::Array{Float64,1}, interimthrespolicy::Array{Int64,2}, bigthresnum::Int64,
							mergedsortidx::Array{Int64,1}, mergedthresholds::Array{Float64,1}, interimthresnum::Int64,
							debt1::Int64, debt2::Int64, res1::Int64, res2::Int64, tempindex::Int64,
							cons1::Float64, cons2::Float64, vdiff::Float64, tempthres::Float64)
	# THIS IS THE DIFFERENT PART. NEED TO MERGE THE DEBT,
	# THRESHOLD AND (IMPLICIT) RESERVES VECTORS IN ONE
    if thresnum==0 # easy to fill if none before
        thresnum=newthresnum
        thresholds[1:thresnum]=newthresholds[1:thresnum]
        threspolicy[1:thresnum ,1]=newthrespolicy[1:thresnum]
        threspolicy[1:thresnum, 2]=idftres
    else # need to merge
        bigthresnum=thresnum+newthresnum
        interimbigthresholds[1:bigthresnum]=[thresholds[1:thresnum]; newthresholds[1:newthresnum]]
        sortperm!(sub(mergedsortidx, 1:bigthresnum), interimbigthresholds[1:bigthresnum], alg=QuickSort)
        interimbigthresholds[1:bigthresnum]=interimbigthresholds[ mergedsortidx[1:bigthresnum] ]
        interimthresnum=0
        # Look at the first interval. Has to compare
        # both inidces are (1)
        debt1=threspolicy[1 ,1]
        res1=threspolicy[1 ,2]
        cons1=consexm[debt1,res1]
        debt2=newthrespolicy[1]
        res2=idftres
        cons2=consexm[debt2,res2]
        vdiff=0.0
        vdiff+=-expvalue[debt1,res1]
        vdiff+=expvalue[debt2,res2]
        # Find threshold and direction
        # Below interimbig[1] we are comparing the first two options on each vector
        if abs(vdiff)<1e-6*(1-bbeta)*valtol # if equal, look consumption
            mergedthresholds[1]=interimbigthresholds[1]
            interimthresnum+=1
            if cons1>cons2
                interimthrespolicy[1, :]=[debt1 res1]
            else
                interimthrespolicy[1, :]=[debt2 res2]
            end
        elseif vdiff>9.0*1e-7*(1-bbeta)*valtol
            if cons2>=cons1 # better present and future
                mergedthresholds[1]=interimbigthresholds[1]
                interimthresnum+=1
                interimthrespolicy[1, :]=[debt2 res2]
            else # Find threshold
                tempthres=-0.5*(cons1+cons2)+0.5*((cons1-cons2)^2+4*(cons1-cons2)/vdiff*(1-bbeta))^0.5
                if tempthres>interimbigthresholds[1] # 1 is preferred until interimbig[1]
                    mergedthresholds[1]=interimbigthresholds[1]
                    interimthresnum+=1
                    interimthrespolicy[1, :]=[debt1 res1]
                elseif tempthres<=mextremes[1] # 2 is always preferred
                    mergedthresholds[1]=interimbigthresholds[1]
                    interimthresnum+=1
                    interimthrespolicy[1, :]=[debt2 res2]
                else # new threshold: 1 is better until it
                    mergedthresholds[1]=tempthres
                    interimthrespolicy[1, :]=[debt1 res1]
                    mergedthresholds[2]=interimbigthresholds[1]
                    interimthrespolicy[2, :]=[debt2 res2]
                    interimthresnum+=2
                end
            end
        else # vdiff is negative Vfut1>Vfut2
            if cons1>=cons2 # better present and future
                mergedthresholds[ ]=interimbigthresholds[1]
                interimthrespolicy[1, :]=[debt1 res1]
                interimthresnum+=1
            else # Find threshold
                tempthres=-0.5*(cons1+cons2)+0.5*((cons1-cons2)^2+4*(cons1-cons2)/vdiff*(1-bbeta))^0.5;
                if tempthres>interimbigthresholds[1] # 2 is preferred
                    mergedthresholds[1]=interimbigthresholds[1]
                    interimthrespolicy[1, :]=[debt2 res2]
                    interimthresnum+=1
                elseif tempthres<=mextremes[1] # 1 is preferred
                    mergedthresholds[1]=interimbigthresholds[1]
                    interimthrespolicy[1, :]=[debt1 res1]
                    interimthresnum+=1
                else # new threshold 2 is better until it
                    mergedthresholds[1]=tempthres
                    interimthrespolicy[1, :]=[debt2 res2]
                    mergedthresholds[2]=interimbigthresholds[1]
                    interimthrespolicy[2, :]=[debt1 res1]
                    interimthresnum+=2
                end
            end
        end

        # Same as before but taking care of the previous interval. 
        # Also uses the indices mergedsortidx. 
        for idbig=2:bigthresnum
            # If too close, discard
            if (interimbigthresholds[idbig]-interimbigthresholds[idbig-1])<1e-14 
                continue
            end
            if mergedsortidx[idbig]>thresnum # this is from 2
                tempindex=findfirst(x->(x<(thresnum+1)), mergedsortidx[idbig:bigthresnum])
                debt1=threspolicy[mergedsortidx[idbig-1+tempindex], 1]
                res1=threspolicy[mergedsortidx[idbig-1+tempindex], 2]
                debt2=newthrespolicy[mergedsortidx[idbig]-thresnum]
                res2=idftres
            else # index is from 1
                tempindex=findfirst(x->(x>thresnum), mergedsortidx[idbig:bigthresnum])
                debt1=threspolicy[mergedsortidx[idbig], 1]
                res1=threspolicy[mergedsortidx[idbig], 2]
                debt2=newthrespolicy[mergedsortidx[idbig-1+tempindex]-thresnum]
                res2=idftres
            end
            cons1=consexm[debt1, res1]
            cons2=consexm[debt2, res2]
            vdiff=0.0
            vdiff+=-expvalue[debt1,res1]
            vdiff+=expvalue[debt2,res2]
            # Find threshold and direction
            if abs(vdiff)<1e-6*(1-bbeta)*valtol # if equal, look consumption
                if cons1>cons2 # append/extend 1
                    if (interimthrespolicy[interimthresnum,1]==debt1) && (interimthrespolicy[interimthresnum,2]==res1)
                        # extend interval
                        mergedthresholds[interimthresnum]=interimbigthresholds[idbig]
                    else # append new interval
                        interimthresnum+=1 # index increased before
                        mergedthresholds[interimthresnum]=interimbigthresholds[idbig]
                        interimthrespolicy[interimthresnum, :]=[debt1 res1]
                    end
                else # append/extend 2
                    if (interimthrespolicy[interimthresnum,1]==debt2) && (interimthrespolicy[interimthresnum,2]==res2)
                        # extend interval
                        mergedthresholds[interimthresnum]=interimbigthresholds[idbig]
                    else # append new interval
                        interimthresnum+=1 # index increased before
                        mergedthresholds[interimthresnum]=interimbigthresholds[idbig]
                        interimthrespolicy[interimthresnum, :]=[debt2 res2]
                    end
                end
            elseif vdiff>9.0*1e-7*(1-bbeta)*valtol
                if cons2>=cons1 # better present and future: append/extend 2
                    if (interimthrespolicy[interimthresnum,1]==debt2) && (interimthrespolicy[interimthresnum,2]==res2)
                        # extend interval
                        mergedthresholds[interimthresnum]=interimbigthresholds[idbig]
                    else # append new interval
                        interimthresnum+=1 # index increased before
                        mergedthresholds[interimthresnum]=interimbigthresholds[idbig]
                        interimthrespolicy[interimthresnum, :]=[debt2 res2]
                    end
                else # Find threshold
                    tempthres=-0.5*(cons1+cons2)+0.5*((cons1-cons2)^2+4*(cons1-cons2)/vdiff*(1-bbeta))^0.5;
                    if tempthres>interimbigthresholds[idbig] # 1 is preferred append/extend 1
                        if (interimthrespolicy[interimthresnum,1]==debt1) && (interimthrespolicy[interimthresnum,2]==res1)
                            # extend interval
                            mergedthresholds[interimthresnum]=interimbigthresholds[idbig]
                        else # append new interval
                            interimthresnum+=1 # index increased before
                            mergedthresholds[interimthresnum]=interimbigthresholds[idbig]
                            interimthrespolicy[interimthresnum, :]=[debt1 res1]
                        end                        
                    elseif tempthres<=interimbigthresholds[idbig-1] # 2 is preferred  append/extend 2
                        if (interimthrespolicy[interimthresnum,1]==debt2) && (interimthrespolicy[interimthresnum,2]==res2)
                            # extend interval
                            mergedthresholds[interimthresnum]=interimbigthresholds[idbig]
                        else # append new interval
                            interimthresnum+=1 # index increased before
                            mergedthresholds[interimthresnum]=interimbigthresholds[idbig]
                            interimthrespolicy[interimthresnum, :]=[debt2 res2]
                        end                        
                    else # new threshold: 1 is better until tempthres, threshold 2 until interimbigthresholds
                        if (interimthrespolicy[interimthresnum,1]==debt1) && (interimthrespolicy[interimthresnum,2]==res1) # append/extend 1 append 2 for sure
                            # extend first interval
                            mergedthresholds[interimthresnum]=tempthres
                            # append second interval
                            interimthresnum+=1 # index increased before
                            mergedthresholds[interimthresnum]=interimbigthresholds[idbig]
                            interimthrespolicy[interimthresnum, :]=[debt2 res2]
                        else # append 2 new intervals
                            interimthresnum+=1 # index increased before
                            mergedthresholds[interimthresnum]=tempthres
                            interimthrespolicy[interimthresnum, :]=[debt1 res1]
                            interimthresnum+=1 # index increased before again                            
                            mergedthresholds[interimthresnum]=interimbigthresholds[idbig]
                            interimthrespolicy[interimthresnum, :]=[debt2 res2]
                        end
                    end
                end
            else # vdiff is negative Vfut1>Vfut2
                if cons1>=cons2 # better present and future: append/extend 1
                    if (interimthrespolicy[interimthresnum,1]==debt1) && (interimthrespolicy[interimthresnum,2]==res1)
                        # extend interval
                        mergedthresholds[interimthresnum]=interimbigthresholds[idbig]
                    else # append new interval
                        interimthresnum+=1 # index increased before
                        mergedthresholds[interimthresnum]=interimbigthresholds[idbig]
                        interimthrespolicy[interimthresnum, :]=[debt1 res1]
                    end                        
                else # Find threshold
                    tempthres=-0.5*(cons1+cons2)+0.5*((cons1-cons2)^2+4*(cons1-cons2)/vdiff*(1-bbeta))^0.5;
                    if tempthres>interimbigthresholds[idbig] # 2 is preferred: append/extend 2
                        if (interimthrespolicy[interimthresnum,1]==debt2) && (interimthrespolicy[interimthresnum,2]==res2)
                            # extend interval
                            mergedthresholds[interimthresnum]=interimbigthresholds[idbig]
                        else # append new interval
                            interimthresnum+=1 # index increased before
                            mergedthresholds[interimthresnum]=interimbigthresholds[idbig]
                            interimthrespolicy[interimthresnum, :]=[debt2 res2]
                        end
                    elseif tempthres<=interimbigthresholds[idbig-1] # 1 is preferred: append/extend 1
                        if (interimthrespolicy[interimthresnum,1]==debt1) && (interimthrespolicy[interimthresnum,2]==res1)
                            # extend interval
                            mergedthresholds[interimthresnum]=interimbigthresholds[idbig]
                        else # append new interval
                            interimthresnum+=1 # index increased before
                            mergedthresholds[interimthresnum]=interimbigthresholds[idbig]
                            interimthrespolicy[interimthresnum, :]=[debt1 res1]
                        end                        
                    else # new threshold: 2 is better until it. Then 1.
                        if (interimthrespolicy[interimthresnum,1]==debt2) && (interimthrespolicy[interimthresnum,2]==res2)
                            # extend interval
                            mergedthresholds[interimthresnum]=tempthres
                            # append first interval
                            interimthresnum+=1 # index increased before
                            mergedthresholds[interimthresnum]=interimbigthresholds[idbig]
                            interimthrespolicy[interimthresnum, :]=[debt1 res1]
                        else # append 2 new intervals
                            interimthresnum+=1 # index increased before
                            mergedthresholds[interimthresnum]=tempthres
                            interimthrespolicy[interimthresnum, :]=[debt2 res2]
                            interimthresnum+=1 # index increased before again                            
                            mergedthresholds[interimthresnum]=interimbigthresholds[idbig]
                            interimthrespolicy[interimthresnum, :]=[debt1 res1]
                        end
                    end
                end
            end
        end # loop over big thresholds ends
        # Finished merging. interimbigthresholds has been run over
        # Fill thresholds, threspolicy[ ,1] and threspolicy[ ,2]
        thresnum=interimthresnum;
        thresholds[1:thresnum]=mergedthresholds[1:thresnum]
        thresholds[thresnum+1]=NaN
        threspolicy[1:thresnum,:]=interimthrespolicy[1:thresnum,:]
        # To make sure we stay in bounds
        threspolicy[thresnum+1, 1]=-1
        threspolicy[thresnum+1, 2]=-1
    end 
    return thresnum
end # Function end