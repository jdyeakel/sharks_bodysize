function modality_analysis(pdist,toothlength)


    lmax=findlocalmaxima(pdist); #
    lmin=findlocalminima(pdist);
    #choose highest probability peaks
    peakprobs = pdist[lmax];
    troughprobs = pdist[lmin];
    sortlmax = lmax[sortperm(peakprobs)];
    sortlmin = lmin[sortperm(troughprobs)];
    sortedpeaks = pdist[sortlmax];
    sortedtroughs = pdist[sortlmin];
    maxpeaka = last(sortedpeaks);
    sortedtooth = toothlength[lmax[sortperm(peakprobs)]];
    maxtooth = last(sortedtooth);
    secondpeak = 0.0;
    secondtooth = 0.0;
    troughtooth = 0.0;
    troughmin = 0.0;
    troughpos = Array{Int64}(undef,0) .* 0;
    troughposbackup = Array{Int64}(undef,0) .* 0;
    trougha1 = Array{Float64}(undef,0) .* 0;
    trougha2 = Array{Float64}(undef,0) .* 0;
    trougha3 = Array{Float64}(undef,0) .* 0;
    trougha4 = Array{Float64}(undef,0) .* 0;
    peakprop = (Array{Float64}(undef,2) .* 0) .- 1;
    if length(sortlmin) > 1
        #Which troughs are in the middle?
        troughpos1 = sortlmin[findall(x->(x < sortlmax[end] && x > sortlmax[end-1]),sortlmin)];
        troughpos2 = sortlmin[findall(x->(x > sortlmax[end] && x < sortlmax[end-1]),sortlmin)];
        troughpos = [troughpos1;troughpos2];
        
        if length(sortlmax) > 2
            troughpos3 = sortlmin[findall(x->(x < sortlmax[end] && x > sortlmax[end-2]),sortlmin)];
            troughpos4 = sortlmin[findall(x->(x > sortlmax[end] && x < sortlmax[end-2]),sortlmin)];
            troughposbackup = [troughpos3;troughpos4];
        end

        trougha1 = pdist[troughpos][findmin(pdist[troughpos])[2]];
        if length(troughposbackup) > 0
            trougha2 = pdist[troughposbackup][findmin(pdist[troughposbackup])[2]];
        end
        
        #need to determine second peak significance
        if length(sortedpeaks) > 1
            #threshold of 25% size of first peak
            # if ((sortedpeaks[length(sortedpeaks)-1] / maxpeaka) > 0.25) && (maxpeaka - (trougha[1]/maxpeaka) > 0.25)
            peakprop[1] = (sortedpeaks[end-1] - trougha1) / (maxpeaka - trougha1);
            if length(troughposbackup) > 0
                peakprop[2] = (sortedpeaks[end-2] - trougha2) / (maxpeaka - trougha2);
            end
            peakpos = 0;
            if peakprop[1] > peakprop[2]
                peakpos = 1
            end
            if peakprop[2] > peakprop[1]
                peakpos = 2
            end
            
            if peakprop[peakpos] > 0.10
                # if length(sortedpeaks) > peakpos
                    secondpeak = sortedpeaks[end-peakpos];
                    secondtooth = sortedtooth[end-peakpos];
                # else 
                #     println(i)
                # end
            end
        end
        
        if length(troughposbackup) > 0
            troughposmin = [troughpos[1], troughposbackup[1]][peakpos];
        else
            troughposmin = troughpos[1];
        end
        
        troughmin = [trougha1,trougha2][peakpos];
        troughtooth = toothlength[troughposmin];
    else
        secondpeak = 0.0;
        secondtooth = 0.0;
        troughtooth = 0.0;
        troughmin = 0.0;
    end
        
    return secondpeak, secondtooth, troughtooth, troughmin, maxtooth
end




