function toothdist_analysis(toothdrop,toothlength)


    jden = (toothdrop[:,1]/sum(toothdrop[:,1]));
    aden = (toothdrop[:,2]/sum(toothdrop[:,2]));

    #MEASURE SOMETHING
    #MEANS
    mvj = dot(toothlength,jden);
    mva = dot(toothlength,aden);

    #MEDIANS
    medianj = toothlength[findall(x->x>0.5,cumsum(jden))[1]];
    mediana = toothlength[findall(x->x>0.5,cumsum(aden))[1]];

    #Quartiles
    quartile25j = toothlength[findall(x->x>0.25,cumsum(jden))[1]];
    quartile25a = toothlength[findall(x->x>0.25,cumsum(aden))[1]];

    quartile75j = toothlength[findall(x->x>0.75,cumsum(jden))[1]];
    quartile75a = toothlength[findall(x->x>0.75,cumsum(aden))[1]];

    
    # meanjuv[r,temp_pos,dist_pos,sigtau_pos,tau_pos] = mvj;
    # meanadult[r,temp_pos,dist_pos,sigtau_pos,tau_pos] = mva;
    
    #NOTE: collect later
    # meanjuv[r,sigtau_pos,tau_pos] = mvj;
    # meanadult[r,sigtau_pos,tau_pos] = mva;
    
    #calculate variance
    meanofxsquared = dot(toothlength.^2,(toothdrop[:,1]/sum(toothdrop[:,1])));
    squareofmeanofx = mvj^2;
    varj = meanofxsquared - squareofmeanofx;
    # varjuv[r,temp_pos,dist_pos,sigtau_pos,tau_pos] = varj;
    #NOTE: collect later
    # varjuv[r,sigtau_pos,tau_pos] = varj;
    
    meanofxsquared = dot(toothlength.^2,(toothdrop[:,2]/sum(toothdrop[:,2])));
    squareofmeanofx = mva^2;
    vara = meanofxsquared - squareofmeanofx;
    # varadult[r,temp_pos,dist_pos,sigtau_pos,tau_pos] = vara;
    #NOTE: collect later
    # varadult[r,sigtau_pos,tau_pos] = vara;

    pj = (toothdrop[:,1]/sum(toothdrop[:,1]));
    pa = (toothdrop[:,2]/sum(toothdrop[:,2]));

    #NOTE: place modality analysis as function
    secondpeakj, secondtoothj, troughtoothj, troughminj, maxtoothj = modality_analysis(pj,toothlength);
    secondpeaka, secondtootha, troughtootha, troughmina, maxtootha = modality_analysis(pa,toothlength);

    

    # lmax=findlocalmaxima(pa); #
    # lmin=findlocalminima(pa);
    # #choose highest probability peaks
    # peakprobs = pa[lmax];
    # troughprobs = pa[lmin];
    # sortlmax = lmax[sortperm(peakprobs)];
    # sortlmin = lmin[sortperm(troughprobs)];
    # sortedpeaks = pa[sortlmax];
    # sortedtroughs = pa[sortlmin];
    # maxpeaka = last(sortedpeaks);
    # sortedtooth = toothlength[lmax[sortperm(peakprobs)]];
    # maxtootha = last(sortedtooth);
    # secondpeaka = 0.0;
    # secondtootha = 0.0;
    # troughtootha = 0.0;
    # troughmina = 0.0;
    # troughpos = Array{Int64}(undef,0) .* 0;
    # troughposbackup = Array{Int64}(undef,0) .* 0;
    # trougha1 = Array{Float64}(undef,0) .* 0;
    # trougha2 = Array{Float64}(undef,0) .* 0;
    # trougha3 = Array{Float64}(undef,0) .* 0;
    # trougha4 = Array{Float64}(undef,0) .* 0;
    # peakprop = (Array{Float64}(undef,2) .* 0) .- 1;
    # if length(sortlmin) > 1
    #     #Which troughs are in the middle?
    #     troughpos1 = sortlmin[findall(x->(x < sortlmax[end] && x > sortlmax[end-1]),sortlmin)];
    #     troughpos2 = sortlmin[findall(x->(x > sortlmax[end] && x < sortlmax[end-1]),sortlmin)];
    #     troughpos = [troughpos1;troughpos2];
        
    #     if length(sortlmax) > 2
    #         troughpos3 = sortlmin[findall(x->(x < sortlmax[end] && x > sortlmax[end-2]),sortlmin)];
    #         troughpos4 = sortlmin[findall(x->(x > sortlmax[end] && x < sortlmax[end-2]),sortlmin)];
    #         troughposbackup = [troughpos3;troughpos4];
    #     end

    #     trougha1 = pa[troughpos][findmin(pa[troughpos])[2]];
    #     if length(troughposbackup) > 0
    #         trougha2 = pa[troughposbackup][findmin(pa[troughposbackup])[2]];
    #     end
        
    #     #need to determine second peak significance
    #     if length(sortedpeaks) > 1
    #         #threshold of 25% size of first peak
    #         # if ((sortedpeaks[length(sortedpeaks)-1] / maxpeaka) > 0.25) && (maxpeaka - (trougha[1]/maxpeaka) > 0.25)
    #         peakprop[1] = (sortedpeaks[end-1] - trougha1) / (maxpeaka - trougha1);
    #         if length(troughposbackup) > 0
    #             peakprop[2] = (sortedpeaks[end-2] - trougha2) / (maxpeaka - trougha2);
    #         end
    #         peakpos = 0;
    #         if peakprop[1] > peakprop[2]
    #             peakpos = 1
    #         end
    #         if peakprop[2] > peakprop[1]
    #             peakpos = 2
    #         end
            
    #         if peakprop[peakpos] > 0.10
    #             # if length(sortedpeaks) > peakpos
    #                 secondpeaka = sortedpeaks[end-peakpos];
    #                 secondtootha = sortedtooth[end-peakpos];
    #             # else 
    #             #     println(i)
    #             # end
    #         end
    #     end
        
    #     if length(troughposbackup) > 0
    #         troughposmin = [troughpos[1], troughposbackup[1]][peakpos];
    #     else
    #         troughposmin = troughpos[1];
    #     end
        
    #     troughmina = [trougha1,trougha2][peakpos];
    #     troughtootha = toothlength[troughposmin];
    # else
    #     secondpeaka = 0.0;
    #     secondtootha = 0.0;
    #     troughtootha = 0.0;
    #     troughmina = 0.0;
    # end
        
        
    
  
    
    
    # lmax=findlocalmaxima(pj); #
    # lmin=findlocalminima(pj);
    # #choose highest probability peaks
    # peakprobs = pj[lmax];
    # troughprobs = pj[lmin];
    # sortlmax = lmax[sortperm(peakprobs)];
    # sortlmin = lmin[sortperm(troughprobs)];
    # sortedpeaks = pj[sortlmax];
    # sortedtroughs = pj[sortlmin];
    # maxpeakj = last(sortedpeaks);
    # sortedtooth = toothlength[lmax[sortperm(peakprobs)]];
    # maxtoothj = last(sortedtooth);
    # secondpeakj = 0.0;
    # secondtoothj = 0.0;
    # troughtoothj = 0.0;
    # troughminj = 0.0;
    # troughpos = Array{Int64}(undef,0) .* 0;
    # troughposbackup = Array{Int64}(undef,0) .* 0;
    # troughj1 = Array{Float64}(undef,0) .* 0;
    # troughj2 = Array{Float64}(undef,0) .* 0;
    # troughj3 = Array{Float64}(undef,0) .* 0;
    # troughj4 = Array{Float64}(undef,0) .* 0;
    # peakprop = (Array{Float64}(undef,2) .* 0) .- 1;
    # if length(sortlmin) > 1
    #     #Which troughs are in the middle?
    #     troughpos1 = sortlmin[findall(x->(x < sortlmax[end] && x > sortlmax[end-1]),sortlmin)];
    #     troughpos2 = sortlmin[findall(x->(x > sortlmax[end] && x < sortlmax[end-1]),sortlmin)];
    #     troughpos = [troughpos1;troughpos2];
        
    #     if length(sortlmax) > 2
    #         troughpos3 = sortlmin[findall(x->(x < sortlmax[end] && x > sortlmax[end-2]),sortlmin)];
    #         troughpos4 = sortlmin[findall(x->(x > sortlmax[end] && x < sortlmax[end-2]),sortlmin)];
    #         troughposbackup = [troughpos3;troughpos4];
    #     end

    #     troughj1 = pj[troughpos][findmin(pj[troughpos])[2]];
    #     if length(troughposbackup) > 0
    #         troughj2 = pj[troughposbackup][findmin(pj[troughposbackup])[2]];
    #     end
        
    #     #need to determine second peak significance
    #     if length(sortedpeaks) > 1
    #         #threshold of 25% size of first peak
    #         # if ((sortedpeaks[length(sortedpeaks)-1] / maxpeaka) > 0.25) && (maxpeaka - (trougha[1]/maxpeaka) > 0.25)
    #         peakprop[1] = (sortedpeaks[end-1] - troughj1) / (maxpeakj - troughj1);
    #         if length(troughposbackup) > 0
    #             peakprop[2] = (sortedpeaks[end-2] - troughj2) / (maxpeakj - troughj2);
    #         end
    #         peakpos = 0;
    #         if peakprop[1] > peakprop[2]
    #             peakpos = 1
    #         end
    #         if peakprop[2] > peakprop[1]
    #             peakpos = 2
    #         end
            
    #         if peakprop[peakpos] > 0.10
    #             # if length(sortedpeaks) > peakpos
    #                 secondpeakj = sortedpeaks[end-peakpos];
    #                 secondtoothj = sortedtooth[end-peakpos];
    #             # else 
    #             #     println(i)
    #             # end
    #         end
    #     end
        
    #     if length(troughposbackup) > 0
    #         troughposmin = [troughpos[1], troughposbackup[1]][peakpos];
    #     else
    #         troughposmin = troughpos[1];
    #     end
        
    #     troughminj = [troughj1,troughj2][peakpos];
    #     troughtoothj = toothlength[troughposmin];
    # else
    #     secondpeakj = 0.0;
    #     secondtoothj = 0.0;
    #     troughtoothj = 0.0;
    #     troughminj = 0.0;
    # end
    




    # pl = lineplot(toothlength,pa);
    # lineplot!(pl,repeat([maxtootha],outer=2),[0,maxpeaka],color=:red)
    # lineplot!(pl,repeat([secondtootha],outer=2),[0,secondpeaka],color=:red)
    # lineplot!(pl,repeat([troughtootha],outer=2),[0,troughmina],color=:blue)
    # 
    # 
    # 
    # pl = lineplot(toothlength,pj);
    # lineplot!(pl,repeat([maxtoothj],outer=2),[0,maxpeakj],color=:red)
    # lineplot!(pl,repeat([secondtoothj],outer=2),[0,secondpeakj],color=:red)
    # lineplot!(pl,repeat([troughtoothj],outer=2),[0,troughminj],color=:blue)
    # 
        
    
    #Record presence of multiple modes
    if secondpeakj == 0.0
        # peakjuv[r,temp_pos,dist_pos,sigtau_pos,tau_pos] = 0;
        # peakjuvquant[r,temp_pos,dist_pos,sigtau_pos,tau_pos] = 0;
        #NOTE: collect later
        # peakjuv[r,sigtau_pos,tau_pos] = 0;
        # peakjuvquant[r,sigtau_pos,tau_pos] = 0;

        peakjuvbin = 0;
        peakjuvdist = 0;

    else
        # peakjuv[r,temp_pos,dist_pos,sigtau_pos,tau_pos] = 1;
        # peakjuvquant[r,temp_pos,dist_pos,sigtau_pos,tau_pos] = abs(secondtoothj-maxtoothj);
        #NOTE: collect later
        # peakjuv[r,sigtau_pos,tau_pos] = 1;
        # peakjuvquant[r,sigtau_pos,tau_pos] = abs(secondtoothj-maxtoothj);

        peakjuvbin = 1;
        peakjuvdist = abs(secondtoothj-maxtoothj);
    end
    if secondpeaka == 0.0
        # peakadult[r,temp_pos,dist_pos,sigtau_pos,tau_pos] = 0;
        # peakadultquant[r,temp_pos,dist_pos,sigtau_pos,tau_pos] = 0;
        #NOTE: collect later
        # peakadult[r,sigtau_pos,tau_pos] = 0;
        # peakadultquant[r,sigtau_pos,tau_pos] = 0;

        peakadultbin = 0;
        peakadultdist = 0;
    else
        # peakadult[r,temp_pos,dist_pos,sigtau_pos,tau_pos] = 1;
        # peakadultquant[r,temp_pos,dist_pos,sigtau_pos,tau_pos] = abs(secondtootha-maxtootha);
        #NOTE: collect later
        # peakadult[r,sigtau_pos,tau_pos] = 1;
        # peakadultquant[r,sigtau_pos,tau_pos] = abs(secondtootha-maxtootha);
        
        peakadultbin = 1;
        peakadultdist = abs(secondtootha-maxtootha);
    end

    modesj = sort([secondtoothj,maxtoothj]);
    modesa = sort([secondtootha,maxtootha]);


    return mvj, mva, varj, vara, pj, pa, peakjuvbin, peakadultbin, peakjuvdist, peakadultdist, modesj, modesa, medianj, mediana, quartile25j, quartile25a, quartile75j, quartile75a





end
