if homedir() == "/home/z840"
    loadfunc = include("$(homedir())/sharks_bodysize/src/loadfuncs.jl");
else
    loadfunc = include("$(homedir())/Dropbox/PostDoc/2018_sharks/src/loadfuncs.jl");
end




#Size at birth (cm)
l0 = 100.0;
#Asymptotic size (cm)
L = 295; #1500; #295.0;

# # Mass from PRECAUDAL LENGTH (Schindler) 
# m0 = 0.00538776*l0^3.102025622731644;
# M = 0.00538776*L^3.102025622731644;

# Mass (kg) from TOTAL LENGTH (Goldman et al. 2006)
m0 = (0.00013*l0^2.4)*1000;
M = (0.00013*L^2.4)*1000;

#Sim params
n0=1000;
gen=1;


mintemp_j = [8,20,9,12];
maxtemp_j = [13,30,17,24];
mintemp_a = [8,20,12,9];
maxtemp_a = [13,30,24,17];

distvec = [200,400,650];
sigtauvec = collect(1:0.5:13);
tauvec = collect(1:1:25);
tauits = length(sigtauvec)*length(tauvec)*length(distvec);


#3 paramters: distance, sigtau, tau
distposvec = repeat(collect(1:3),inner = length(sigtauvec)*length(tauvec));
sigtauposvec = repeat(collect(1:length(sigtauvec)),inner = length(tauvec),outer=length(distvec));
tauposvec = repeat(collect(1:length(tauvec)),outer=length(distvec)*length(sigtauvec));
paramposvec_pre = [distposvec sigtauposvec tauposvec];
#temperature | distance | sigtauvec | tauposvec
paramposvec = [repeat(collect(1:4),inner = size(paramposvec_pre)[1]) repeat(paramposvec_pre,outer=4)];

reps = 25;
paramvec_pre = repeat(paramposvec,outer = reps);
paramvec = [repeat(collect(1:reps),inner=size(paramposvec)[1]) paramvec_pre];

its = size(paramvec)[1];


meanjuv = SharedArray{Float64}(reps,4,3,length(sigtauvec),length(tauvec));
meanadult = SharedArray{Float64}(reps,4,3,length(sigtauvec),length(tauvec));
varjuv = SharedArray{Float64}(reps,4,3,length(sigtauvec),length(tauvec));
varadult = SharedArray{Float64}(reps,4,3,length(sigtauvec),length(tauvec));
peakjuv = SharedArray{Int64}(reps,4,3,length(sigtauvec),length(tauvec));
peakadult = SharedArray{Int64}(reps,4,3,length(sigtauvec),length(tauvec));
peakadult2 = SharedArray{Int64}(reps,4,3,length(sigtauvec),length(tauvec));

peakjuvquant = SharedArray{Float64}(reps,4,3,length(sigtauvec),length(tauvec));
peakadultquant = SharedArray{Float64}(reps,4,3,length(sigtauvec),length(tauvec));

estmeanjuvG1 = SharedArray{Float64}(reps,4,3,length(sigtauvec),length(tauvec),2);
estvarjuvG1 = SharedArray{Float64}(reps,4,3,length(sigtauvec),length(tauvec),2);
estweightjuvG1 = SharedArray{Float64}(reps,4,3,length(sigtauvec),length(tauvec),2);
estmeanjuvG2 = SharedArray{Float64}(reps,4,3,length(sigtauvec),length(tauvec),2);
estvarjuvG2 = SharedArray{Float64}(reps,4,3,length(sigtauvec),length(tauvec),2);
estweightjuvG2 = SharedArray{Float64}(reps,4,3,length(sigtauvec),length(tauvec),2);

estmeanadultG1 = SharedArray{Float64}(reps,4,3,length(sigtauvec),length(tauvec),2);
estvaradultG1 = SharedArray{Float64}(reps,4,3,length(sigtauvec),length(tauvec),2);
estweightadultG1 = SharedArray{Float64}(reps,4,3,length(sigtauvec),length(tauvec),2);
estmeanadultG2 = SharedArray{Float64}(reps,4,3,length(sigtauvec),length(tauvec),2);
estvaradultG2 = SharedArray{Float64}(reps,4,3,length(sigtauvec),length(tauvec),2);
estweightadultG2 = SharedArray{Float64}(reps,4,3,length(sigtauvec),length(tauvec),2);

probunijuv = SharedArray{Float64}(reps,4,3,length(sigtauvec),length(tauvec));
probuniadult = SharedArray{Float64}(reps,4,3,length(sigtauvec),length(tauvec));


@time @sync @distributed for i=1:its
    # println(i)
    #set parameters
    pos = paramvec[i,:];
    r = pos[1];
    temp_pos = pos[2];
    dist_pos = pos[3];
    sigtau_pos = pos[4];
    tau_pos = pos[5];
    
    if temp_pos == 4 && dist_pos == 2
    
        #Temperature range (high latitude: 8-13; Low latitude 22-30)
        #Juvenile site
        tempmin1 = mintemp_j[temp_pos]+273.15; tempmax1 = maxtemp_j[temp_pos]+273.15;
        #Adult site
        tempmin2 = mintemp_a[temp_pos]+273.15; tempmax2 = maxtemp_j[temp_pos]+273.15;

        tempvec1 = Array{Float64}(undef,0);
        tempvec2 = Array{Float64}(undef,0);
        if tempmin1 == tempmax1
            tempvec1 = repeat([tempmin1],inner=100);
        else
            tempvec1 = collect(tempmin1:((tempmax1-tempmin1)/(100-1)):tempmax1);
        end;
        if tempmin2 == tempmax2
            tempvec2 = repeat([tempmin2],inner=100);
        else
            tempvec2 = collect(tempmin2:((tempmax2-tempmin2)/(100-1)):tempmax2);
        end;

        #distance between site (km * 1000)
        distance = distvec[dist_pos]*1000; #3.779e6/10; #1500
        #Shark velocity (m/s)
        velocity = 1;
        D = 1;

        #Steepness of juvenile migration (function of mass)
        sigtau = sigtauvec[sigtau_pos];
        #Steepness of adult migration (function of time)
        tau = tauvec[tau_pos];

        #save data
        filename = "data/sharks_3oceans2/simdata.jld";
        indices = [r,temp_pos,dist_pos,sigtau_pos,tau_pos];
        namespace = smartpath(filename,indices);
        
        @load namespace mass1 mass2 toothdrop toothlength1 toothlength2;
        
        
        #MEASURE SOMETHING
        mvj = dot(toothlength1[1,:],(toothdrop[:,1]/sum(toothdrop[:,1])));
        mva = dot(toothlength1[1,:],(toothdrop[:,2]/sum(toothdrop[:,2])));
        
        meanjuv[r,temp_pos,dist_pos,sigtau_pos,tau_pos] = mvj;
        meanadult[r,temp_pos,dist_pos,sigtau_pos,tau_pos] = mva;
        
        #calculate variance
        meanofxsquared = dot(toothlength1[1,:].^2,(toothdrop[:,1]/sum(toothdrop[:,1])));
        squareofmeanofx = mvj^2;
        varj = meanofxsquared - squareofmeanofx;
        varjuv[r,temp_pos,dist_pos,sigtau_pos,tau_pos] = varj;
        
        meanofxsquared = dot(toothlength1[1,:].^2,(toothdrop[:,2]/sum(toothdrop[:,2])));
        squareofmeanofx = mva^2;
        vara = meanofxsquared - squareofmeanofx;
        varadult[r,temp_pos,dist_pos,sigtau_pos,tau_pos] = vara;
        
        #Detect multiple modes (yes, no)
        #How many modes, and rank by significance
        #Where are the multiple modes (position of each? two most significant?)
        pj = (toothdrop[:,1]/sum(toothdrop[:,1]));
        pa = (toothdrop[:,2]/sum(toothdrop[:,2]));
        
        lmax=findlocalmaxima(pa); #
        lmin=findlocalminima(pa);
        #choose highest probability peaks
        peakprobs = pa[lmax];
        troughprobs = pa[lmin];
        sortlmax = lmax[sortperm(peakprobs)];
        sortlmin = lmin[sortperm(troughprobs)];
        sortedpeaks = pa[sortlmax];
        sortedtroughs = pa[sortlmin];
        maxpeaka = last(sortedpeaks);
        sortedtooth = toothlength1[1,:][lmax[sortperm(peakprobs)]];
        maxtootha = last(sortedtooth);
        secondpeaka = 0.0;
        secondtootha = 0.0;
        troughtootha = 0.0;
        troughmina = 0.0;
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

            trougha1 = pa[troughpos][findmin(pa[troughpos])[2]];
            if length(troughposbackup) > 0
                trougha2 = pa[troughposbackup][findmin(pa[troughposbackup])[2]];
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
                        secondpeaka = sortedpeaks[end-peakpos];
                        secondtootha = sortedtooth[end-peakpos];
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
            
            troughmina = [trougha1,trougha2][peakpos];
            troughtootha = toothlength1[1,:][troughposmin];
        else
            secondpeaka = 0.0;
            secondtootha = 0.0;
            troughtootha = 0.0;
            troughmina = 0.0;
        end
            
            
        
        # pl = lineplot(toothlength1[1,:],pa);
        # lineplot!(pl,repeat([maxtootha],outer=2),[0,maxpeaka],color=:red)
        # lineplot!(pl,repeat([secondtootha],outer=2),[0,secondpeaka],color=:red)
        # lineplot!(pl,repeat([troughtootha],outer=2),[0,troughmina],color=:blue)
        
        
        # 
        # filename = "density.pdf";
        # namespace = smartpath(filename)
        # R"""
        # pdf($namespace,width=6,height=10)
        # par(mfrow=c(2,1))
        # plot(density($foundteetha))
        # plot($(toothlength1[1,:]),$pa,type='b')
        # dev.off()
        # """
        # 
        
        
        lmax=findlocalmaxima(pj); #
        lmin=findlocalminima(pj);
        #choose highest probability peaks
        peakprobs = pj[lmax];
        troughprobs = pj[lmin];
        sortlmax = lmax[sortperm(peakprobs)];
        sortlmin = lmin[sortperm(troughprobs)];
        sortedpeaks = pj[sortlmax];
        sortedtroughs = pj[sortlmin];
        maxpeakj = last(sortedpeaks);
        sortedtooth = toothlength1[1,:][lmax[sortperm(peakprobs)]];
        maxtoothj = last(sortedtooth);
        secondpeakj = 0.0;
        secondtoothj = 0.0;
        troughtoothj = 0.0;
        troughminj = 0.0;
        troughpos = Array{Int64}(undef,0) .* 0;
        troughposbackup = Array{Int64}(undef,0) .* 0;
        troughj1 = Array{Float64}(undef,0) .* 0;
        troughj2 = Array{Float64}(undef,0) .* 0;
        troughj3 = Array{Float64}(undef,0) .* 0;
        troughj4 = Array{Float64}(undef,0) .* 0;
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

            troughj1 = pj[troughpos][findmin(pj[troughpos])[2]];
            if length(troughposbackup) > 0
                troughj2 = pj[troughposbackup][findmin(pj[troughposbackup])[2]];
            end
            
            #need to determine second peak significance
            if length(sortedpeaks) > 1
                #threshold of 25% size of first peak
                # if ((sortedpeaks[length(sortedpeaks)-1] / maxpeaka) > 0.25) && (maxpeaka - (trougha[1]/maxpeaka) > 0.25)
                peakprop[1] = (sortedpeaks[end-1] - troughj1) / (maxpeakj - troughj1);
                if length(troughposbackup) > 0
                    peakprop[2] = (sortedpeaks[end-2] - troughj2) / (maxpeakj - troughj2);
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
                        secondpeakj = sortedpeaks[end-peakpos];
                        secondtoothj = sortedtooth[end-peakpos];
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
            
            troughminj = [troughj1,troughj2][peakpos];
            troughtoothj = toothlength1[1,:][troughposmin];
        else
            secondpeakj = 0.0;
            secondtoothj = 0.0;
            troughtoothj = 0.0;
            troughminj = 0.0;
        end
        
        # 
        # 
        # 
        # pl = lineplot(toothlength1[1,:],pj);
        # lineplot!(pl,repeat([maxtoothj],outer=2),[0,maxpeakj],color=:red)
        # lineplot!(pl,repeat([secondtoothj],outer=2),[0,secondpeakj],color=:red)
        # lineplot!(pl,repeat([troughtoothj],outer=2),[0,troughminj],color=:blue)
        # 
            
        
        #Record presence of multiple modes
        if secondpeakj == 0.0
            peakjuv[r,temp_pos,dist_pos,sigtau_pos,tau_pos] = 0;
            peakjuvquant[r,temp_pos,dist_pos,sigtau_pos,tau_pos] = 0;
        else
            peakjuv[r,temp_pos,dist_pos,sigtau_pos,tau_pos] = 1;
            peakjuvquant[r,temp_pos,dist_pos,sigtau_pos,tau_pos] = abs(secondtoothj-maxtoothj);
        end
        if secondpeaka == 0.0
            peakadult[r,temp_pos,dist_pos,sigtau_pos,tau_pos] = 0;
            peakadultquant[r,temp_pos,dist_pos,sigtau_pos,tau_pos] = 0;
        else
            peakadult[r,temp_pos,dist_pos,sigtau_pos,tau_pos] = 1;
            peakadultquant[r,temp_pos,dist_pos,sigtau_pos,tau_pos] = abs(secondtootha-maxtootha);
        end
        
        #NOTE - TODO
        #1) run 1-mode or 2-mode depending on above ^^
        #2) save parameters
        #3) in a separate analysis, determine likelihood surface
        
        #Gaussian multi-mixture version
        #generate data
        numteethj = Int64.(floor.(pj*100000));
        foundteethj = Array{Float64}(undef,sum(numteethj));
        let tic = 1
            for i=1:length(numteethj)
                newteethj = repeat([toothlength1[1,i]],outer=numteethj[i]) .+ rand(Normal(0,0.01),numteethj[i]);
                foundteethj[tic:((tic-1) + numteethj[i])] = newteethj;
                tic += numteethj[i];
            end
        end
        
        numteetha = Int64.(floor.(pa*100000));
        foundteetha = Array{Float64}(undef,sum(numteetha));
        let tic = 1
            for i=1:length(numteetha)
                newteetha = repeat([toothlength1[1,i]],outer=numteetha[i]) .+ rand(Normal(0,0.01),numteetha[i]);
                foundteetha[tic:((tic-1) + numteetha[i])] = newteetha;
                tic += numteetha[i];
            end
        end
        #We are going to COMPARE likelihoods
        #Juvenile parameterization
        # if peakjuv[r,temp_pos,dist_pos,sigtau_pos,tau_pos] == 0
            R"""
            library(mclust)
            x.gmmj1 <- Mclust($foundteethj,G=1,verbose=FALSE)
            # x.gmm.2 = Mclust($foundteetha,G=2,verbose=FALSE)
            estmeanjG1 <- as.numeric(x.gmmj1$parameters$mean)
            estvarjG1 <- as.numeric(x.gmmj1$parameters$mean)
            estweightjG1 <-  as.numeric(x.gmmj1$parameters$pro)
            """
        # else
            R"""
            library(mclust)
            # x.gmm.1 = Mclust($foundteetha,G=1,verbose=FALSE)
            x.gmmj2 <- Mclust($foundteethj,G=2,verbose=FALSE)
            estmeanjG2 <- as.numeric(x.gmmj2$parameters$mean)
            estvarjG2 <- as.numeric(x.gmmj2$parameters$variance$sigmasq)
            estweightjG2 <- as.numeric(x.gmmj2$parameters$pro)
            pdj = 1-pchisq(logLik(x.gmmj2)[1] - logLik(x.gmmj1)[1], df=x.gmmj2$df - x.gmmj1$df);
            """
        # end
        #Adult parameterization
        # if peakadult[r,temp_pos,dist_pos,sigtau_pos,tau_pos] == 0
            R"""
            library(mclust)
            x.gmma1 <- Mclust($foundteetha,G=1,verbose=FALSE)
            # x.gmm.2 = Mclust($foundteetha,G=2,verbose=FALSE)
            estmeanaG1 <- as.numeric(x.gmma1$parameters$mean)
            estvaraG1 <- as.numeric(x.gmma1$parameters$mean)
            estweightaG1 <-  as.numeric(x.gmma1$parameters$pro)
            """
        # else
            R"""
            library(mclust)
            # x.gmm.1 = Mclust($foundteetha,G=1,verbose=FALSE)
            x.gmma2 <- Mclust($foundteetha,G=2,verbose=FALSE)
            estmeanaG2 <- as.numeric(x.gmma2$parameters$mean)
            estvaraG2 <- as.numeric(x.gmma2$parameters$variance$sigmasq)
            estweightaG2 <- as.numeric(x.gmma2$parameters$pro)
            pda = 1-pchisq(logLik(x.gmma2)[1] - logLik(x.gmma1)[1], df=x.gmma2$df - x.gmma1$df);
            """
        # end
        #likelihood ratio test for 1 vs. 2 components
        # likediff = logLik(x.gmm.2)-logLik(x.gmm.1)
        # pvalue = 1-pchisq(likediff, df=3)
        #If pvalue is < 0.05, accept 2 modes
        #If pvalue is > 0.05, accept 1 mode
        
        #or likelihood ratio test?
        # likeratiotest = exp(logLik(x.gmm.1)-logLik(x.gmm.2))
        #if < c then reject the null gmm.1
        #if > c then accept the null gmm.1
        # summary(x.gmm)
        
        @rget estmeanjG1;
        @rget estvarjG1;
        @rget estweightjG1;
        @rget estmeanjG2;
        @rget estvarjG2;
        @rget estweightjG2;
        @rget pdj
        estmeanjuvG1[r,temp_pos,dist_pos,sigtau_pos,tau_pos,:] .= estmeanjG1;
        estvarjuvG1[r,temp_pos,dist_pos,sigtau_pos,tau_pos,:] .= estvarjG1;
        estweightjuvG1[r,temp_pos,dist_pos,sigtau_pos,tau_pos,:] .= estweightjG1;
        estmeanjuvG2[r,temp_pos,dist_pos,sigtau_pos,tau_pos,:] .= estmeanjG2;
        estvarjuvG2[r,temp_pos,dist_pos,sigtau_pos,tau_pos,:] .= estvarjG2;
        estweightjuvG2[r,temp_pos,dist_pos,sigtau_pos,tau_pos,:] .= estweightjG2;
        probunijuv[r,temp_pos,dist_pos,sigtau_pos,tau_pos] = pdj;
        
        @rget estmeanaG1;
        @rget estvaraG1;
        @rget estweightaG1;
        @rget estmeanaG2;
        @rget estvaraG2;
        @rget estweightaG2;
        @rget pda;
        estmeanadultG1[r,temp_pos,dist_pos,sigtau_pos,tau_pos,:] .= estmeanaG1;
        estvaradultG1[r,temp_pos,dist_pos,sigtau_pos,tau_pos,:] .= estvaraG1;
        estweightadultG1[r,temp_pos,dist_pos,sigtau_pos,tau_pos,:] .= estweightaG1;
        estmeanadultG2[r,temp_pos,dist_pos,sigtau_pos,tau_pos,:] .= estmeanaG2;
        estvaradultG2[r,temp_pos,dist_pos,sigtau_pos,tau_pos,:] .= estvaraG2;
        estweightadultG2[r,temp_pos,dist_pos,sigtau_pos,tau_pos,:] .= estweightaG2;
        probuniadult[r,temp_pos,dist_pos,sigtau_pos,tau_pos] = pda;
        
        # @rget likeratiotest
        
        # #accept null hypothesis of one mode
        # if likeratiotest > 0.05
        #     peakadult2[r,temp_pos,dist_pos,sigtau_pos,tau_pos] = 0;
        # #else reject null hypothesis of one mode
        # else
        #     peakadult2[r,temp_pos,dist_pos,sigtau_pos,tau_pos] = 1;
        # end
        # #JUVENILE PEAKS
        
        
        
    end #end if loop
    
    
    
end #end its

filename = "data/sharks_3oceans2/parameterest_data.jld";
namespace = smartpath(filename);
@save namespace meanjuv meanadult peakjuvquant peakadultquant estmeanjuvG1 estvarjuvG1 estweightjuvG1 estmeanjuvG2 estvarjuvG2 estweightjuvG2 estmeanadultG1 estvaradultG1 estweightadultG1 estmeanadultG2 estvaradultG2 estweightadultG2 probunijuv probuniadult;




#Import Site Data
data = 
#Construct likelihood surface
#Right now, we will construct a likelihood for *both* adult and juv distributions... agnostic to data source.
#Not sure (doubt) if likelihoods can be directly compared since they are from 2 different likelihood functions? Kind of but not, because the unimodal likelihood is the same as the bimodal likelihood except one weight is zero. So maybe the same and directly comparable!

#NOTE: TESING A SINGLE DATASET
#Does the max likelihood point back to the correct parameters???
#save data
r=1
temp_pos = 4;
dist_pos = 2;
sigtau_pos = 15; 
tau_pos = 20;
println([sigtauvec[sigtau_pos],tauvec[tau_pos]])
filename = "data/sharks_3oceans2/simdata.jld";
indices = [r,temp_pos,dist_pos,sigtau_pos,tau_pos];
namespace = smartpath(filename,indices);
@load namespace mass1 mass2 toothdrop toothlength1 toothlength2;
pj = (toothdrop[:,1]/sum(toothdrop[:,1]));
pa = (toothdrop[:,2]/sum(toothdrop[:,2]));
numteethj = Int64.(floor.(pj*1000));
foundteethj = Array{Float64}(undef,sum(numteethj));
let tic = 1
    for i=1:length(numteethj)
        newteethj = repeat([toothlength1[1,i]],outer=numteethj[i]) .+ rand(Normal(0,0.01),numteethj[i]);
        foundteethj[tic:((tic-1) + numteethj[i])] = newteethj;
        tic += numteethj[i];
    end
end
numteetha = Int64.(floor.(pa*1000));
foundteetha = Array{Float64}(undef,sum(numteetha));
let tic = 1
    for i=1:length(numteetha)
        newteetha = repeat([toothlength1[1,i]],outer=numteetha[i]) .+ rand(Normal(0,0.01),numteetha[i]);
        foundteetha[tic:((tic-1) + numteetha[i])] = newteetha;
        tic += numteetha[i];
    end
end


mnegllG1_juv = Array{Float64}(undef,length(sigtauvec),length(tauvec));
mnegllG2_juv = Array{Float64}(undef,length(sigtauvec),length(tauvec));
mnegllG1_adult = Array{Float64}(undef,length(sigtauvec),length(tauvec));
mnegllG2_adult = Array{Float64}(undef,length(sigtauvec),length(tauvec));

mR_juv = Array{Float64}(undef,length(sigtauvec),length(tauvec));
mR_adult = Array{Float64}(undef,length(sigtauvec),length(tauvec));

compllG2 = Array{Float64}(undef,length(sigtauvec),length(tauvec));

temp_pos = 4;
dist_pos = 2;
for i = 1:length(sigtauvec)
    sigtau_pos = copy(i);
    for j = 1:length(tauvec)
        tau_pos = copy(j);
        
        negllG1_juv = Array{Float64}(undef,reps);
        negllG2_juv = Array{Float64}(undef,reps);
        negllG1_adult = Array{Float64}(undef,reps);
        negllG2_adult = Array{Float64}(undef,reps);
        for r=1:reps
            
            
            # if jweight[1] > 1
                #Unimodal
                data = copy(foundteethj);
                #JUVENILE LIKELIHOOD
                jmean = estmeanjuvG1[r,temp_pos,dist_pos,sigtau_pos,tau_pos,1];
                jvar = estvarjuvG1[r,temp_pos,dist_pos,sigtau_pos,tau_pos,1];
                jweight = estweightjuvG1[r,temp_pos,dist_pos,sigtau_pos,tau_pos,1];
                
                mu1 = jmean[1];
                sigma1 = sqrt(jvar[1]);
                # for x=1:length(data)
                #     x1 = data[x];
                #     logllvec_juv[x] = log(exp(-(x1-mu1)^2/(2*sigma1^2))/(sqrt(2*pi)*sigma1));
                # end
                negllG1_juv[r] = sum(((data .- mu1).^2) ./ (2*sigma1^2)) + length(data)*(log(sigma1) + 0.5*log(2*pi));
                # negll_juv[r] = -sum(logllvec_juv,dims=2);
            # else
                jmean = estmeanjuvG2[r,temp_pos,dist_pos,sigtau_pos,tau_pos,:];
                jvar = estvarjuvG2[r,temp_pos,dist_pos,sigtau_pos,tau_pos,:];
                jweight = estweightjuvG2[r,temp_pos,dist_pos,sigtau_pos,tau_pos,:];
                
                
                #Bimodal
                #JUVENILE LIKELIHOOD
                mu1 = jmean[1];
                mu2 = jmean[2];
                sigma1 = sqrt(jvar[1]);
                sigma2 = sqrt(jvar[2]);
                w1 = jweight[1];
                w2 = jweight[2];
                logllvec_juv = Array{Float64}(undef,length(data));
                for x=1:length(data)
                    x1 = data[x];
                    logllvec_juv[x] = log((exp(-(x1-mu1)^2/(2*sigma1^2))*w1)/(sqrt(2*pi)*sigma1*(w1 + w2)) + (exp(-(x1-mu2)^2/(2*sigma2^2))*w2)/(sqrt(2*pi)*sigma2*(w1 + w2)));
                end
                negllG2_juv[r] = -sum(logllvec_juv);
            # end
            # if aweight[1] > 1
            data = copy(foundteetha);
                #Unimodal
                #ADULT LIKELIHOOD
                amean = estmeanadultG1[r,temp_pos,dist_pos,sigtau_pos,tau_pos,1];
                avar = estvaradultG1[r,temp_pos,dist_pos,sigtau_pos,tau_pos,1];
                aweight = estweightadultG1[r,temp_pos,dist_pos,sigtau_pos,tau_pos,1];
                mu1 = amean[1];
                sigma1 = sqrt(avar[1]);
                negllG1_adult[r] = sum(((data .- mu1).^2) ./ (2*sigma1^2)) + length(data)*(log(sigma1) + 0.5*log(2*pi));
                # for x=1:length(data)
                #     x1 = data[x];
                #     logllvec_adult[x] = log(exp(-(x1-mu1)^2/(2*sigma1^2))/(sqrt(2*pi)*sigma1));
                # end
                # negll_adult[r] = -sum(logllvec_adult);
            # else
                #Bimodal
                #Adult LIKELIHOOD
                amean = estmeanadultG2[r,temp_pos,dist_pos,sigtau_pos,tau_pos,:];
                avar = estvaradultG2[r,temp_pos,dist_pos,sigtau_pos,tau_pos,:];
                aweight = estweightadultG2[r,temp_pos,dist_pos,sigtau_pos,tau_pos,:];
                mu1 = amean[1];
                mu2 = amean[2];
                sigma1 = sqrt(avar[1]);
                sigma2 = sqrt(avar[2]);
                w1 = aweight[1];
                w2 = aweight[2];
                logllvec_adult = Array{Float64}(undef,length(data));
                for x=1:length(data)
                    x1 = data[x];
                    logllvec_adult[x] = log((exp(-(x1-mu1)^2/(2*sigma1^2))*w1)/(sqrt(2*pi)*sigma1*(w1 + w2)) + (exp(-(x1-mu2)^2/(2*sigma2^2))*w2)/(sqrt(2*pi)*sigma2*(w1 + w2)));
                end
                negllG2_adult[r] = -sum(logllvec_adult);
            # end
        end
        #sum to get likelihood, take mean across reps
        mnegllG1_juv[i,j] = mean(negllG1_juv);
        mnegllG2_juv[i,j] = mean(negllG2_juv);
        mnegllG1_adult[i,j] = mean(negllG1_adult);
        mnegllG2_adult[i,j] = mean(negllG2_adult);
        
        compllG2[i,j] = mean(negllG2_juv .+ negllG2_adult);
        
        
        mR_juv[i,j] = mean(2*(negllG1_juv .- negllG2_juv));
        mR_adult[i,j] = mean(2*(negllG1_adult .- negllG2_adult));
    end
end

findmin(mnegllG1_juv)
findmin(mnegllG2_juv)

findmin(mnegllG1_adult)
findmin(mnegllG2_adult)



minl = minimum([mnegllG2_juv;mnegll_adult])
maxl = maximum([mnegllG2_juv;mnegll_adult])
filename = "figures/likelihood_test.pdf";
namespace = smartpath(filename);
R"""
library(RColorBrewer)
library(fields)
pal = colorRampPalette(rev(brewer.pal(11,'Spectral')))(50)
pdf($namespace,width=16,height=6)
par(mfrow=c(1,2))
image.plot(x = $sigtauvec,y = $tauvec, z = $mnegll_juv,zlim=c($minl,$maxl),xlab='Juvenile migration window',ylab='Adult migration window',col=pal)
image.plot(x = $sigtauvec,y = $tauvec, z = $mnegll_adult,zlim=c($minl,$maxl),xlab='Juvenile migration window',ylab='Adult migration window',col=pal)
dev.off()
"""

#Build a composite likelihood by adding the log-likelihoods
# comp_llG2 = mnegllG2_juv .+ mnegllG2_adult;
# comp_llG1 = mnegllG1_juv .+ mnegllG1_adult;
minpoint = findmin(compllG2)[2]
filename = "figures/compositelikelihood.pdf";
namespace = smartpath(filename);
R"""
library(RColorBrewer)
library(fields)
pal = colorRampPalette(rev(brewer.pal(11,'Spectral')))(50)
pdf($namespace,width=8,height=6)
par(mfrow=c(1,1))
image.plot(x = $sigtauvec,y = $tauvec, z = $compllG2,xlab='Juvenile migration window',ylab='Adult migration window',col=pal)
points($(sigtauvec[minpoint[1]]),$(tauvec[minpoint[2]]),pch=16,col='white')
points($(sigtauvec[sigtau_pos]),$(tauvec[tau_pos]),pch=1,col='black')
# image.plot(x = $sigtauvec,y = $tauvec, z = $mnegll_adult,zlim=c($minl,$maxl),xlab='Juvenile migration window',ylab='Adult migration window',col=pal)
dev.off()
"""



#Generating a dataset for every sigtau/tau value
#Generating likelihoods for G1 and G2
#Calculating likelihood ratio

mnegllG1_juv = Array{Float64}(undef,length(sigtauvec),length(tauvec));
mnegllG2_juv = Array{Float64}(undef,length(sigtauvec),length(tauvec));
mnegllG1_adult = Array{Float64}(undef,length(sigtauvec),length(tauvec));
mnegllG2_adult = Array{Float64}(undef,length(sigtauvec),length(tauvec));

mR_juv = Array{Float64}(undef,length(sigtauvec),length(tauvec));
mR_adult = Array{Float64}(undef,length(sigtauvec),length(tauvec));

# r=1
temp_pos = 4;
dist_pos = 2;
for i = 1:length(sigtauvec)
    for j = 1: length(tauvec)
        sigtau_pos = copy(i); 
        tau_pos = copy(j);

# temp_pos = 4;
# dist_pos = 2;
# for i = 1:length(sigtauvec)
#     sigtau_pos = copy(i);
#     for j = 1:length(tauvec)
#         tau_pos = copy(j);
        
        negllG1_juv = Array{Float64}(undef,reps);
        negllG2_juv = Array{Float64}(undef,reps);
        negllG1_adult = Array{Float64}(undef,reps);
        negllG2_adult = Array{Float64}(undef,reps);
        for r=1:reps
            # println([sigtauvec[sigtau_pos],tauvec[tau_pos]])
            filename = "data/sharks_3oceans2/simdata.jld";
            indices = [r,temp_pos,dist_pos,sigtau_pos,tau_pos];
            namespace = smartpath(filename,indices);
            @load namespace mass1 mass2 toothdrop toothlength1 toothlength2;
            pj = (toothdrop[:,1]/sum(toothdrop[:,1]));
            pa = (toothdrop[:,2]/sum(toothdrop[:,2]));
            numteethj = Int64.(floor.(pj*15));
            foundteethj = Array{Float64}(undef,sum(numteethj));
            let tic = 1
                for i=1:length(numteethj)
                    newteethj = repeat([toothlength1[1,i]],outer=numteethj[i]) .+ rand(Normal(0,0.2),numteethj[i]);
                    foundteethj[tic:((tic-1) + numteethj[i])] = newteethj;
                    tic += numteethj[i];
                end
            end
            numteetha = Int64.(floor.(pa*15));
            foundteetha = Array{Float64}(undef,sum(numteetha));
            let tic = 1
                for i=1:length(numteetha)
                    newteetha = repeat([toothlength1[1,i]],outer=numteetha[i]) .+ rand(Normal(0,0.2),numteetha[i]);
                    foundteetha[tic:((tic-1) + numteetha[i])] = newteetha;
                    tic += numteetha[i];
                end
            end
            
            
            data = copy(foundteethj);

            
            # if jweight[1] > 1
                #Unimodal
                #JUVENILE LIKELIHOOD
                jmean = estmeanjuvG1[r,temp_pos,dist_pos,sigtau_pos,tau_pos,1];
                jvar = estvarjuvG1[r,temp_pos,dist_pos,sigtau_pos,tau_pos,1];
                jweight = estweightjuvG1[r,temp_pos,dist_pos,sigtau_pos,tau_pos,1];
                
                mu1 = jmean[1];
                sigma1 = sqrt(jvar[1]);
                # for x=1:length(data)
                #     x1 = data[x];
                #     logllvec_juv[x] = log(exp(-(x1-mu1)^2/(2*sigma1^2))/(sqrt(2*pi)*sigma1));
                # end
                negllG1_juv[r] = sum(((data .- mu1).^2) ./ (2*sigma1^2)) + length(data)*(log(sigma1) + 0.5*log(2*pi));
                # negll_juv[r] = -sum(logllvec_juv,dims=2);
            # else
                jmean = estmeanjuvG2[r,temp_pos,dist_pos,sigtau_pos,tau_pos,:];
                jvar = estvarjuvG2[r,temp_pos,dist_pos,sigtau_pos,tau_pos,:];
                jweight = estweightjuvG2[r,temp_pos,dist_pos,sigtau_pos,tau_pos,:];
                
                #Bimodal
                #JUVENILE LIKELIHOOD
                mu1 = jmean[1];
                mu2 = jmean[2];
                sigma1 = sqrt(jvar[1]);
                sigma2 = sqrt(jvar[2]);
                w1 = jweight[1];
                w2 = jweight[2];
                logllvec_juv = Array{Float64}(undef,length(data));
                for x=1:length(data)
                    x1 = data[x];
                    logllvec_juv[x] = log((exp(-(x1-mu1)^2/(2*sigma1^2))*w1)/(sqrt(2*pi)*sigma1*(w1 + w2)) + (exp(-(x1-mu2)^2/(2*sigma2^2))*w2)/(sqrt(2*pi)*sigma2*(w1 + w2)));
                end
                negllG2_juv[r] = -sum(logllvec_juv);
            # end
            # if aweight[1] > 1
            
            data = copy(foundteetha);
            
                #Unimodal
                #ADULT LIKELIHOOD
                amean = estmeanadultG1[r,temp_pos,dist_pos,sigtau_pos,tau_pos,1];
                avar = estvaradultG1[r,temp_pos,dist_pos,sigtau_pos,tau_pos,1];
                aweight = estweightadultG1[r,temp_pos,dist_pos,sigtau_pos,tau_pos,1];
                mu1 = amean[1];
                sigma1 = sqrt(avar[1]);
                negllG1_adult[r] = sum(((data .- mu1).^2) ./ (2*sigma1^2)) + length(data)*(log(sigma1) + 0.5*log(2*pi));
                # for x=1:length(data)
                #     x1 = data[x];
                #     logllvec_adult[x] = log(exp(-(x1-mu1)^2/(2*sigma1^2))/(sqrt(2*pi)*sigma1));
                # end
                # negll_adult[r] = -sum(logllvec_adult);
            # else
                #Bimodal
                #Adult LIKELIHOOD
                amean = estmeanadultG2[r,temp_pos,dist_pos,sigtau_pos,tau_pos,:];
                avar = estvaradultG2[r,temp_pos,dist_pos,sigtau_pos,tau_pos,:];
                aweight = estweightadultG2[r,temp_pos,dist_pos,sigtau_pos,tau_pos,:];
                mu1 = amean[1];
                mu2 = amean[2];
                sigma1 = sqrt(avar[1]);
                sigma2 = sqrt(avar[2]);
                w1 = aweight[1];
                w2 = aweight[2];
                logllvec_adult = Array{Float64}(undef,length(data));
                for x=1:length(data)
                    x1 = data[x];
                    logllvec_adult[x] = log((exp(-(x1-mu1)^2/(2*sigma1^2))*w1)/(sqrt(2*pi)*sigma1*(w1 + w2)) + (exp(-(x1-mu2)^2/(2*sigma2^2))*w2)/(sqrt(2*pi)*sigma2*(w1 + w2)));
                end
                negllG2_adult[r] = -sum(logllvec_adult);
            # end
        end
        #sum to get likelihood, take mean across reps
        mnegllG1_juv[i,j] = mean(negllG1_juv);
        mnegllG2_juv[i,j] = mean(negllG2_juv);
        mnegllG1_adult[i,j] = mean(negllG1_adult);
        mnegllG2_adult[i,j] = mean(negllG2_adult);
        
        
        mR_juv[i,j] = mean(2*(negllG1_juv .- negllG2_juv));
        mR_adult[i,j] = mean(2*(negllG1_adult .- negllG2_adult));
    end
end



findmin(mnegllG1_juv)
findmin(mnegllG2_juv)

findmin(mnegllG1_adult)
findmin(mnegllG2_adult)




# minl = minimum([mnegll_juv;mnegll_adult])
# maxl = maximum([mnegll_juv;mnegll_adult])
filename = "figures/likelihoodratio.pdf";
namespace = smartpath(filename);
R"""
library(RColorBrewer)
library(fields)
pal = colorRampPalette(rev(brewer.pal(11,'Spectral')))(50)
pdf($namespace,width=16,height=6)
par(mfrow=c(1,2))
image.plot(x = $sigtauvec,y = $tauvec, z = $mR_juv,zlim=c(0,max($mR_adult)),xlab='Juvenile migration window',ylab='Adult migration window',col=pal)
image.plot(x = $sigtauvec,y = $tauvec, z = $mR_adult,zlim=c(0,max($mR_adult)),xlab='Juvenile migration window',ylab='Adult migration window',col=pal)
dev.off()
"""

filename = "figures/probunimodal.pdf";
namespace = smartpath(filename);
R"""
library(RColorBrewer)
library(fields)
pal = colorRampPalette(rev(brewer.pal(11,'Spectral')))(50)
pdf($namespace,width=16,height=6)
par(mfrow=c(1,2))
image.plot(x = $sigtauvec,y = $tauvec, z = $(mean(probunijuv[:,4,2,:,:],dims=1)[1,:,:]),zlim=c(0,1),xlab='Juvenile migration window',ylab='Adult migration window',col=pal)
image.plot(x = $sigtauvec,y = $tauvec, z = $(mean(probuniadult[:,4,2,:,:],dims=1)[1,:,:]),zlim=c(0,1),xlab='Juvenile migration window',ylab='Adult migration window',col=pal)
dev.off()
"""


#IF we apply MClustG=2 across EVERYTHING, plot difference in means (modes) of mixed normal
temp_pos = 4;
dist_pos = 2;
diffmod_juv = Array{Float64}(undef,length(sigtauvec),length(tauvec));
diffmod_adult = Array{Float64}(undef,length(sigtauvec),length(tauvec));
for i = 1:length(sigtauvec)
    sigtau_pos = copy(i);
    for j=1:length(tauvec)
        tau_pos = copy(j);
        jdiff = Array{Float64}(undef,reps);
        adiff = Array{Float64}(undef,reps);
        for r = 1:reps
            juvmodes = estmeanjuv[r,temp_pos,dist_pos,sigtau_pos,tau_pos,:];
            adultmodes = estmeanadult[r,temp_pos,dist_pos,sigtau_pos,tau_pos,:];
            jdiff[r] = sqrt((juvmodes[1]-juvmodes[2])^2);
            adiff[r] = sqrt((adultmodes[1]-adultmodes[2])^2);
        end
        diffmod_juv[i,j] = mean(jdiff);
        diffmod_adult[i,j] = mean(adiff);
    end
end


# SHOWCASE A SINGLE TREATMENT!
filename = "figures/fig_means_peaks3.pdf";
namespace = smartpath(filename);
i = 4; #temp regime
j = 2; #dist regime
mmeanjuv = mean(meanjuv,dims=1)[1,:,:,:,:];
mmeanadult = mean(meanadult,dims=1)[1,:,:,:,:];
minsize = minimum([mmeanjuv[i,j,:,:]; mmeanadult[i,j,:,:]]);
maxsize = maximum([mmeanjuv[i,j,:,:]; mmeanadult[i,j,:,:]]);
R"""
library(RColorBrewer)
library(fields)
pal = colorRampPalette(rev(brewer.pal(11,'Spectral')))(50)
pdf($namespace,width=10,height=8)
layout(matrix(c(1,2,3,4),2,2,byrow=TRUE))
par(oma = c(2, 2, 2, 1), mar = c(4, 5, 2, 5)) #,mai=c(0.6,0.6,0,0.1)
image.plot(x=$sigtauvec,y=$tauvec,z=$(mmeanjuv[i,j,:,:]),xlab='Juvenile migration window',ylab='Adult migration window',col=pal,zlim=c($minsize,$maxsize))
image.plot(x=$sigtauvec,y=$tauvec,z=$(mmeanadult[i,j,:,:]),xlab='Juvenile migration window',ylab='Adult migration window',col=pal,zlim=c($minsize,$maxsize))
"""
# minmoddiff = minimum(([diffmod_juv[3:25,3:25];diffmod_adult[3:25,3:25]]));
minmoddiff = 3;
maxmoddiff = maximum(([diffmod_juv[2:25,2:25];diffmod_adult[2:25,2:25]]));
R"""
pal = c('white',colorRampPalette((brewer.pal(11,'YlGnBu')))(50))
image.plot(x=$sigtauvec,y=$tauvec,z=$((diffmod_juv)),xlab='Juvenile migration window',ylab='Adult migration window',col=pal,zlim=c($minmoddiff,$maxmoddiff))
image.plot(x=$sigtauvec,y=$tauvec,z=$((diffmod_adult)),xlab='Juvenile migration window',ylab='Adult migration window',col=pal,zlim=c($minmoddiff,$maxmoddiff))
dev.off()
"""



#Build a COMPOSITE LIKELIHOOD by ADDING THE LOG-LIKELIHOODS






filename = "figures/fig_3oceans_means2.pdf";
namespace = smartpath(filename);
toothmatrix_juv = (meanjuv[1,1,:,:]);
toothmatrix_adult = meanadult[1,1,:,:];
minsize = minimum([vec(toothmatrix_juv) vec(toothmatrix_adult)]);
maxsize = maximum([vec(toothmatrix_juv) vec(toothmatrix_adult)]);
R"""
library(RColorBrewer)
library(fields)
pal = colorRampPalette(rev(brewer.pal(11,'Spectral')))(50)
pdf($namespace,width=14,height=6)
layout(matrix(c(1,2),1,2,byrow=TRUE))
par(oma = c(2, 2, 2, 1), mar = c(4, 5, 2, 5)) #,mai=c(0.6,0.6,0,0.1)
image.plot(x=$sigtauvec,y=$tauvec,z=$(toothmatrix_juv),xlab='Juvenile migration window',ylab='Adult migration window',main='Juvenile site tooth means',col=pal,zlim=c($minsize,$maxsize))
image.plot(x=$sigtauvec,y=$tauvec,z=$(toothmatrix_adult),xlab='Juvenile migration window',ylab='Adult migration window',main='Adult site tooth means',col=pal,zlim=c($minsize,$maxsize))
dev.off()
"""
filename = "figures/fig_3oceans_SD.pdf";
namespace = smartpath(filename);
toothmatrix_juv = sqrt.(varjuv[1,1,:,:]);
toothmatrix_adult = sqrt.(varadult[1,1,:,:]);
minsd = minimum([vec(toothmatrix_juv) vec(toothmatrix_adult)]);
maxsd = maximum([vec(toothmatrix_juv) vec(toothmatrix_adult)]);
R"""
library(RColorBrewer)
library(fields)
pal = colorRampPalette(rev(brewer.pal(11,'Spectral')))(50)
pdf($namespace,width=14,height=6)
layout(matrix(c(1,2),1,2,byrow=TRUE))
par(oma = c(2, 2, 2, 1), mar = c(4, 5, 2, 5)) #,mai=c(0.6,0.6,0,0.1)
image.plot(x=$sigtauvec,y=$tauvec,z=$(toothmatrix_juv),xlab='Juvenile migration window',ylab='Adult migration window',main='Juvenile site tooth SD',col=pal,zlim=c($minsd,$maxsd))
image.plot(x=$sigtauvec,y=$tauvec,z=$(toothmatrix_adult),xlab='Juvenile migration window',ylab='Adult migration window',main='Adult site tooth SD',col=pal,zlim=c($minsd,$maxsd))
dev.off()
"""

filename = "figures/fig_3oceans_peaks2.pdf";
namespace = smartpath(filename);
R"""
library(RColorBrewer)
library(fields)
pal = colorRampPalette(rev(brewer.pal(11,'Spectral')))(50)
pdf($namespace,width=14,height=6)
layout(matrix(c(1,2),1,2,byrow=TRUE))
par(oma = c(2, 2, 2, 1), mar = c(4, 5, 2, 5)) #,mai=c(0.6,0.6,0,0.1)
image.plot(x=$sigtauvec,y=$tauvec,z=$(peakjuv[3,1,:,:]),xlab='Juvenile migration window',ylab='Adult migration window',main='Juvenile site tooth peaks',col=c('white','black'))
image.plot(x=$sigtauvec,y=$tauvec,z=$(peakadult[3,1,:,:]),xlab='Juvenile migration window',ylab='Adult migration window',main='Adult site tooth peaks',col=c('white','black'))
dev.off()
"""


#Take means across reps
mpeakjuv = mean(peakjuv,dims=1)[1,:,:,:,:];
mpeakadult = mean(peakadult,dims=1)[1,:,:,:,:];

filename = "figures/fig_mpeaksjuv2.pdf";
namespace = smartpath(filename);
R"""
library(RColorBrewer)
library(fields)
pal = colorRampPalette((brewer.pal(11,'Greys')))(10)
pdf($namespace,width=12,height=10)
layout(matrix(seq(1,12), 4, 3, byrow = TRUE), 
   widths=rep(1,9), heights=rep(1,9))
par(oma = c(4, 5, 1, 1), mar = c(1, 0, 0, 1))
"""
for i=1:4
    for j=1:3
        R"""
        image(x=$sigtauvec,y=$tauvec,z=$(mpeakjuv[i,j,:,:]),xlab='Juvenile migration window',ylab='Adult migration window',main='Juvenile site tooth peaks',col=pal,zlim=c(0,1))
        """
    end
end
R"""
dev.off()
"""

filename = "figures/fig_mpeaksadult2.pdf";
namespace = smartpath(filename);
R"""
library(RColorBrewer)
library(fields)
pal = colorRampPalette((brewer.pal(11,'Greys')))(10)
pdf($namespace,width=12,height=10)
layout(matrix(seq(1,12), 4, 3, byrow = TRUE), 
   widths=rep(1,9), heights=rep(1,9))
par(oma = c(4, 5, 1, 1), mar = c(1, 0, 0, 1))
"""
for i=1:4
    for j=1:3
        R"""
        image(x=$sigtauvec,y=$tauvec,z=$(mpeakadult[i,j,:,:]),xlab='Juvenile migration window',ylab='Adult migration window',main='Adult site tooth peaks',col=pal,zlim=c(0,1))
        """
    end
end
R"""
dev.off()
"""

#QUANTITATIVE PEAKS

#Take means across reps
mpeakjuvquant = mean(peakjuvquant,dims=1)[1,:,:,:,:];
mpeakadultquant = mean(peakadultquant,dims=1)[1,:,:,:,:];
maxpeakquant = maximum([mpeakjuvquant;mpeakadultquant]);

filename = "figures/fig_mpeaksjuvquant2.pdf";
namespace = smartpath(filename);
R"""
library(RColorBrewer)
library(fields)
pal = c('white',colorRampPalette((brewer.pal(11,'YlGnBu')))(20))
pdf($namespace,width=12,height=10)
layout(matrix(seq(1,12), 4, 3, byrow = TRUE), 
   widths=rep(1,9), heights=rep(1,9))
par(oma = c(4, 5, 1, 1), mar = c(1, 0, 0, 1))
"""
for i=1:4
    for j=1:3
        R"""
        image(x=$sigtauvec,y=$tauvec,z=$(mpeakjuvquant[i,j,:,:]),xlab='Juvenile migration window',ylab='Adult migration window',main='Juvenile site tooth peaks',col=pal,zlim=c(0,$maxpeakquant))
        """
    end
end
R"""
dev.off()
"""

filename = "figures/fig_mpeaksadultquant2.pdf";
namespace = smartpath(filename);
R"""
library(RColorBrewer)
library(fields)
pal = c('white',colorRampPalette((brewer.pal(11,'YlGnBu')))(20))
pdf($namespace,width=12,height=10)
layout(matrix(seq(1,12), 4, 3, byrow = TRUE), 
   widths=rep(1,9), heights=rep(1,9))
par(oma = c(4, 5, 1, 1), mar = c(1, 0, 0, 1))
"""
for i=1:4
    for j=1:3
        R"""
        image(x=$sigtauvec,y=$tauvec,z=$(mpeakadultquant[i,j,:,:]),xlab='Juvenile migration window',ylab='Adult migration window',main='Adult site tooth peaks',col=pal,zlim=c(0,$maxpeakquant))
        """
    end
end
R"""
dev.off()
"""


# SHOWCASE A SINGLE TREATMENT!
filename = "figures/fig_means_peaks2.pdf";
namespace = smartpath(filename);
i = 4; #temp regime
j = 2; #dist regime
mmeanjuv = mean(meanjuv,dims=1)[1,:,:,:,:];
mmeanadult = mean(meanadult,dims=1)[1,:,:,:,:];
minsize = minimum([mmeanjuv[i,j,:,:]; mmeanadult[i,j,:,:]]);
maxsize = maximum([mmeanjuv[i,j,:,:]; mmeanadult[i,j,:,:]]);
R"""
library(RColorBrewer)
library(fields)
pal = colorRampPalette(rev(brewer.pal(11,'Spectral')))(50)
pdf($namespace,width=10,height=8)
layout(matrix(c(1,2,3,4),2,2,byrow=TRUE))
par(oma = c(2, 2, 2, 1), mar = c(4, 5, 2, 5)) #,mai=c(0.6,0.6,0,0.1)
image.plot(x=$sigtauvec,y=$tauvec,z=$(mmeanjuv[i,j,:,:]),xlab='Juvenile migration window',ylab='Adult migration window',col=pal,zlim=c($minsize,$maxsize))
image.plot(x=$sigtauvec,y=$tauvec,z=$(mmeanadult[i,j,:,:]),xlab='Juvenile migration window',ylab='Adult migration window',col=pal,zlim=c($minsize,$maxsize))
"""
mpeakjuvquant = mean(peakjuvquant,dims=1)[1,:,:,:,:];
mpeakadultquant = mean(peakadultquant,dims=1)[1,:,:,:,:];
maxpeakquant = maximum([mpeakjuvquant;mpeakadultquant]);

R"""
pal = c('white',colorRampPalette((brewer.pal(11,'YlGnBu')))(50))
image.plot(x=$sigtauvec,y=$tauvec,z=$(mpeakjuvquant[i,j,:,:]),xlab='Juvenile migration window',ylab='Adult migration window',col=pal,zlim=c(0,$maxpeakquant))
image.plot(x=$sigtauvec,y=$tauvec,z=$(mpeakadultquant[i,j,:,:]),xlab='Juvenile migration window',ylab='Adult migration window',col=pal,zlim=c(0,$maxpeakquant))
dev.off()
"""







# image.plot(x=$sigtauvec,y=$tauvec,z=$(peakadult2[3,1,:,:]),xlab='Juvenile migration window',ylab='Adult migration window',main='Adult site tooth peaks',col=c('white','black'))



#CHECKING VARIANCE
pj = (toothdrop[:,1]/sum(toothdrop[:,1]));
pa = (toothdrop[:,2]/sum(toothdrop[:,2]));


numteethj = Int64.(floor.(pj*100000));
foundteethj = Array{Float64}(undef,sum(numteethj));
let tic = 1
    for i=1:length(numteethj)
        newteethj = repeat([toothlength1[1,i]],outer=numteethj[i]);
        foundteethj[tic:((tic-1) + numteethj[i])] = newteethj;
        tic += numteethj[i];
    end
end

numteetha = Int64.(floor.(pa*100000));
foundteetha = Array{Float64}(undef,sum(numteetha));
let tic = 1
    for i=1:length(numteetha)
        newteetha = repeat([toothlength1[1,i]],outer=numteetha[i]);
        foundteetha[tic:((tic-1) + numteetha[i])] = newteetha;
        tic += numteetha[i];
    end
end
