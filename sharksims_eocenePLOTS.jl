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


mintemp_j = 12;
maxtemp_j = 24;
mintemp_a = 9;
maxtemp_a = 17;

distvec = 400;
sigtauvec = collect(0.5:0.5:25);
tauvec = collect(1:1:50);
tauits = length(sigtauvec)*length(tauvec);


#3 paramters: distance, sigtau, tau
# distposvec = repeat(collect(1:3),inner = length(sigtauvec)*length(tauvec));
sigtauposvec = repeat(collect(1:length(sigtauvec)),inner = length(tauvec));
tauposvec = repeat(collect(1:length(tauvec)),outer=length(sigtauvec));
paramposvec_pre = [sigtauposvec tauposvec];
#temperature | distance | sigtauvec | tauposvec
# paramposvec = [repeat(collect(1:4),inner = size(paramposvec_pre)[1]) repeat(paramposvec_pre,outer=4)];

reps = 25;
paramvec_pre = repeat(paramposvec_pre,outer = reps);
paramvec = [repeat(collect(1:reps),inner=size(paramposvec_pre)[1]) paramvec_pre];

its = size(paramvec)[1];


meanjuv = SharedArray{Float64}(reps,length(sigtauvec),length(tauvec));
meanadult = SharedArray{Float64}(reps,length(sigtauvec),length(tauvec));
varjuv = SharedArray{Float64}(reps,length(sigtauvec),length(tauvec));
varadult = SharedArray{Float64}(reps,length(sigtauvec),length(tauvec));
peakjuv = SharedArray{Int64}(reps,length(sigtauvec),length(tauvec));
peakadult = SharedArray{Int64}(reps,length(sigtauvec),length(tauvec));
peakadult2 = SharedArray{Int64}(reps,length(sigtauvec),length(tauvec));

peakjuvquant = SharedArray{Float64}(reps,length(sigtauvec),length(tauvec));
peakadultquant = SharedArray{Float64}(reps,length(sigtauvec),length(tauvec));

@time @sync @distributed for i=1:its
    # println(i)
    #set parameters
    pos = paramvec[i,:];
    r = pos[1];
    # temp_pos = pos[2];
    # dist_pos = pos[3];
    sigtau_pos = pos[2];
    tau_pos = pos[3];
    
    #Temperature range (high latitude: 8-13; Low latitude 22-30)
    #Juvenile site
    tempmin1 = mintemp_j+273.15; tempmax1 = maxtemp_j+273.15;
    #Adult site
    tempmin2 = mintemp_a+273.15; tempmax2 = maxtemp_a+273.15;

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
    distance = distvec*1000; #3.779e6/10; #1500
    #Shark velocity (m/s)
    velocity = 1;
    D = 1;

    #Steepness of juvenile migration (function of mass)
    sigtau = sigtauvec[sigtau_pos];
    #Steepness of adult migration (function of time)
    tau = tauvec[tau_pos];

    #save data
    filename = "data/sharks_eocene/simdata.jld";
    indices = [r,sigtau_pos,tau_pos];
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
    
    # #Gaussian multi-mixture version
    # #generate data
    # numteetha = Int64.(floor.(pa*100000));
    # foundteetha = Array{Float64}(undef,sum(numteetha));
    # let tic = 1
    #     for i=1:length(numteetha)
    #         newteetha = repeat([toothlength1[1,i]],outer=numteetha[i]) .+ rand(Normal(0,0.5),numteetha[i]);
    #         foundteetha[tic:((tic-1) + numteetha[i])] = newteetha;
    #         tic += numteetha[i];
    #     end
    # end
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
    # R"""
    # library(mclust)
    # x.gmm.1 = Mclust($foundteetha,G=1,verbose=FALSE)
    # x.gmm.2 = Mclust($foundteetha,G=2,verbose=FALSE)
    # #likelihood ratio test for 1 vs. 2 components
    # # likediff = logLik(x.gmm.2)-logLik(x.gmm.1)
    # # pvalue = 1-pchisq(likediff, df=3)
    # #If pvalue is < 0.05, accept 2 modes
    # #If pvalue is > 0.05, accept 1 mode
    # 
    # #or likelihood ratio test?
    # likeratiotest = exp(logLik(x.gmm.1)-logLik(x.gmm.2))
    # #if < c then reject the null gmm.1
    # #if > c then accept the null gmm.1
    # # summary(x.gmm)
    # """
    # @rget likeratiotest
    # 
    # #accept null hypothesis of one mode
    # if likeratiotest > 0.05
    #     peakadult2[r,temp_pos,dist_pos,sigtau_pos,tau_pos] = 0;
    # #else reject null hypothesis of one mode
    # else
    #     peakadult2[r,temp_pos,dist_pos,sigtau_pos,tau_pos] = 1;
    # end
    # #JUVENILE PEAKS
    # 
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
    
end

# Eocene Figure
filename = "figures/fig_means_peaks.pdf";
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







image.plot(x=$sigtauvec,y=$tauvec,z=$(peakadult2[3,1,:,:]),xlab='Juvenile migration window',ylab='Adult migration window',main='Adult site tooth peaks',col=c('white','black'))



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
