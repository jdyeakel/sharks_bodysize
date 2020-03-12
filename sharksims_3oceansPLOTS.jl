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
sigtauvec = collect(1:10);
tauvec = collect(5:25);
tauits = length(sigtauvec)*length(tauvec)*length(distvec);


#3 paramters: distance, sigtau, tau
distposvec = repeat(collect(1:3),inner = length(sigtauvec)*length(tauvec));
sigtauposvec = repeat(collect(1:length(sigtauvec)),inner = length(tauvec),outer=length(distvec));
tauposvec = repeat(collect(1:length(tauvec)),outer=length(distvec)*length(sigtauvec));
paramposvec_pre = [distposvec sigtauposvec tauposvec];
#temperature | distance | sigtauvec | tauposvec
paramposvec = [repeat(collect(1:4),inner = size(paramposvec_pre)[1]) repeat(paramposvec_pre,outer=4)];

reps = 10;
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


@time @sync @distributed for i=1:its
    # println(i)
    #set parameters
    pos = paramvec[i,:];
    r = pos[1];
    temp_pos = pos[2];
    dist_pos = pos[3];
    sigtau_pos = pos[4];
    tau_pos = pos[5];
    
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
    filename = "data/sharks_3oceans/simdata.jld";
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
    else
        peakjuv[r,temp_pos,dist_pos,sigtau_pos,tau_pos] = 1;
    end
    if secondpeaka == 0.0
        peakadult[r,temp_pos,dist_pos,sigtau_pos,tau_pos] = 0;
    else
        peakadult[r,temp_pos,dist_pos,sigtau_pos,tau_pos] = 1;
    end
    
end

tempposvec = repeat(collect(1:4),inner=3)
distposvec = repeat(collect(1:3),outer=4)
ymin = 12;
ymax = 20;

#Inset histogram
filename = "figures/fig_mean.pdf";
namespace = smartpath(filename);
R"""
library(RColorBrewer)
pal=brewer.pal(3,'Set1')
pdf($namespace,width=12,height=10)
layout(matrix(seq(1,9), 3, 3, byrow = TRUE), 
   widths=rep(1,9), heights=rep(1,9))
par(oma = c(4, 5, 1, 1), mar = c(1, 0, 0, 1))

xt = c('n','n','n','n','n','n','s','s','s');
yt = c('s','n','n','s','n','n','s','n','n');

#t1 d1
plot(x=$sigtauvec,y=$(meanjuv[1,1,:,1]),pch=16,xlab='Juvenile migration window',ylab='Mean tooth size',ylim=c($ymin,$ymax),xaxt = 'n',yaxt='s',las=1)
lines(x = $sigtauvec,y=$(meanjuv[1,1,:,1]))
points(x=$sigtauvec,y=$(meanjuv[1,1,:,10]),pch=17)
lines(x = $sigtauvec,y=$(meanjuv[1,1,:,10]))
points(x=$sigtauvec,y=$(meanjuv[1,1,:,21]),pch=18)
lines(x = $sigtauvec,y=$(meanjuv[1,1,:,21]))

points(x=$sigtauvec,y=$(meanadult[1,1,:,1]),pch=16,col=pal[1])
lines(x = $sigtauvec,y=$(meanadult[1,1,:,1]),col=pal[1])
points(x=$sigtauvec,y=$(meanadult[1,1,:,10]),pch=17,col=pal[1])
lines(x = $sigtauvec,y=$(meanadult[1,1,:,10]),col=pal[1])
points(x=$sigtauvec,y=$(meanadult[1,1,:,21]),pch=18,col=pal[1])
lines(x = $sigtauvec,y=$(meanadult[1,1,:,21]),col=pal[1])
"""
for i=2:9
    R"""
    plot(x=$sigtauvec,y=$(meanjuv[tempposvec[i],distposvec[i],:,1]),pch=16,xlab='Juvenile migration window',ylab='Mean tooth size',ylim=c($ymin,$ymax),,xaxt = xt[$i],yaxt=yt[$i],las=1)
    lines(x = $sigtauvec,y=$(meanjuv[tempposvec[i],distposvec[i],:,1]))
    points(x=$sigtauvec,y=$(meanjuv[tempposvec[i],distposvec[i],:,10]),pch=17)
    lines(x = $sigtauvec,y=$(meanjuv[tempposvec[i],distposvec[i],:,10]))
    points(x=$sigtauvec,y=$(meanjuv[tempposvec[i],distposvec[i],:,21]),pch=18)
    lines(x = $sigtauvec,y=$(meanjuv[tempposvec[i],distposvec[i],:,21]))

    points(x=$sigtauvec,y=$(meanadult[tempposvec[i],distposvec[i],:,1]),pch=16,col=pal[1])
    lines(x = $sigtauvec,y=$(meanadult[tempposvec[i],distposvec[i],:,1]),col=pal[1])
    points(x=$sigtauvec,y=$(meanadult[tempposvec[i],distposvec[i],:,10]),pch=17,col=pal[1])
    lines(x = $sigtauvec,y=$(meanadult[tempposvec[i],distposvec[i],:,10]),col=pal[1])
    points(x=$sigtauvec,y=$(meanadult[tempposvec[i],distposvec[i],:,21]),pch=18,col=pal[1])
    lines(x = $sigtauvec,y=$(meanadult[tempposvec[i],distposvec[i],:,21]),col=pal[1])
    """
end
R"""
title(xlab='Juvenile migration window',ylab='Mean tooth size',outer=TRUE,line=2.5,cex.lab=2)
dev.off()
"""


filename = "figures/fig_3oceans_means.pdf";
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

filename = "figures/fig_3oceans_peaks.pdf";
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

filename = "figures/fig_mpeaksjuv.pdf";
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

filename = "figures/fig_mpeaksadult.pdf";
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
