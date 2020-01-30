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


mintemp_j = [8,20,9];
maxtemp_j = [13,30,17];
mintemp_a = [8,20,12];
maxtemp_a = [13,30,24];

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
paramposvec = [repeat(collect(1:3),inner = size(paramposvec_pre)[1]) repeat(paramposvec_pre,outer=3)];



its = size(paramposvec)[1];


meanjuv = SharedArray{Float64}(3,3,length(sigtauvec),length(tauvec));
meanadult = SharedArray{Float64}(3,3,length(sigtauvec),length(tauvec));
varjuv = SharedArray{Float64}(3,3,length(sigtauvec),length(tauvec));
varadult = SharedArray{Float64}(3,3,length(sigtauvec),length(tauvec));
peakjuv = SharedArray{Int64}(3,3,length(sigtauvec),length(tauvec));
peakadult = SharedArray{Int64}(3,3,length(sigtauvec),length(tauvec));

@time @sync @distributed for i=1:its
    
    #set parameters
    pos = paramposvec[i,:];
    temp_pos = pos[1];
    dist_pos = pos[2];
    sigtau_pos = pos[3];
    tau_pos = pos[4];
    
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
    indices = [temp_pos,dist_pos,sigtau_pos,tau_pos];
    namespace = smartpath(filename,indices);
    
    @load namespace mass1 mass2 epsilonvec clock popstate toothdrop state toothlength1 toothlength2;
    
    
    #MEASURE SOMETHING
    mvj = dot(toothlength1[1,:],(toothdrop[:,1]/sum(toothdrop[:,1])));
    mva = dot(toothlength1[1,:],(toothdrop[:,2]/sum(toothdrop[:,2])));
    
    meanjuv[temp_pos,dist_pos,sigtau_pos,tau_pos] = mvj;
    meanadult[temp_pos,dist_pos,sigtau_pos,tau_pos] = mva;
    
    #calculate variance
    meanofxsquared = dot(toothlength1[1,:].^2,(toothdrop[:,1]/sum(toothdrop[:,1])));
    squareofmeanofx = mvj^2;
    varj = meanofxsquared - squareofmeanofx;
    varjuv[temp_pos,dist_pos,sigtau_pos,tau_pos] = varj;
    
    meanofxsquared = dot(toothlength1[1,:].^2,(toothdrop[:,2]/sum(toothdrop[:,2])));
    squareofmeanofx = mva^2;
    vara = meanofxsquared - squareofmeanofx;
    varadult[temp_pos,dist_pos,sigtau_pos,tau_pos] = vara;
    
    #Detect multiple modes (yes, no)
    #How many modes, and rank by significance
    #Where are the multiple modes (position of each? two most significant?)
    pj = (toothdrop[:,1]/sum(toothdrop[:,1]));
    pa = (toothdrop[:,2]/sum(toothdrop[:,2]));
    
    lm=findlocalmaxima(pa)
    #choose highest probability peaks
    peakprobs = pa[lm];
    sortedpeaks = pa[lm[sortperm(peakprobs)]];
    maxpeaka = last(sortedpeaks);
    sortedtooth = toothlength1[1,:][lm[sortperm(peakprobs)]];
    maxtootha = last(sortedtooth);
    secondpeaka = 0.0;
    secondtootha = 0.0;
    #need to determine second peak significance
    if length(sortedpeaks) > 1
        #threshold of 25% size of first peak
        if (sortedpeaks[length(sortedpeaks)-1] / maxpeak) > 0.25
            secondpeaka = sortedpeaks[length(sortedpeaks)-1];
            secondtootha = sortedtooth[length(sortedpeaks)-1];
        end
    end
    
    # pl = lineplot(toothlength1[1,:],pa)
    # lineplot!(pl,repeat([maxtooth],outer=2),[0,maxpeak],color=:red)
    # lineplot!(pl,repeat([secondtooth],outer=2),[0,secondpeak],color=:red)
    
    lm=findlocalmaxima(pj)
    #choose highest probability peaks
    peakprobs = pa[lm];
    sortedpeaks = pa[lm[sortperm(peakprobs)]];
    maxpeakj = last(sortedpeaks);
    sortedtooth = toothlength1[1,:][lm[sortperm(peakprobs)]];
    maxtoothj = last(sortedtooth);
    secondpeakj = 0.0;
    secondtoothj = 0.0;
    #need to determine second peak significance
    if length(sortedpeaks) > 1
        #threshold of 25% size of first peak
        if (sortedpeaks[length(sortedpeaks)-1] / maxpeak) > 0.25
            secondpeakj = sortedpeaks[length(sortedpeaks)-1];
            secondtoothj = sortedtooth[length(sortedpeaks)-1];
        end
    end
    
    if secondpeakj == 0.0
        peakjuv[temp_pos,dist_pos,sigtau_pos,tau_pos] = 0;
    else
        peakjuv[temp_pos,dist_pos,sigtau_pos,tau_pos] = 1;
    end
    if secondpeaka == 0.0
        peakadult[temp_pos,dist_pos,sigtau_pos,tau_pos] = 0;
    else
        peakadult[temp_pos,dist_pos,sigtau_pos,tau_pos] = 1;
    end
end

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
