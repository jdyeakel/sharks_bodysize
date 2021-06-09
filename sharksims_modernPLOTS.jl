if homedir() == "/home/z840"
    loadfunc = include("$(homedir())/sharks_bodysize/src/loadfuncs.jl");
else
    loadfunc = include("$(homedir())/Dropbox/PostDoc/2018_sharks/src/loadfuncs.jl");
end

filename_settings = "data/sharks_modern/simsettings.jld";
namespace = smartpath(filename_settings);
@load namespace l0 L n0 gen mintemp_j maxtemp_j mintemp_a maxtemp_a distvec sigtauvec tauvec reps paramvec its


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
    filename = "data/sharks_modern/simdata.jld";
    indices = [r,sigtau_pos,tau_pos];
    namespace = smartpath(filename,indices);
    
    @load namespace mass1 mass2 toothdrop toothlength1 toothlength2;
    toothlength = toothlength1[1,:];
    
    mvj, mva, varj, vara, pj, pa, peakjuvbin, peakadultbin, peakjuvdist, peakadultdist, modesj, modesa = toothdist_analysis(toothdrop,toothlength);

    meanjuv[r,sigtau_pos,tau_pos] = mvj;
    meanadult[r,sigtau_pos,tau_pos] = mva;

    varjuv[r,sigtau_pos,tau_pos] = varj;
    varadult[r,sigtau_pos,tau_pos] = vara;

    peakjuv[r,sigtau_pos,tau_pos] = peakjuvbin;
    peakadult[r,sigtau_pos,tau_pos] = peakadultbin;

    peakjuvquant[r,sigtau_pos,tau_pos] = peakjuvdist;
    peakadultquant[r,sigtau_pos,tau_pos] = peakadultdist;

end

# Eocene Figure
filename = "figures/fig_means_peaks2_modern.pdf";
namespace = smartpath(filename);
i = 4; #temp regime
j = 2; #dist regime
mmeanjuv = mean(meanjuv,dims=1)[1,:,:];
mmeanadult = mean(meanadult,dims=1)[1,:,:];
minsize = minimum([mmeanjuv[:,:]; mmeanadult[:,:]]);
maxsize = maximum([mmeanjuv[:,:]; mmeanadult[:,:]]);
R"""
library(RColorBrewer)
library(fields)
pal = colorRampPalette(rev(brewer.pal(11,'Spectral')))(50)
pdf($namespace,width=10,height=8)
layout(matrix(c(1,2,3,4),2,2,byrow=TRUE))
par(oma = c(2, 2, 2, 1), mar = c(4, 5, 2, 5)) #,mai=c(0.6,0.6,0,0.1)
image.plot(x=$sigtauvec,y=$tauvec,z=$(mmeanjuv[:,:]),xlab='Juvenile migration window',ylab='Adult migration window',col=pal,zlim=c($minsize,$maxsize))
image.plot(x=$sigtauvec,y=$tauvec,z=$(mmeanadult[:,:]),xlab='Juvenile migration window',ylab='Adult migration window',col=pal,zlim=c($minsize,$maxsize))
"""
mpeakjuvquant = mean(peakjuvquant,dims=1)[1,:,:];
mpeakadultquant = mean(peakadultquant,dims=1)[1,:,:];
maxpeakquant = maximum([mpeakjuvquant;mpeakadultquant]);

R"""
pal = c('white',colorRampPalette((brewer.pal(11,'YlGnBu')))(50))
image.plot(x=$sigtauvec,y=$tauvec,z=$(mpeakjuvquant[:,:]),xlab='Juvenile migration window',ylab='Adult migration window',col=pal,zlim=c(0,$maxpeakquant))
image.plot(x=$sigtauvec,y=$tauvec,z=$(mpeakadultquant[:,:]),xlab='Juvenile migration window',ylab='Adult migration window',col=pal,zlim=c(0,$maxpeakquant))
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
