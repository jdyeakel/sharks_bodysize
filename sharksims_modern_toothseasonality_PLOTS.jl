if homedir() == "/home/z840"
    loadfunc = include("$(homedir())/sharks_bodysize/src/loadfuncs.jl");
else
    loadfunc = include("$(homedir())/Dropbox/PostDoc/2018_sharks/src/loadfuncs.jl");
end

filename_settings = "data/sharks_modern2_toothloss_maxwarm/simsettings.jld";
# filename_settings = "data/sharks_modern2_toothloss_maxwarm/simsettings.jld";
namespace = smartpath(filename_settings);
@load namespace l0 L n0 gen mintemp_j maxtemp_j mintemp_a maxtemp_a distvec toothlossratepercentchangevec sigtauvec tauvec ltl paramvec its


meanjuv = SharedArray{Float64}(ltl,length(sigtauvec),length(tauvec));
meanadult = SharedArray{Float64}(ltl,length(sigtauvec),length(tauvec));


# varjuv = SharedArray{Float64}(reps,length(sigtauvec),length(tauvec));
# varadult = SharedArray{Float64}(reps,length(sigtauvec),length(tauvec));
# peakjuv = SharedArray{Int64}(reps,length(sigtauvec),length(tauvec));
# peakadult = SharedArray{Int64}(reps,length(sigtauvec),length(tauvec));
# peakadult2 = SharedArray{Int64}(reps,length(sigtauvec),length(tauvec));

# peakjuvquant = SharedArray{Float64}(reps,length(sigtauvec),length(tauvec));
# peakadultquant = SharedArray{Float64}(reps,length(sigtauvec),length(tauvec));




@time @sync @distributed for i=1:its
    # println(i)
    #set parameters
   #set parameters
   #set parameters
   pos = paramvec[i,:];
   toothlossratepercentchange_pos = pos[1];
   # temp_pos = pos[2];
   # dist_pos = pos[3];
   sigtau_pos = pos[2];
   tau_pos = pos[3];
   

   #Does the file exist? If not, continue
   

    
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
    # dist = distvec[dist_pos];
    distance = distvec*1000; #3.779e6/10; #1500
    #Shark velocity (m/s)
    velocity = 1;
    D = 1;

    toothlossratepercentchange = toothlossratepercentchangevec[toothlossratepercentchange_pos];
    #Steepness of juvenile migration (function of mass)
    sigtau = sigtauvec[sigtau_pos];
    #Steepness of adult migration (function of time)
    tau = tauvec[tau_pos];

    #save data
    filename = "data/sharks_modern2_toothloss_maxwarm/simdata.jld";
    # filename = "data/sharks_modern2_toothloss_maxwarm/simdata.jld";
    indices = [toothlossratepercentchange_pos,sigtau_pos,tau_pos];
    namespace = smartpath(filename,indices);
    
    @load namespace mass1 mass2 toothdrop toothlength1 toothlength2;
    toothlength = toothlength1[1,:];
    
    mvj, mva, varj, vara, pj, pa, peakjuvbin, peakadultbin, peakjuvdist, peakadultdist, modesj, modesa = toothdist_analysis(toothdrop,toothlength);

    meanjuv[toothlossratepercentchange_pos,sigtau_pos,tau_pos] = mvj;
    meanadult[toothlossratepercentchange_pos,sigtau_pos,tau_pos] = mva;

    # varjuv[r,sigtau_pos,tau_pos] = varj;
    # varadult[r,sigtau_pos,tau_pos] = vara;

    # peakjuv[r,sigtau_pos,tau_pos] = peakjuvbin;
    # peakadult[r,sigtau_pos,tau_pos] = peakadultbin;

    # peakjuvquant[r,sigtau_pos,tau_pos] = peakjuvdist;
    # peakadultquant[r,sigtau_pos,tau_pos] = peakadultdist;

end

#Calculate matrix differencese across distance
matrixdiff_j = Array{Float64}(undef,ltl);
matrixdiff_a = Array{Float64}(undef,ltl);
for i=1:ltl
    matrixdiff_j[i] = (mean(vec(meanjuv[(i),:,:] ./ meanjuv[1,:,:])) .- 1)*100;
    matrixdiff_a[i] = (mean(vec(meanadult[(i),:,:] ./ meanadult[1,:,:])) .- 1)*100;
    # matrixdiff_j[i] = sum(sqrt.((meanjuv[7,:,:] .- meanjuv[(i),:,:]).^2));
    # matrixdiff_a[i] = sum(sqrt.((meanadult[7,:,:] .- meanadult[(i),:,:]).^2));
end

# distvec[1] .- distvec

ply = lineplot((toothlossratepercentchangevec .- toothlossratepercentchangevec[1]),matrixdiff_j)
lineplot!(ply,(toothlossratepercentchangevec .- toothlossratepercentchangevec[1]),7matrixdiff_a)



# Eocene Figure
filename = "figures/fig_toothrateseasonality_modern_maxwarm.pdf";
# filename = "figures/fig_toothrateseasonality_modern_maxwarm.pdf";
namespace = smartpath(filename);
R"""
library(RColorBrewer)
library(fields)
pal = colorRampPalette(rev(brewer.pal(11,'Spectral')))(50)
pdf($namespace,width=6,height=5)
plot($(toothlossratepercentchangevec[1:11] .- toothlossratepercentchangevec[1])*100,$(matrixdiff_j[1:11]),type='l',lty=2,xlab='Loss rate seasonality (%)',ylab='Mean percent change (%)',xlim=c(0,10),ylim=c(-0.5,0.5))
lines($(toothlossratepercentchangevec[1:11] .- toothlossratepercentchangevec[1])*100,$(matrixdiff_a[1:11]),lty = 1)
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
