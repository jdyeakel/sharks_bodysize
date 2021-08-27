if homedir() == "/home/z840"
    loadfunc = include("$(homedir())/sharks_bodysize/src/loadfuncs.jl");
    filename = "SandTiger_all.csv";
    namespace = string("$(homedir())/sharks_bodysize/",filename);
    data = CSV.read(namespace,header=true,DataFrame)
else
    loadfunc = include("$(homedir())/Dropbox/PostDoc/2018_sharks/src/loadfuncs.jl");
    filename = "SandTiger_all.csv";
    namespace = string("$(homedir())/Dropbox/PostDoc/2018_sharks/",filename);
    data = CSV.read(namespace,header=true,DataFrame)
end


num = length(names(data));

# filename = "SandTiger_all.csv";
# namespace = string("$(homedir())/Dropbox/PostDoc/2018_sharks/",filename);
# namespace = string("$(homedir())/sharks_bodysize/",filename);
# data = CSV.read(namespace,header=true,DataFrame)

filename_settings = "data/sharks_eocene_lowlatitude/simsettings.jld";
namespace = smartpath(filename_settings);
@load namespace l0 L n0 gen mintemp_j maxtemp_j mintemp_a maxtemp_a distvec sigtauvec tauvec reps paramvec its

mchij = SharedArray{Float64}(num,reps,length(sigtauvec),length(tauvec));
mchia = SharedArray{Float64}(num,reps,length(sigtauvec),length(tauvec));
modchij = SharedArray{Float64}(num,reps,length(sigtauvec),length(tauvec));
modchia = SharedArray{Float64}(num,reps,length(sigtauvec),length(tauvec));
moddistchij = SharedArray{Float64}(num,reps,length(sigtauvec),length(tauvec));
moddistchia = SharedArray{Float64}(num,reps,length(sigtauvec),length(tauvec));
sddistchij = SharedArray{Float64}(num,reps,length(sigtauvec),length(tauvec));
sddistchia = SharedArray{Float64}(num,reps,length(sigtauvec),length(tauvec));

@time @sync @distributed for i=1:num

    measures = Array{Float64}(data[!,i][findall(!ismissing,data[!,i])]);
    # U = kde(measures);
    # lineplot(U.x,U.density)

    #Now walk through the simulation results

    

    # sigtauposvec = repeat(collect(1:length(sigtauvec)),inner = length(tauvec));
    # tauposvec = repeat(collect(1:length(tauvec)),outer=length(sigtauvec));
    # paramposvec_pre = [sigtauposvec tauposvec];

    # paramvec_pre = repeat(paramposvec_pre,outer = reps);
    # paramvec = [repeat(collect(1:reps),inner=size(paramposvec_pre)[1]) paramvec_pre];

    # its = size(paramvec)[1];

    for j=1:its
        pos = paramvec[j,:];
        r = pos[1];
        # temp_pos = pos[2];
        # dist_pos = pos[3];
        sigtau_pos = pos[2];
        tau_pos = pos[3];

        #Steepness of juvenile migration (function of mass)
        sigtau = sigtauvec[sigtau_pos];
        #Steepness of adult migration (function of time)
        tau = tauvec[tau_pos];
    
        #save data
        filename = "data/sharks_eocene_lowlatitude/simdata.jld";
        indices = [r,sigtau_pos,tau_pos];
        namespace = smartpath(filename,indices);
        
        @load namespace mass1 mass2 toothdrop toothlength1 toothlength2;
        toothlength = toothlength1[1,:];
        
        #plot it!
        # simdensity_juv = toothdrop[:,1]./sum(toothdrop[:,1]);
        # simdensity_adult = toothdrop[:,2]./sum(toothdrop[:,2]);
        # U = kde(measures);
        # maxy = maximum([simdensity_juv; simdensity_adult; U.density]);
        # ply = lineplot(toothlength,simdensity_juv,xlim = [0, 40],ylim=[0,maxy],color=:green)
        # lineplot!(ply,toothlength,simdensity_adult,color=:green)
        # lineplot!(ply,U.x,U.density,color=:blue)

        meanchij, meanchia, modechij, modechia, modedistchij, modedistchia, sdchij, sdchia = empirical_sim_comparison(toothdrop,toothlength,measures)

        mchij[i,r,sigtau_pos,tau_pos] = meanchij;
        mchia[i,r,sigtau_pos,tau_pos] = meanchia;
        modchij[i,r,sigtau_pos,tau_pos] = modechij;
        modchia[i,r,sigtau_pos,tau_pos] = modechia;
        moddistchij[i,r,sigtau_pos,tau_pos] = modedistchij;
        moddistchia[i,r,sigtau_pos,tau_pos] = modedistchia;

        sddistchij[i,r,sigtau_pos,tau_pos] = sdchij;
        sddistchia[i,r,sigtau_pos,tau_pos] = sdchia;
    end


end

#take means across reps
mcj = mean(mchij,dims=2)[:,1,:,:];
mca = mean(mchia,dims=2)[:,1,:,:];
mdj = mean(modchij,dims=2)[:,1,:,:];
mda = mean(modchia,dims=2)[:,1,:,:];
distj = mean(moddistchij,dims=2)[:,1,:,:];
dista = mean(moddistchia,dims=2)[:,1,:,:];
msdj = mean(sddistchij,dims=2)[:,1,:,:];
msda = mean(sddistchia,dims=2)[:,1,:,:];

ndata = names(data);

filename = "data/sharks_eocene_lowlatitude/analysisdata.jld";
namespace = smartpath(filename);
@save namespace mcj mca mdj mda distj dista

binmatrixj = Array{Int64}(undef,num,length(sigtauvec),length(tauvec));
binmatrixa = Array{Int64}(undef,num,length(sigtauvec),length(tauvec));
cof = 0.5
for i=1:5
    binmatrixj[i,:,:] = (mcj .< cof)[i,:,:] .* (mdj .< cof)[i,:,:] .* (distj .< cof)[i,:,:] .* (msdj .< cof)[i,:,:];
    binmatrixa[i,:,:] = (mca .< cof)[i,:,:] .* (mda .< cof)[i,:,:] .* (dista .< cof)[i,:,:] .* (msda .< cof)[i,:,:];
end
qmatrixj = mcj .+ mdj .+ distj .+ msdj;
qmatrixa = mca .+ mda .+ dista .+ msda;
bfcoordsj = Array{Float64}(undef,num,2);
bfcoordsa = Array{Float64}(undef,num,2);
for i=1:num
    cartj = findmin(qmatrixj[i,:,:]); carta = findmin(qmatrixa[i,:,:]);
    coordsj = [sigtauvec[cartj[2][1]],tauvec[cartj[2][2]]]; coordsa = [sigtauvec[carta[2][1]],tauvec[carta[2][2]]];
    bfcoordsj[i,:] = coordsj; bfcoordsa[i,:] = coordsa;
end


##################
# i=1; M = binmatrixj[i,:,:];
# filename_data = "data/sharks_eocene2/simdata.jld";
# measures = Array{Float64}(data[!,i][findall(!ismissing,data[!,i])]);
# heatmap(M)
# datadensity, toothlength, scaledsimdensity = plotcompare(M,filename_data,measures);
# ply = lineplot(datadensity.x,datadensity.density/maximum(datadensity.density))
# lineplot!(ply,toothlength,scaledsimdensity,color=:red)
###################


filename = "figures/fig_empirical_comp_lowlatitude2.pdf";
namespace = smartpath(filename);
filename_data = "data/sharks_eocene_lowlatitude/simdata.jld";
Mj = binmatrixj[1,:,:]; Ma = binmatrixa[1,:,:];
qMj = qmatrixj[1,:,:]; qMa = qmatrixa[1,:,:];
measures = Array{Float64}(data[!,1][findall(!ismissing,data[!,1])]);
r=1;
datadensityj, toothlengthj, scaledsimdensityj = plotcompare(Mj,qMj,filename_data,measures,r,"juv");
datadensitya, toothlengtha, scaledsimdensitya = plotcompare(Ma,qMa,filename_data,measures,r,"adult");
R"""
library(fields)
library(RColorBrewer)
pal = brewer.pal(5,'Set1')
ncol = c('black','black','black','white','white')
pdf($namespace,width=12,height=15)
par(mfrow=c(5,4))
image(x=$sigtauvec,y=$tauvec,z=($(binmatrixj[1,:,:])),col=c('white','black'),xlab='Juvenile migration window',ylab='Adult migration window',main='Juvenile site')
points($(bfcoordsj[1,1]),$(bfcoordsj[1,2]),pch=21,col='white',bg=pal[1],cex=2)
text(10,38,$(ndata[1]),col=ncol[1])
image(x=$sigtauvec,y=$tauvec,z=($(binmatrixa[1,:,:])),col=c('white','black'),xlab='Juvenile migration window',ylab='Adult migration window',main='Adult site')
points($(bfcoordsa[1,1]),$(bfcoordsa[1,2]),pch=21,col='white',bg=pal[1],cex=2)
text(12.5,38,$(ndata[1]),col=ncol[1])
plot($(datadensityj.x),$(datadensityj.density/maximum(datadensityj.density)),type='l',xlab='Tooth length (mm)',ylab='Scaled density',main='Juvenile site',col=pal[1],lwd=2,xlim=c(0,40))
lines($toothlengthj,$scaledsimdensityj,lty=1,col='#00000020')
"""
for r=2:reps
    datadensityj, toothlengthj, scaledsimdensityj = plotcompare(Mj,qMj,filename_data,measures,r,"juv");
    R"""
    lines($toothlengthj,$scaledsimdensityj,lty=1,col='#00000020')
    """
end
R"""
plot($(datadensitya.x),$(datadensitya.density/maximum(datadensitya.density)),type='l',xlab='Tooth length (mm)',ylab='Scaled density',main='Adult site',col=pal[1],lwd=2,xlim=c(0,40))
lines($toothlengtha,$scaledsimdensitya,lty=1,col='#00000020')
"""
for r=2:reps
    datadensitya, toothlengtha, scaledsimdensitya = plotcompare(Ma,qMa,filename_data,measures,r,"adult");
    R"""
    lines($toothlengtha,$scaledsimdensitya,lty=1,col='#00000020')
    """
end
for i=2:num
    Mj = binmatrixj[i,:,:]; Ma = binmatrixa[i,:,:];
    qMj = qmatrixj[i,:,:]; qMa = qmatrixa[i,:,:];
    measures = Array{Float64}(data[!,i][findall(!ismissing,data[!,i])]);
    r=1;
    datadensityj, toothlengthj, scaledsimdensityj = plotcompare(Mj,qMj,filename_data,measures,r,"juv");
    datadensitya, toothlengtha, scaledsimdensitya = plotcompare(Ma,qMa,filename_data,measures,r,"adult");
    R"""
    image(x=$sigtauvec,y=$tauvec,z=($(Mj)),col=c('white','black'),xlab='Juvenile migration window',ylab='Adult migration window')
    text(10,38,$(ndata[i]),col=ncol[$i])
    points($(bfcoordsj[i,1]),$(bfcoordsj[i,2]),pch=21,col='white',bg=pal[$i],cex=2)
    image(x=$sigtauvec,y=$tauvec,z=($(Ma)),col=c('white','black'),xlab='Juvenile migration window',ylab='Adult migration window')
    text(10,38,$(ndata[i]),col=ncol[$i])
    points($(bfcoordsa[i,1]),$(bfcoordsa[i,2]),pch=21,col='white',bg=pal[$i],cex=2)
    plot($(datadensityj.x),$(datadensityj.density/maximum(datadensityj.density)),type='l',xlab='Tooth length (mm)',ylab='Scaled density',col=pal[$i],lwd=2,xlim=c(0,40))
    lines($toothlengthj,$scaledsimdensityj,lty=1,col='#00000020')
    """
    for r=2:reps
        datadensityj, toothlengthj, scaledsimdensityj = plotcompare(Mj,qMj,filename_data,measures,r,"juv");
        R"""
        lines($toothlengthj,$scaledsimdensityj,lty=1,col='#00000020')
        """
    end
    R"""
    plot($(datadensitya.x),$(datadensitya.density/maximum(datadensitya.density)),type='l',xlab='Tooth length (mm)',ylab='Scaled density',col=pal[$i],lwd=2,xlim=c(0,40))
    lines($toothlengtha,$scaledsimdensitya,lty=1,col='#00000020')
    """
    for r=2:reps
        datadensitya, toothlengtha, scaledsimdensitya = plotcompare(Ma,qMa,filename_data,measures,r,"adult");
        R"""
        lines($toothlengtha,$scaledsimdensitya,lty=1,col='#00000020')
        """
    end
end
R"dev.off()"



