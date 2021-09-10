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

# filename = "SandTiger_all.csv";
# namespace = string("$(homedir())/Dropbox/PostDoc/2018_sharks/",filename);
# namespace = string("$(homedir())/sharks_bodysize/",filename);
# data = CSV.read(namespace,header=true,DataFrame)

filename_settings = "data/sharks_modern/simsettings.jld";
namespace = smartpath(filename_settings);
@load namespace l0 L n0 gen mintemp_j maxtemp_j mintemp_a maxtemp_a distvec sigtauvec tauvec reps paramvec its

mchij = SharedArray{Float64}(reps,length(sigtauvec),length(tauvec));
mchia = SharedArray{Float64}(reps,length(sigtauvec),length(tauvec));
modchij = SharedArray{Float64}(reps,length(sigtauvec),length(tauvec));
modchia = SharedArray{Float64}(reps,length(sigtauvec),length(tauvec));
moddistchij = SharedArray{Float64}(reps,length(sigtauvec),length(tauvec));
moddistchia = SharedArray{Float64}(reps,length(sigtauvec),length(tauvec));
sddistchij = SharedArray{Float64}(reps,length(sigtauvec),length(tauvec));
sddistchia = SharedArray{Float64}(reps,length(sigtauvec),length(tauvec));


measures = Array{Float64}(data[!,5][findall(!ismissing,data[!,5])]);
# U = kde(measures);
# lineplot(U.x,U.density)

#Now walk through the simulation results




#3 paramters: distance, sigtau, tau
# distposvec = repeat(collect(1:3),inner = length(sigtauvec)*length(tauvec));
sigtauposvec = repeat(collect(1:length(sigtauvec)),inner = length(tauvec));
tauposvec = repeat(collect(1:length(tauvec)),outer=length(sigtauvec));
paramposvec_pre = [sigtauposvec tauposvec];
#temperature | distance | sigtauvec | tauposvec
# paramposvec = [repeat(collect(1:4),inner = size(paramposvec_pre)[1]) repeat(paramposvec_pre,outer=4)];

paramvec_pre = repeat(paramposvec_pre,outer = reps);
paramvec = [repeat(collect(1:reps),inner=size(paramposvec_pre)[1]) paramvec_pre];

its = size(paramvec)[1];

@time @sync @distributed for j=1:its
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
    filename = "data/sharks_modern/simdata.jld";
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

    meanchij, meanchia, modechij, modechia, modedistchij, modedistchia, sdchij, sdchia  = empirical_sim_comparison(toothdrop,toothlength,measures)

    mchij[r,sigtau_pos,tau_pos] = meanchij;
    mchia[r,sigtau_pos,tau_pos] = meanchia;
    modchij[r,sigtau_pos,tau_pos] = modechij;
    modchia[r,sigtau_pos,tau_pos] = modechia;
    moddistchij[r,sigtau_pos,tau_pos] = modedistchij;
    moddistchia[r,sigtau_pos,tau_pos] = modedistchia;
    sddistchij[r,sigtau_pos,tau_pos] = sdchij;
    sddistchia[r,sigtau_pos,tau_pos] = sdchia;
end



#take means across reps
mcj = mean(mchij,dims=1)[1,:,:];
mca = mean(mchia,dims=1)[1,:,:];
mdj = mean(modchij,dims=1)[1,:,:];
mda = mean(modchia,dims=1)[1,:,:];
distj = mean(moddistchij,dims=1)[1,:,:];
dista = mean(moddistchia,dims=1)[1,:,:];
msdj = mean(sddistchij,dims=2)[1,:,:];
msda = mean(sddistchia,dims=2)[1,:,:];

ndata = names(data);

filename = "data/sharks_modern/analysisdata2.jld";
namespace = smartpath(filename);
@save namespace mcj mca mdj mda distj dista

binmatrixj = Array{Int64}(undef,length(sigtauvec),length(tauvec));
binmatrixa = Array{Int64}(undef,length(sigtauvec),length(tauvec));
cof = 0.5

binmatrixj = (mcj .< cof) .* (mdj .< cof) .* (distj .< cof) .* (msdj .< cof);
binmatrixa = (mca .< cof) .* (mda .< cof) .* (dista .< cof) .* (msda .< cof);

qmatrixj = mcj .+ mdj .+ distj .+ msdj;
qmatrixa = mca .+ mda .+ dista .+ msda;

bfcoordsj = Array{Float64}(undef,2);
bfcoordsa = Array{Float64}(undef,2);
cartj = findmin(qmatrixj); carta = findmin(qmatrixa);
bfvaluej = cartj[1];
bfvaluea = carta[1];
coordsj = [sigtauvec[cartj[2][1]],tauvec[cartj[2][2]]]; coordsa = [sigtauvec[carta[2][1]],tauvec[carta[2][2]]];
bfcoordsj = coordsj; bfcoordsa = coordsa;


##################
M = binmatrixj;
qM = qmatrixj
filename_data = "data/sharks_modern/simdata.jld";
measures = Array{Float64}(data[!,5][findall(!ismissing,data[!,5])]);
heatmap(M')
r=1;
site_type = "juv";
datadensity, toothlength, scaledsimdensity = plotcompare(M,qM,filename_data,measures,r,site_type);
ply = lineplot(datadensity.x,datadensity.density/maximum(datadensity.density))
lineplot!(ply,toothlength,scaledsimdensity,color=:red)
###################


filename = "figures/fig_empirical_comp_modern2.pdf";
namespace = smartpath(filename);
filename_data = "data/sharks_modern/simdata.jld";
Mj = binmatrixj[:,:]; Ma = binmatrixa[:,:];
qMj = qmatrixj[:,:]; qMa = qmatrixa[:,:];
zmin = minimum([qMj;qMa]);
zmax = maximum([qMj;qMa]);
measures = Array{Float64}(data[!,5][findall(!ismissing,data[!,5])]);
r=1;
datadensityj, toothlengthj, scaledsimdensityj = plotcompare(Mj,qMj,filename_data,measures,r,"juv");
datadensitya, toothlengtha, scaledsimdensitya = plotcompare(Ma,qMa,filename_data,measures,r,"adult");
R"""
library(fields)
library(RColorBrewer)
pal = brewer.pal(5,'Set1')
ncol = c('black','black','black','white','white')
palq = brewer.pal(11,'Spectral')
pdf($namespace,width=12,height=3)
par(mfrow=c(1,4))
image(x=$sigtauvec,y=$tauvec,z=($(qmatrixj[:,:])),zlim=c($zmin,$zmax),col=palq,xlab='Juvenile migration window',ylab='Adult migration window',main='Juvenile site')
points($(bfcoordsj[1]),$(bfcoordsj[2]),pch=21,col='white',bg=pal[5],cex=2)
text(10,36,paste($(ndata[5]),': ',round($(bfvaluej),2),sep=''),col=ncol[5])
image(x=$sigtauvec,y=$tauvec,z=($(qmatrixa[:,:])),zlim=c($zmin,$zmax),col=palq,xlab='Juvenile migration window',ylab='Adult migration window',main='Adult site')
points($(bfcoordsa[1]),$(bfcoordsa[2]),pch=21,col='white',bg=pal[5],cex=2)
text(10,36,paste($(ndata[5]),': ',round($(bfvaluea),2),sep=''),col=ncol[5])
plot($(datadensityj.x),$(datadensityj.density/maximum(datadensityj.density)),type='l',xlab='Tooth length (mm)',ylab='Scaled density',main='Juvenile site',col=pal[5],lwd=2,xlim=c(0,40))
lines($toothlengthj,$scaledsimdensityj,lty=1,col='#00000020')
"""
for r=2:reps
    datadensityj, toothlengthj, scaledsimdensityj = plotcompare(Mj,qMj,filename_data,measures,r,"juv");
    R"""
    lines($toothlengthj,$scaledsimdensityj,lty=1,col='#00000020')
    """
end
R"""
plot($(datadensitya.x),$(datadensitya.density/maximum(datadensitya.density)),type='l',xlab='Tooth length (mm)',ylab='Scaled density',main='Adult site',col=pal[5],lwd=2,xlim=c(0,40))
lines($toothlengtha,$scaledsimdensitya,lty=1,col='#00000020')
"""
for r=2:reps
    datadensitya, toothlengtha, scaledsimdensitya = plotcompare(Ma,qMa,filename_data,measures,r,"adult");
    R"""
    lines($toothlengtha,$scaledsimdensitya,lty=1,col='#00000020')
    """
end
R"dev.off()"


# for i=2:num
#     Mj = binmatrixj[i,:,:]; Ma = binmatrixa[i,:,:];
#     qMj = qmatrixj[i,:,:]; qMa = qmatrixa[i,:,:];
#     measures = Array{Float64}(data[!,i][findall(!ismissing,data[!,i])]);
#     r=1;
#     datadensityj, toothlengthj, scaledsimdensityj = plotcompare(Mj,qMj,filename_data,measures,r);
#     datadensitya, toothlengtha, scaledsimdensitya = plotcompare(Ma,qMa,filename_data,measures,r);
#     R"""
#     image(x=$sigtauvec,y=$tauvec,z=($(Mj)),col=c('white','black'),xlab='Juvenile migration window',ylab='Adult migration window')
#     text(12.5,48,$(ndata[i]),col=ncol[$i])
#     points($(bfcoordsj[i,1]),$(bfcoordsj[i,2]),pch=21,col='white',bg=pal[$i],cex=2)
#     image(x=$sigtauvec,y=$tauvec,z=($(Ma)),col=c('white','black'),xlab='Juvenile migration window',ylab='Adult migration window')
#     text(12.5,48,$(ndata[i]),col=ncol[$i])
#     points($(bfcoordsa[i,1]),$(bfcoordsa[i,2]),pch=21,col='white',bg=pal[$i],cex=2)
#     plot($(datadensityj.x),$(datadensityj.density/maximum(datadensityj.density)),type='l',xlab='Tooth length (mm)',ylab='Scaled density',col=pal[$i],lwd=2)
#     lines($toothlengthj,$scaledsimdensityj,lty=1,col='#00000020')
#     """
#     for r=2:reps
#         datadensityj, toothlengthj, scaledsimdensityj = plotcompare(Mj,qMj,filename_data,measures,r);
#         R"""
#         lines($toothlengthj,$scaledsimdensityj,lty=1,col='#00000020')
#         """
#     end
#     R"""
#     plot($(datadensitya.x),$(datadensitya.density/maximum(datadensitya.density)),type='l',xlab='Tooth length (mm)',ylab='Scaled density',col=pal[$i],lwd=2)
#     lines($toothlengtha,$scaledsimdensitya,lty=1,col='#00000020')
#     """
#     for r=2:reps
#         datadensitya, toothlengtha, scaledsimdensitya = plotcompare(Ma,qMa,filename_data,measures,r);
#         R"""
#         lines($toothlengtha,$scaledsimdensitya,lty=1,col='#00000020')
#         """
#     end
# end
# R"dev.off()"



