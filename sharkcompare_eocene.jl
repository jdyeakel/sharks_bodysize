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

#Trim data so only the lowlatitude sites are analyzed
#1 Banks Island (3)
#2 Antarctica = la meseta (4)

data = data[:,3:4];
num = length(names(data));

filename_settings = "data/sharks_eocene2/simsettings.jld";
namespace = smartpath(filename_settings);
@load namespace l0 L n0 gen mintemp_j maxtemp_j mintemp_a maxtemp_a distvec sigtauvec tauvec reps paramvec its

#Convert sigtauvec to mass units
M = (0.00013*L^2.4)*1000;
sigtauvecmass = ((sigtauvec .* M) ./ 50) ./ 1000; #in KG


# reps = 25;
# sigtauvec = collect(0.5:0.5:25); # collect(0.5:0.5:25)
# tauvec = collect(1:1:50); #collect(1:1:50);
# tauits = length(sigtauvec)*length(tauvec);

mchij = SharedArray{Float64}(num,reps,length(sigtauvec),length(tauvec));
mchia = SharedArray{Float64}(num,reps,length(sigtauvec),length(tauvec));
modchij = SharedArray{Float64}(num,reps,length(sigtauvec),length(tauvec));
modchia = SharedArray{Float64}(num,reps,length(sigtauvec),length(tauvec));
moddistchij = SharedArray{Float64}(num,reps,length(sigtauvec),length(tauvec));
moddistchia = SharedArray{Float64}(num,reps,length(sigtauvec),length(tauvec));
sddistchij = SharedArray{Float64}(num,reps,length(sigtauvec),length(tauvec));
sddistchia = SharedArray{Float64}(num,reps,length(sigtauvec),length(tauvec));

mediandistchij = SharedArray{Float64}(num,reps,length(sigtauvec),length(tauvec));
mediandistchia = SharedArray{Float64}(num,reps,length(sigtauvec),length(tauvec));
quartile25distchij = SharedArray{Float64}(num,reps,length(sigtauvec),length(tauvec));
quartile25distchia = SharedArray{Float64}(num,reps,length(sigtauvec),length(tauvec));
quartile75distchij = SharedArray{Float64}(num,reps,length(sigtauvec),length(tauvec));
quartile75distchia = SharedArray{Float64}(num,reps,length(sigtauvec),length(tauvec));

@time @sync @distributed for i=1:num

    measures = Array{Float64}(data[!,i][findall(!ismissing,data[!,i])]);
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
        filename = "data/sharks_eocene2/simdata.jld";
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

        meanchij, meanchia, modechij, modechia, modedistchij, modedistchia, sdchij, sdchia, medianchij, medianchia, quartile25chij, quartile25chia, quartile75chij, quartile75chia  = empirical_sim_comparison(toothdrop,toothlength,measures);

        mchij[i,r,sigtau_pos,tau_pos] = meanchij;
        mchia[i,r,sigtau_pos,tau_pos] = meanchia;
        modchij[i,r,sigtau_pos,tau_pos] = modechij;
        modchia[i,r,sigtau_pos,tau_pos] = modechia;
        moddistchij[i,r,sigtau_pos,tau_pos] = modedistchij;
        moddistchia[i,r,sigtau_pos,tau_pos] = modedistchia;
        sddistchij[i,r,sigtau_pos,tau_pos] = sdchij;
        sddistchia[i,r,sigtau_pos,tau_pos] = sdchia;

        mediandistchij[i,r,sigtau_pos,tau_pos] = medianchij;
        mediandistchia[i,r,sigtau_pos,tau_pos] = medianchia;
        quartile25distchij[i,r,sigtau_pos,tau_pos] = quartile25chij;
        quartile25distchia[i,r,sigtau_pos,tau_pos] = quartile25chia;
        quartile75distchij[i,r,sigtau_pos,tau_pos] = quartile75chij;
        quartile75distchia[i,r,sigtau_pos,tau_pos] = quartile75chia;
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

medj = mean(mediandistchij,dims=2)[:,1,:,:];
meda = mean(mediandistchia,dims=2)[:,1,:,:];
q25j = mean(quartile25distchij,dims=2)[:,1,:,:];
q25a = mean(quartile25distchia,dims=2)[:,1,:,:];
q75j = mean(quartile75distchij,dims=2)[:,1,:,:];
q75a = mean(quartile75distchia,dims=2)[:,1,:,:];

# ndata = names(data);
ndata = ["Banks Island","Seymour Island"];

filename = "data/sharks_eocene2/analysisdata3.jld";
namespace = smartpath(filename);
@save namespace mcj mca mdj mda distj dista medj meda q25j q25a q75j q75a

binmatrixj = Array{Int64}(undef,num,length(sigtauvec),length(tauvec));
binmatrixa = Array{Int64}(undef,num,length(sigtauvec),length(tauvec));
cof = 0.5
for i=1:num
    binmatrixj[i,:,:] = (mcj .< cof)[i,:,:] .* (mdj .< cof)[i,:,:] .* (distj .< cof)[i,:,:] .* (msdj .< cof)[i,:,:] .* (medj .< cof)[i,:,:] .* (q25j .< cof)[i,:,:] .* (q75j .< cof)[i,:,:];
    binmatrixa[i,:,:] = (mca .< cof)[i,:,:] .* (mda .< cof)[i,:,:] .* (dista .< cof)[i,:,:] .* (msda .< cof)[i,:,:] .* (meda .< cof)[i,:,:] .* (q25a .< cof)[i,:,:] .* (q75a .< cof)[i,:,:];
end
qmatrixj = mcj .+ mdj .+ distj .+ msdj .+ medj .+ q25j .+ q75j;
qmatrixa = mca .+ mda .+ dista .+ msda .+ meda .+ q25a .+ q75a;
bfcoordsj = Array{Float64}(undef,num,2);
bfcoordsa = Array{Float64}(undef,num,2);
bfvaluej = Array{Float64}(undef,num);
bfvaluea = Array{Float64}(undef,num);
for i=1:num
    cartj = findmin(qmatrixj[i,:,:]); carta = findmin(qmatrixa[i,:,:]);
    bfvaluej[i] = cartj[1];
    bfvaluea[i] = carta[1];
    coordsj = [sigtauvecmass[cartj[2][1]],tauvec[cartj[2][2]]]; coordsa = [sigtauvecmass[carta[2][1]],tauvec[carta[2][2]]];
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


filename = "figures/fig_empirical_comp_highlatitude3_rev2.pdf";
namespace = smartpath(filename);
filename_data = "data/sharks_eocene2/simdata.jld";
Mj = binmatrixj[1,:,:]; Ma = binmatrixa[1,:,:];
qMj = qmatrixj[1,:,:]; qMa = qmatrixa[1,:,:];
zmin = minimum([qmatrixj[1,:,:];qmatrixa[1,:,:]]);
zmax = maximum([qmatrixj[1,:,:];qmatrixa[1,:,:]]);
measures = Array{Float64}(data[!,1][findall(!ismissing,data[!,1])]);
r=1;
datadensityj, toothlengthj, scaledsimdensityj = plotcompare(Mj,qMj,filename_data,measures,r,"juv");
datadensitya, toothlengtha, scaledsimdensitya = plotcompare(Ma,qMa,filename_data,measures,r,"adult");
R"""
library(fields)
library(RColorBrewer)
pal = brewer.pal(5,'Set1')
pal = c('#6E7F40','#3251E8')
palq = brewer.pal(10,'Spectral')
ncol = c('black','black','black','white','white')
pdf($namespace,width=12,height=5)
par(mfrow=c($num,4))
par(list(oma = c(3, 4, 0, 0), mar = c(1, 1, 1, 2)))
image(x=$sigtauvecmass,y=$tauvec,z=($(qmatrixj[1,:,:])),zlim=c($zmin,$zmax),col=palq,xlab='',ylab='',xaxt='n',yaxt='n')
points($(bfcoordsj[1,1]),$(bfcoordsj[1,2]),pch=21,col='white',bg=pal[1],cex=5)
# text(12.5,48,paste($(ndata[1]),': ',round($(bfvaluej[1]),2),sep=''),col=ncol[1])
axis(side=2,at =NULL,mgp=c(3, 0.75, 0),las=2)
mtext('Adult migration window',side=2,outer=TRUE,adj=0.5,padj=-1.5,cex=1.2)

image(x=$sigtauvecmass,y=$tauvec,z=($(qmatrixa[1,:,:])),zlim=c($zmin,$zmax),col=palq,xlab='',ylab='',xaxt='n',yaxt='n')
points($(bfcoordsa[1,1]),$(bfcoordsa[1,2]),pch=21,col='white',bg=pal[1],cex=5)
# text(12.5,48,paste($(ndata[1]),': ',round($(bfvaluea[1]),2),sep=''),col=ncol[1])
axis(side=2,at =NULL,mgp=c(3, 0.75, 0),las=2)

plot($(datadensityj.x),$(datadensityj.density/maximum(datadensityj.density)),type='l',xlab='',ylab='',col=pal[1],lwd=3,xlim=c(0,40),xaxt='n',yaxt='n')
text(42, 0.92, paste($(ndata[1]),"\n","Error = ",round($(bfvaluej[1]),2),'*',sep=''), pos = 2, cex=1.2)
# axis(side=2,at =NULL,mgp=c(3, 0.75, 0),las=2)
mtext('Scaled density',side=2,outer=TRUE,adj=0.5,padj=41.5,cex=1.2)
lines($toothlengthj,$scaledsimdensityj,lty=1,col='#00000020')
"""
for r=2:reps
    datadensityj, toothlengthj, scaledsimdensityj = plotcompare(Mj,qMj,filename_data,measures,r,"juv");
    R"""
    lines($toothlengthj,$scaledsimdensityj,lty=1,col='#00000020')
    """
end
R"""
plot($(datadensitya.x),$(datadensitya.density/maximum(datadensitya.density)),type='l',xlab='',ylab='',col=pal[1],lwd=3,xlim=c(0,40),xaxt='n',yaxt='n')
lines($toothlengtha,$scaledsimdensitya,lty=1,col='#00000020')
text(42, 0.92, paste($(ndata[1]),"\n","Error = ",round($(bfvaluea[1]),2),sep=''), pos = 2, cex=1.2)
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
    zmin = minimum([qMj;qMa]);
    zmax = maximum([qMj;qMa]);
    measures = Array{Float64}(data[!,i][findall(!ismissing,data[!,i])]);
    r=1;
    datadensityj, toothlengthj, scaledsimdensityj = plotcompare(Mj,qMj,filename_data,measures,r,"juv");
    datadensitya, toothlengtha, scaledsimdensitya = plotcompare(Ma,qMa,filename_data,measures,r,"adult");
    R"""
    image(x=$sigtauvecmass,y=$tauvec,z=($(qMj)),zlim=c($zmin,$zmax),col=palq,xlab='Juvenile migration window',ylab='Adult migration window',xaxt='n',yaxt='n')
    # text(12.5,48,paste($(ndata[i]),': ',round($(bfvaluej[i]),2),sep=''),col=ncol[$i])
    """
    if i==num
        R"""
        axis(side=1,at =NULL,mgp=c(3, 0.75, 0))
        axis(side=2,at =NULL,mgp=c(3, 0.75, 0),las=2)
        # mtext('Juvenile migration window',side=1,outer=TRUE,adj=0.18,padj=1.5,cex=1.2)
        """
    end
    R"""
    points($(bfcoordsj[i,1]),$(bfcoordsj[i,2]),pch=21,col='white',bg=pal[$i],cex=5)
    
    image(x=$sigtauvecmass,y=$tauvec,z=($(qMa)),zlim=c($zmin,$zmax),col=palq,xlab='',ylab='',xaxt='n',yaxt='n')
    # text(12.5,48,paste($(ndata[i]),': ',round($(bfvaluea[i]),2),sep=''),col=ncol[$i])
    points($(bfcoordsa[i,1]),$(bfcoordsa[i,2]),pch=21,col='white',bg=pal[$i],cex=5)
    """
    if i==num
        R"""
        axis(side=2,at =NULL,mgp=c(3, 0.75, 0),las=2)
        axis(side=1,at =NULL,mgp=c(3, 0.75, 0))
        # mtext('Juvenile migration window',side=1,outer=TRUE,adj=0.32,padj=-9,cex=1.5)
        """
    end
    R"""
    plot($(datadensityj.x),$(datadensityj.density/maximum(datadensityj.density)),type='l',xlab='',ylab='',col=pal[$i],lwd=3,xlim=c(0,40),xaxt='n',yaxt='n')
    lines($toothlengthj,$scaledsimdensityj,lty=1,col='#00000020')
    text(42, 0.92, paste($(ndata[i]),"\n","Error = ",round($(bfvaluej[i]),2),'*',sep=''), pos = 2,cex=1.2)
    """
    for r=2:reps
        datadensityj, toothlengthj, scaledsimdensityj = plotcompare(Mj,qMj,filename_data,measures,r,"juv");
        R"""
        lines($toothlengthj,$scaledsimdensityj,lty=1,col='#00000020')
        """
    end
    if i==num
        R"""
        axis(side=1,at =NULL,mgp=c(3, 0.75, 0))
        # mtext('Tooth length (mm)',side=1,outer=TRUE,adj=0.78,padj=1.5,cex=1.2)
        """
    end
    R"""
    plot($(datadensitya.x),$(datadensitya.density/maximum(datadensitya.density)),type='l',xlab='',ylab='',col=pal[$i],lwd=3,xlim=c(0,40),xaxt='n',yaxt='n')
    lines($toothlengtha,$scaledsimdensitya,lty=1,col='#00000020')
    text(42, 0.92, paste($(ndata[i]),"\n","Error = ",round($(bfvaluea[i]),2),sep=''), pos = 2,cex=1.2)
    """
    for r=2:reps
        datadensitya, toothlengtha, scaledsimdensitya = plotcompare(Ma,qMa,filename_data,measures,r,"adult");
        R"""
        lines($toothlengtha,$scaledsimdensitya,lty=1,col='#00000020')
        """
    end
    if i==num
        R"""
        axis(side=1,at =NULL,mgp=c(3, 0.75, 0))
        # mtext('Tooth length (mm)',side=1,outer=TRUE,adj=0.32,padj=-9,cex=1.5)
        """
    end
end
R"""
dev.off()
"""

#DOUBLE CHECK THE DIMENSIONS OF QMATRIX
qmatrixat = [repeat([0],20)';qmatrixa]
qmatrixat = [repeat([0],20)';qmatrixat]
qmatrixat = [repeat([0],20)';qmatrixat]
filename = "figures/fig_qmatrix_test.pdf";
namespace = smartpath(filename);
R"""
pdf($namespace,width=5,height=5)
image(x=$(collect(1:size(qmatrixat)[1])),y=$tauvec,z=($(qmatrixat)))
dev.off()
"""

