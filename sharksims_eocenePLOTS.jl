if homedir() == "/home/z840"
    loadfunc = include("$(homedir())/sharks_bodysize/src/loadfuncs.jl");
else
    loadfunc = include("$(homedir())/Dropbox/PostDoc/2018_sharks/src/loadfuncs.jl");
end

# Highlatitude
# foldername = "sharks_eocene2";
# figfilename = "eocene_highlatitude";

# LowLatitude
# foldername = "sharks_eocene_lowlatitude";
# figfilename = "eocene_lowlatitude";

# Modern
foldername = "sharks_modern2";
figfilename = "modern";

filename_settings = string("data/",foldername,"/simsettings.jld");
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
    filename = string("data/",foldername,"/simdata.jld");
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

filename = string("data/",foldername,"/simdata.jld");
#Upper Left
indices = [1,1,length(tauvec)]
tlj1, dj1 = plotdensityreturn(filename,indices,"juv");
tla1, da1 = plotdensityreturn(filename,indices,"adult");

#Upper Right
indices = [1,length(sigtauvec),length(tauvec)]
tlj2, dj2 = plotdensityreturn(filename,indices,"juv");
tla2, da2 = plotdensityreturn(filename,indices,"adult");

#Lower Left
indices = [1,1,1]
tlj3, dj3 = plotdensityreturn(filename,indices,"juv");
tla3, da3 = plotdensityreturn(filename,indices,"adult");

#Lower Right
indices = [1,length(sigtauvec),1]
tlj4, dj4 = plotdensityreturn(filename,indices,"juv");
tla4, da4 = plotdensityreturn(filename,indices,"adult");



# Eocene Figure
filename = string("figures/fig_means_peaks_",figfilename,".pdf");
namespace = smartpath(filename);
# i = 4; #temp regime
# j = 2; #dist regime
mmeanjuv = mean(meanjuv,dims=1)[1,:,:];
mmeanadult = mean(meanadult,dims=1)[1,:,:];
minsize = minimum([mmeanjuv[:,:]; mmeanadult[:,:]]);
maxsize = maximum([mmeanjuv[:,:]; mmeanadult[:,:]]);
R"""
library(RColorBrewer)
library(wesanderson)
library(fields)
# pal = colorRampPalette(rev(brewer.pal(11,'Spectral')))(50)
# pal <- wes_palette("FantasticFox1", 100, type = "continuous")
pal = colorRampPalette((brewer.pal(9,'YlGnBu')))(50)
palpoints = brewer.pal(9,'Set1'); palpoints[6]=palpoints[9];
pdf($namespace,width=10,height=8)
par(list(oma = c(2, 3, 0, 0), mar = c(1, 3, 1, 2)))
par(list(new=TRUE, plt=c(0.02, 0.35, .63, 0.98)))
image(x=$sigtauvec,y=$tauvec,z=$(mmeanjuv[:,:]),xlab='',ylab='',col=pal,zlim=c($minsize,$maxsize),xaxt='n',yaxt='n')
points(head($sigtauvec,1),tail($tauvec,1),pch=21,bg=palpoints[1],col='black',cex=3)
text(head($sigtauvec,1)*4,tail($tauvec,1)*0.9,'I',cex=1.5)
points(tail($sigtauvec,1),tail($tauvec,1),pch=21,bg=palpoints[2],col='black',cex=3)
text(tail($sigtauvec,1)*(0.95),tail($tauvec,1)*0.9,'II',cex=1.5)
points(head($sigtauvec,1),head($tauvec,1),pch=21,bg=palpoints[3],col='black',cex=3)
text(head($sigtauvec,1)*4.4,head($tauvec,1)*4.4,'III',cex=1.5)
points(tail($sigtauvec,1),head($tauvec,1),pch=21,bg=palpoints[4],col='black',cex=3)
text(tail($sigtauvec,1)*0.95,head($tauvec,1)*4.4,'IV',cex=1.5)
axis(side=2,at =NULL,mgp=c(3, 0.75, 0),las=2)
# title(ylab = 'Adult migration window', cex.lab = 1, line = 0.8)
mtext('Adult migration window',side=2,outer=TRUE,adj=0.6,padj=-1.5,cex=1.5)
par(list(new=TRUE, plt=c(0.36, 0.68, .63, 0.98)))
image.plot(x=$sigtauvec,y=$tauvec,z=$(mmeanadult[:,:]),xlab='',ylab='',col=pal,zlim=c($minsize,$maxsize),yaxt='n',xaxt='n',legend.line=2.5)
points(head($sigtauvec,1),tail($tauvec,1),pch=21,bg=palpoints[5],col='black',cex=3)
text(head($sigtauvec,1)*4,tail($tauvec,1)*0.9,'I',cex=1.5,col='white')
points(tail($sigtauvec,1),tail($tauvec,1),pch=21,bg=palpoints[6],col='black',cex=3)
text(tail($sigtauvec,1)*(0.95),tail($tauvec,1)*0.9,'II',cex=1.5)
points(head($sigtauvec,1),head($tauvec,1),pch=21,bg=palpoints[7],col='black',cex=3)
text(head($sigtauvec,1)*4.4,head($tauvec,1)*4.4,'III',cex=1.5,col='white')
points(tail($sigtauvec,1),head($tauvec,1),pch=21,bg=palpoints[8],col='black',cex=3)
text(tail($sigtauvec,1)*0.95,head($tauvec,1)*4.4,'IV',cex=1.5)

"""
mpeakjuvquant = mean(peakjuvquant,dims=1)[1,:,:];
mpeakadultquant = mean(peakadultquant,dims=1)[1,:,:];
maxpeakquant = maximum([mpeakjuvquant;mpeakadultquant]);

R"""
pal = c('white',colorRampPalette((brewer.pal(9,'YlGnBu')))(50))
par(list(new=TRUE, plt=c(0.02, 0.35, .26, 0.61)))
image(x=$sigtauvec,y=$tauvec,z=$(mpeakjuvquant[:,:]),xlab='',ylab='',col=pal,zlim=c(0,$maxpeakquant),yaxt='n',xaxt='n')
points(head($sigtauvec,1),tail($tauvec,1),pch=21,bg=palpoints[1],col='black',cex=3)
text(head($sigtauvec,1)*4,tail($tauvec,1)*0.9,'I',cex=1.5)
points(tail($sigtauvec,1),tail($tauvec,1),pch=21,bg=palpoints[2],col='black',cex=3)
text(tail($sigtauvec,1)*(0.95),tail($tauvec,1)*0.9,'II',cex=1.5)
points(head($sigtauvec,1),head($tauvec,1),pch=21,bg=palpoints[3],col='black',cex=3)
text(head($sigtauvec,1)*4.4,head($tauvec,1)*4.4,'III',cex=1.5)
points(tail($sigtauvec,1),head($tauvec,1),pch=21,bg=palpoints[4],col='black',cex=3)
text(tail($sigtauvec,1)*0.95,head($tauvec,1)*4.4,'IV',cex=1.5,col='white') #,col='white'
axis(side=1,at =NULL,mgp=c(3, 0.75, 0))
axis(side=2,at =NULL,mgp=c(3, 0.75, 0),las=2)
mtext('Juvenile migration window',side=1,outer=TRUE,adj=0.32,padj=-9,cex=1.5)
mtext('Mean',side=1,outer=TRUE,adj=0.73,padj=-63,cex=1)
mtext(expression(paste(Delta,' mode')),side=1,outer=TRUE,adj=0.74,padj=-39,cex=1)
# title(xlab = 'Juvenile migration window', cex.lab = 1, line = 1.75)
# mtext('Adult migration window',side=2,outer=FALSE,adj=0.5,padj=-4)
par(list(new=TRUE, plt=c(0.36, 0.68, .26, 0.61)))
image.plot(x=$sigtauvec,y=$tauvec,z=$(mpeakadultquant[:,:]),xlab='',ylab='',col=pal,zlim=c(0,$maxpeakquant),yaxt='n',xaxt='n')
points(head($sigtauvec,1),tail($tauvec,1),pch=21,bg=palpoints[5],col='black',cex=3)
text(head($sigtauvec,1)*4,tail($tauvec,1)*0.9,'I',cex=1.5)
points(tail($sigtauvec,1),tail($tauvec,1),pch=21,bg=palpoints[6],col='black',cex=3)
text(tail($sigtauvec,1)*(0.95),tail($tauvec,1)*0.9,'II',cex=1.5)
points(head($sigtauvec,1),head($tauvec,1),pch=21,bg=palpoints[7],col='black',cex=3)
text(head($sigtauvec,1)*4.4,head($tauvec,1)*4.4,'III',cex=1.5)
points(tail($sigtauvec,1),head($tauvec,1),pch=21,bg=palpoints[8],col='black',cex=3)
text(tail($sigtauvec,1)*0.95,head($tauvec,1)*4.4,'IV',cex=1.5,col='white')
axis(side=1,at =NULL,mgp=c(3, 0.75, 0))
# title(xlab = 'Juvenile migration window', cex.lab = 1, line = 1.75)

par(list(new=TRUE, plt=c(0.01, 0.2, .04, 0.18)))
plot($tlj1, $dj1,xlim=c(0,40),type='l',lwd=2,yaxt='n',xaxt='n',ylab='',xlab='',axes=FALSE,col=palpoints[1])
polygon(c($tlj1,rev($tlj1)),c($dj1,rep(0,length($dj1))),col=paste(palpoints[1],50,sep=''),border=NA)
text(0,max($dj1),'I')
axis(side=1,at =NULL,mgp=c(3, 0.5, 0),las=1,cex.axis=0.8)


par(list(new=TRUE, plt=c(0.2, 0.4, .04, 0.18)))
plot($tlj2, $dj2,xlim=c(0,40),type='l',lwd=2,yaxt='n',xaxt='n',ylab='',xlab='',axes=FALSE,col=palpoints[2])
polygon(c($tlj2,rev($tlj2)),c($dj2,rep(0,length($dj2))),col=paste(palpoints[2],50,sep=''),border=NA)
text(0,max($dj2),'II')
axis(side=1,at =NULL,mgp=c(3, 0.5, 0),las=1,cex.axis=0.8)

par(list(new=TRUE, plt=c(0.4, 0.6, .04, 0.18)))
plot($tlj3, $dj3,xlim=c(0,40),type='l',lwd=2,yaxt='n',xaxt='n',ylab='',xlab='',axes=FALSE,col=palpoints[3])
polygon(c($tlj3,rev($tlj3)),c($dj3,rep(0,length($dj3))),col=paste(palpoints[3],50,sep=''),border=NA)
text(0,max($dj3),'III')
axis(side=1,at =NULL,mgp=c(3, 0.5, 0),las=1,cex.axis=0.8)

par(list(new=TRUE, plt=c(0.6, 0.8, .04, 0.18)))
plot($tlj4, $dj4,xlim=c(0,40),type='l',lwd=2,yaxt='n',xaxt='n',ylab='',xlab='',axes=FALSE,col=palpoints[4])
polygon(c($tlj4,rev($tlj4)),c($dj4,rep(0,length($dj4))),col=paste(palpoints[4],50,sep=''),border=NA)
axis(side=1,at =NULL,mgp=c(3, 0.5, 0),las=1,cex.axis=0.8)
text(0,max($dj4),'IV')
mtext('Juvenile site tooth length',side=1,outer=TRUE,adj=0.4,padj=-0.1,cex=-0.5)


par(list(new=TRUE, plt=c(0.77, 0.97, 0.84, 0.98)))
plot($tla1, $da1,xlim=c(0,40),type='l',lwd=2,yaxt='n',xaxt='n',ylab='',xlab='',axes=FALSE,col=palpoints[5])
polygon(c($tla1,rev($tla1)),c($da1,rep(0,length($da1))),col=paste(palpoints[5],50,sep=''),border=NA)
text(40,max($da1),'I')
axis(side=1,at =NULL,mgp=c(3, 0.5, 0),las=1,cex.axis=0.8)

par(list(new=TRUE, plt=c(0.77, 0.97, 0.64, 0.78)))
plot($tla2, $da2,xlim=c(0,40),type='l',lwd=2,yaxt='n',xaxt='n',ylab='',xlab='',axes=FALSE,col=palpoints[6])
polygon(c($tla2,rev($tla2)),c($da2,rep(0,length($da2))),col=paste(palpoints[6],50,sep=''),border=NA)
text(40,max($da2),'II')
axis(side=1,at =NULL,mgp=c(3, 0.5, 0),las=1,cex.axis=0.8)

par(list(new=TRUE, plt=c(0.77, 0.97, 0.44, 0.58)))
plot($tla3, $da3,xlim=c(0,40),type='l',lwd=2,yaxt='n',xaxt='n',ylab='',xlab='',axes=FALSE,col=palpoints[7])
polygon(c($tla3,rev($tla3)),c($da3,rep(0,length($da3))),col=paste(palpoints[7],50,sep=''),border=NA)
text(40,max($da3),'III')
axis(side=1,at =NULL,mgp=c(3, 0.5, 0),las=1,cex.axis=0.8)

par(list(new=TRUE, plt=c(0.77, 0.97, 0.24, 0.38)))
plot($tla4, $da4,xlim=c(0,40),type='l',lwd=2,yaxt='n',xaxt='n',ylab='',axes=FALSE,col=palpoints[8],xlab='')
polygon(c($tla4,rev($tla4)),c($da4,rep(0,length($da4))),col=paste(palpoints[8],50,sep=''),border=NA)
text(40,max($da4),'IV')
title(xlab='Adult site tooth length',line=1.4)
axis(side=1,at =NULL,mgp=c(3, 0.5, 0),las=1,cex.axis=0.8)


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
