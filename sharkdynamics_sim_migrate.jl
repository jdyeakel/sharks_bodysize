using Distributed
@everywhere using RCall
@everywhere using Distributions
@everywhere using LinearAlgebra
@everywhere using SharedArrays
@everywhere include("$(homedir())/Dropbox/PostDoc/2018_sharks/src/popgen_migrate_g.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2018_sharks/src/ts.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2018_sharks/src/F.jl")
R"""
logoing_func<-function(logo, x, y, size){
  dims<-dim(logo)[1:2] #number of x-y pixels for the logo (aspect ratio)
  AR<-dims[1]/dims[2]
  par(usr=c(0, 1, 0, 1))
  rasterImage(logo, x-(size/2), y-(AR*size/2), x+(size/2), y+(AR*size/2), interpolate=TRUE)
}
"""

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

#Temperature range (high latitude: 8-13; Low latitude 22-30)
#Juvenile site
tempmin1 = 11+273.15; tempmax1 = 19+273.15;
#Adult site
tempmin2 = 12+273.15; tempmax2 = 20+273.15;

tempvec1 = Array{Float64}(undef,0);
tempvec2 = Array{Float64}(undef,0);
if tempmin1 == tempmax1
    tempvec1 = repeat([tempmin1],inner=100);
else
    tempvec1 = collect(tempmin1:((tempmax1-tempmin1)/(100-1)):tempmax1);
end
if tempmin2 == tempmax2
    tempvec2 = repeat([tempmin2],inner=100);
else
    tempvec2 = collect(tempmin2:((tempmax2-tempmin2)/(100-1)):tempmax2);
end

#distance between site (km * 1000)
distance = 650*1000; #3.779e6/10; #1500
#Shark velocity (m/s)
velocity = 1;
D = 1;

#Steepness of juvenile migration (function of mass)
sigtau = 8.;
#Steepness of adult migration (function of time)
tau = 20.;


mass1,
mass2,
epsilonvec,
clock,
popstate,
toothdrop,
state = popgen_migrate_g(m0,M,tempvec1,tempvec2,n0,gen,distance,velocity,D,sigtau,tau);
tpop = sum(popstate,dims=2);
toothlength1 = 2.13337 .+ (0.187204 .* mass1.^(0.416667)); #mm
toothlength2 = 2.13337 .+ (0.187204 .* mass2.^(0.416667)); #mm

println(dot(toothlength2[1,:],(toothdrop[:,1]/sum(toothdrop[:,1]))))
println(dot(toothlength2[1,:],(toothdrop[:,2]/sum(toothdrop[:,2]))))



namespace = "$(homedir())/Dropbox/PostDoc/2018_sharks/figures/fig_toothdist_650km_largewindow.pdf";
sharklogo = "$(homedir())/Dropbox/PostDoc/2018_sharks/sharktooth.png";
R"""
library(png)
library(RCurl)
library(RColorBrewer)
shark_logo <-  readPNG($sharklogo);
pal = brewer.pal(3,'Set1')
pdf($namespace, width = 8, height = 4)
layout(matrix(c(1,2), 1, 2, byrow = TRUE), 
             widths=c(1,1), heights=c(1,1))
# par(mfrow=c(2,2))
par(oma = c(2,2,1.5,1.6), mar = c(2,2,1,1))
# par(mfrow=c(1,2))
plot($(vec(toothlength1[1,:])),$(toothdrop[:,1])/sum($(toothdrop[:,1])),type='h',log='',lwd=3,ylim=c(0.0001,max($(toothdrop[:,1])/sum($(toothdrop[:,1])))),xlim=c(0,max(c($(toothlength1[1,:]),$(toothlength2[1,:])))),xlab='Body size',ylab='Proportion (tooth drop)',col=pal[1])
logoing_func(shark_logo, x=0.90, y=0.90, size=0.15)
plot($(vec(toothlength2[1,:])),$(toothdrop[:,2])/sum($(toothdrop[:,2])),type='h',log='',lwd=3,ylim=c(0.0001,max($(toothdrop[:,2])/sum($(toothdrop[:,2])))),xlim=c(0,max(c($(toothlength1[1,:]),$(toothlength2[1,:])))),xlab='Body size',ylab='Proportion (tooth drop)',col=pal[2])
logoing_func(shark_logo, x=0.90, y=0.90, size=0.15)
title(xlab='Tooth height (mm)',ylab='Frequency',outer=TRUE,line=0.70,cex.lab=1.25)
dev.off()
"""






namespace = "$(homedir())/Dropbox/PostDoc/2018_sharks/figures/fig_toothdist_group.pdf";
R"""
library(png)
library(RCurl)
library(RColorBrewer)
shark_logo <-  readPNG($sharklogo);
pal = brewer.pal(3,'Set1')
pdf($namespace, width = 5, height = 5)
par(mfrow=c(1,1))
plot($(vec(toothlength2[1,:])),$(vec(sum(toothdrop,dims=2)))/$(sum(vec(sum(toothdrop,dims=2)))),type='h',log='',lwd=3,ylim=c(0.0001,max($(vec(sum(toothdrop,dims=2)))/$(sum(vec(sum(toothdrop,dims=2)))))),xlim=c(0,max(c($(toothlength1[1,:]),$(toothlength2[1,:])))),xlab='Body size',ylab='Proportion (tooth drop)',col=pal[3])
logoing_func(shark_logo, x=0.90, y=0.90, size=0.15)
dev.off()
"""
println(dot(toothlength2[1,:],(toothdrop[:,1]/sum(toothdrop[:,1]))))
println(dot(toothlength2[1,:],(toothdrop[:,2]/sum(toothdrop[:,2]))))




namespace = "$(homedir())/Dropbox/PostDoc/2018_sharks/figures/fig_population.pdf";
sharklogo = "$(homedir())/Dropbox/PostDoc/2018_sharks/sharktooth.png";
R"""
library(png)
library(RCurl)
library(RColorBrewer)
shark_logo <-  readPNG($sharklogo);
pal = brewer.pal(3,'Set1')
pdf($namespace, width = 12, height = 10)
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE), 
             widths=c(1,1,1,1), heights=c(1,1,1,1))
# par(mfrow=c(2,2))
par(oma = c(2,2,1.5,1.6), mar = c(2,2,1,1))
plot($clock/60/60/24/365,$(popstate[:,1]),type='l',ylim=c(1,max($tpop)),log='',xlab='years',ylab='Population',col=pal[1])
lines($clock/60/60/24/365,$(popstate[:,2]),col=pal[2])
lines($clock/60/60/24/365,$(tpop),col=pal[3])
# par(mfrow=c(1,2))
plot($(vec(toothlength1[1,:])),$(toothdrop[:,1])/sum($(toothdrop[:,1])),type='h',log='',lwd=3,ylim=c(0.0001,max($(toothdrop[:,1])/sum($(toothdrop[:,1])))),xlim=c(0,max(c($(toothlength1[1,:]),$(toothlength2[1,:])))),xlab='Body size',ylab='Proportion (tooth drop)',col=pal[1])
logoing_func(shark_logo, x=0.90, y=0.90, size=0.15)
plot($(vec(toothlength2[1,:])),$(toothdrop[:,2])/sum($(toothdrop[:,2])),type='h',log='',lwd=3,ylim=c(0.0001,max($(toothdrop[:,2])/sum($(toothdrop[:,2])))),xlim=c(0,max(c($(toothlength1[1,:]),$(toothlength2[1,:])))),xlab='Body size',ylab='Proportion (tooth drop)',col=pal[2])
logoing_func(shark_logo, x=0.90, y=0.90, size=0.15)
dev.off()
"""


R"""
plot($(vec(toothlength2[1,:])),$(vec(sum(toothdrop,dims=2)))/$(sum(vec(sum(toothdrop,dims=2)))),type='h',log='x',lwd=5,xlab='Body size',ylab='Proportion (tooth drop)',col=pal[2])
logoing_func(shark_logo, x=0.90, y=0.90, size=0.15)
"""







clockyrs = clock ./ 60 ./ 60 ./ 24 ./ 365;
inds = Array{Int64}(undef,0);
massint = Int64.(round.(mass[1,:]/1000));
for i=1:length(state)
    n = repeat([massint[i]],outer=state[i]);
    append!(inds,n);
end


R"""
library(png)
library(RCurl)
sharkurl<-"http://phylopic.org/assets/images/submissions/545d45f0-0dd1-4cfd-aad6-2b835223ea0d.1024.png"
shark_logo <-  readPNG(getURLContent(sharkurl))
layout(matrix(c(1,2,3), 1, 3, byrow = TRUE), 
   widths=c(1,1,1), heights=c(0.5,0.5,0.5))
par(oma = c(2, 1, 1, 1), mar = c(3, 4, 0, 1))
plot($clockyrs,$popstate,type='l',xlab='',ylab='')
logoing_func(shark_logo, x=0.90, y=0.10, size=0.15)
title(1,xlab='Years',line=2.5,xpd=NA)
title(2,ylab='Shark population',line=2.5,xpd=NA)
barplot($(toothdrop),names.arg=round($(mass[1,:])/1000,0),xlab='');
title(1,xlab='Size class',line=2.5,xpd=NA)
title(2,ylab='Number of teeth',line=2.5,xpd=NA)
hist($(inds),breaks=30,col='gray',main='')
title(1,xlab='Number of sharks',line=2.5,xpd=NA)
#dev.off()
"""


# accumstate = sum(savestate,dims=1);
# laststate = findall(iszero,vec(accumstate))[1];
R"""
par=mfrow=c(2,1)
barplot($(accumstate[1:laststate]),names.arg=round($(mass[1,:])/1000,0));
"""



poprep = "juv";
temp = [8,22];
tempdiff = [1,7];
namespace = string("$(homedir())/Dropbox/PostDoc/2018_sharks/figures/",poprep,"_migrate.pdf");
R"""
library(png)
library(RCurl)
sharkurl<-"http://phylopic.org/assets/images/submissions/545d45f0-0dd1-4cfd-aad6-2b835223ea0d.1024.png"
shark_logo <-  readPNG(getURLContent(sharkurl))
pdf($namespace,width=8,height=7)
layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE), 
   widths=c(1,1,1,1), heights=c(0.5,0.5,0.5,0.5))
par(oma = c(2, 1, 1, 1), mar = c(3, 4, 0, 1))
"""
#iterate
for i=1:2
    for j=1:2
        #Temperature range (high latitude: 8-13; Low latitude 22-30)
        tempmin = temp[i]+273.15; tempmax = (temp[i]+tempdiff[j]) + 273.15;
        tempvec = Array{Float64}(undef,0);
        if tempmin == tempmax
            tempvec = repeat([tempmin],inner=100);
        else
            tempvec = collect(tempmin:((tempmax-tempmin)/(100-1)):tempmax);
        end

        #How many iterations to save to calculate steady state
        savebin=1000;
        
        mass,
        epsilonvec,
        clock,
        popstate,
        savestate,
        toothdrop,
        state = popgen(m0,M,tempvec,n0,savebin,gen,poprep);

        clockyrs = clock ./ 60 ./ 60 ./ 24 ./ 365;
        inds = Array{Int64}(undef,0);
        massint = Int64.(round.(mass[1,:]/1000));
        for i=1:length(state)
            n = repeat([massint[i]],outer=state[i]);
            append!(inds,n);
        end

        R"""
        barplot($(toothdrop),names.arg=round($(mass[1,:])/1000,0),xlab='',ylab='',main='');
        title(xlab='Size class',line=2.5,xpd=NA)
        title(ylab='Number of teeth',line=2.5,xpd=NA)
        logoing_func(shark_logo, x=0.90, y=0.90, size=0.15)
        mtext(paste('mintemp = ',$(temp[i])),side=3,padj=5,adj=1.13)
        mtext(paste('maxtemp = ',$(temp[i]+tempdiff[j])),side=3,padj=6.5,adj=1.15)
        """
    end
end
R"dev.off()"
