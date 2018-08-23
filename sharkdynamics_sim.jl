@everywhere using RCall
@everywhere using Distributions
@everywhere include("$(homedir())/Dropbox/PostDoc/2018_sharks/src/popgen.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2018_sharks/src/ts.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2018_sharks/src/F.jl")


#Body size at birth and maturity
m0 = 136.0;
M = 500000.0;

#Reproduction params
aalpha = 2.0;
bbeta = 0.1;
maturity=1.0;

#Sim params
n0=1000;
gen=20;

#Temperature range
tempmin = 8+273.15; tempmax = 13+273.15;
tempvec = Array{Float64}();
if tempmin == tempmax
    tempvec = repeat([tempmin],inner=100);
else
    tempvec = collect(tempmin:((tempmax-tempmin)/(100-1)):tempmax);
end

#How many iterations to save to calculate steady state
savebin=1000;

epsilonvec,
popstate,
savestate,
toothdrop = popgen(m0,M,tempvec,aalpha,bbeta,maturity,n0,savebin,gen);


namespace = string("$(homedir())/Dropbox/PostDoc/2018_sharks/figures/toothdrop.pdf");
R"""
pdf($namespace,height=5,width=6)
barplot($(toothdrop),names.arg=$(epsilonvec*M)/1000,xlab='Teeth by body mass (Kg)');
dev.off()
"""


R"""
plot($popstate)
"""

accumstate = sum(savestate,1);
laststate = find(iszero,accumstate)[1];
R"""
par=mfrow=c(2,1)
barplot($(accumstate[1:laststate]),names.arg=$((epsilonvec*M)[1:laststate])/1000);
"""
