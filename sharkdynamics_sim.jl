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
n0=100;
gen=20;

#Temperature range
tempmin = 280; tempmax = 300;
tempvec = collect(tempmin:((tempmax-tempmin)/(100-1)):tempmax);

#How many iterations to save to calculate steady state
savebin=1000;


popstate,
savestate,
toothdrop = popgen(m0,M,tempvec,aalpha,bbeta,maturity,n0,savebin,gen);
