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
tempmin = 17+273.15; tempmax = 23+273.15;
tempvec = Array{Float64}();
if tempmin == tempmax
    tempvec = repeat([tempmin],inner=100);
else
    tempvec = collect(tempmin:((tempmax-tempmin)/(100-1)):tempmax);
end

#How many iterations to save to calculate steady state
savebin=1000;


popstate,
savestate,
toothdrop = popgen(m0,M,tempvec,aalpha,bbeta,maturity,n0,savebin,gen);

R"""
plot($popstate)
"""
