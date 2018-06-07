@everywhere using RCall
@everywhere include("$(homedir())/Dropbox/PostDoc/2018_sharks/src/popgen.jl")



m0 = 136;
M = 500000;
#the range of temperature experienced
tempmin = 280; tempmax = 300;
tempvec = collect(tempmin:((tempmax-tempmin)/(100-1)):tempmax);
