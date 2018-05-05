using Distributions
using HDF5
using JLD
using RCall


#Initialize variables
n0 = 100;

agestep = 100;
agevec = collect(0:1/agestep:1);


a = 1;
m0 = 1;
M = 1;
eta = 3/4;

function ts(espilon)
	ts = ln((1-(m0/M)^(1-eta)/(1-epsilon^(1-eta))))*((M^(1-eta)/(a*(1-eta))));
	return ts
end

tvec = mapslices(ts,agevec);



