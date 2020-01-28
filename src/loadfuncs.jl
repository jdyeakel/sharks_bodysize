using Distributed

@everywhere using Distributions
@everywhere using RCall
# @everywhere using HDF5
@everywhere using JLD2
@everywhere using DataFrames
@everywhere using LinearAlgebra
@everywhere using SharedArrays
# @everywhere using EcologicalNetworks



if homedir() == "/home/z840"
    
    @everywhere include("$(homedir())/2018_sharks/src/popgen_migrate_g.jl")
    @everywhere include("$(homedir())/2018_sharks/src/ts.jl")
    @everywhere include("$(homedir())/2018_sharks/src/F.jl")
    @everywhere include("$(homedir())/2018_sharks/src/smartpath.jl")
    
else
    
    @everywhere include("$(homedir())/Dropbox/PostDoc/2018_sharks/src/popgen_migrate_g.jl")
    @everywhere include("$(homedir())/Dropbox/PostDoc/2018_sharks/src/ts.jl")
    @everywhere include("$(homedir())/Dropbox/PostDoc/2018_sharks/src/F.jl")
    @everywhere include("$(homedir())/Dropbox/PostDoc/2018_sharks/src/smartpath.jl")
    
end
