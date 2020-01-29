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
    
    @everywhere include("$(homedir())/sharks_bodysize/src/popgen_migrate_g.jl")
    @everywhere include("$(homedir())/sharks_bodysize/src/ts.jl")
    @everywhere include("$(homedir())/sharks_bodysize/src/F.jl")
    @everywhere include("$(homedir())/sharks_bodysize/src/smartpath.jl")
    @everywhere include("$(homedir())/sharks_bodysize/src/findlocalmaxima.jl")
    
else
    
    @everywhere include("$(homedir())/Dropbox/PostDoc/2018_sharks/src/popgen_migrate_g.jl")
    @everywhere include("$(homedir())/Dropbox/PostDoc/2018_sharks/src/ts.jl")
    @everywhere include("$(homedir())/Dropbox/PostDoc/2018_sharks/src/F.jl")
    @everywhere include("$(homedir())/Dropbox/PostDoc/2018_sharks/src/smartpath.jl")
    @everywhere include("$(homedir())/Dropbox/PostDoc/2018_sharks/src/findlocalmaxima.jl")

end
