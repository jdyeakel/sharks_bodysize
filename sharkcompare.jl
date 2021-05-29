if homedir() == "/home/z840"
    loadfunc = include("$(homedir())/sharks_bodysize/src/loadfuncs.jl");
else
    loadfunc = include("$(homedir())/Dropbox/PostDoc/2018_sharks/src/loadfuncs.jl");
end

filename = "SandTiger_all.csv";
namespace = string("$(homedir())/Dropbox/PostDoc/2018_sharks/",filename);
data = CSV.read(namespace,header=true,DataFrame)

num = length(names(data));

for i=1:num

    measures = Array{Float64}(data[!,i][findall(!ismissing,data[!,i])]);


    #Now walk through the simulation results

    sigtauvec = collect(0.5:0.5:25); # collect(0.5:0.5:25)
    tauvec = collect(1:1:50); #collect(1:1:50);
    tauits = length(sigtauvec)*length(tauvec);


    #3 paramters: distance, sigtau, tau
    # distposvec = repeat(collect(1:3),inner = length(sigtauvec)*length(tauvec));
    sigtauposvec = repeat(collect(1:length(sigtauvec)),inner = length(tauvec));
    tauposvec = repeat(collect(1:length(tauvec)),outer=length(sigtauvec));
    paramposvec_pre = [sigtauposvec tauposvec];
    #temperature | distance | sigtauvec | tauposvec
    # paramposvec = [repeat(collect(1:4),inner = size(paramposvec_pre)[1]) repeat(paramposvec_pre,outer=4)];

    reps = 25;
    paramvec_pre = repeat(paramposvec_pre,outer = reps);
    paramvec = [repeat(collect(1:reps),inner=size(paramposvec_pre)[1]) paramvec_pre];

    its = size(paramvec)[1];

    for i=1:its
        pos = paramvec[i,:];
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
        

        empirical_sim_comparison(toothdrop,toothlength1,measures)

