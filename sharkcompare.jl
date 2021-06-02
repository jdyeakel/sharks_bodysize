if homedir() == "/home/z840"
    loadfunc = include("$(homedir())/sharks_bodysize/src/loadfuncs.jl");
    filename = "SandTiger_all.csv";
    namespace = string("$(homedir())/sharks_bodysize/",filename);
    data = CSV.read(namespace,header=true,DataFrame)
else
    loadfunc = include("$(homedir())/Dropbox/PostDoc/2018_sharks/src/loadfuncs.jl");
    filename = "SandTiger_all.csv";
    namespace = string("$(homedir())/Dropbox/PostDoc/2018_sharks/",filename);
    data = CSV.read(namespace,header=true,DataFrame)
end

# filename = "SandTiger_all.csv";
# namespace = string("$(homedir())/Dropbox/PostDoc/2018_sharks/",filename);
# namespace = string("$(homedir())/sharks_bodysize/",filename);
# data = CSV.read(namespace,header=true,DataFrame)

num = length(names(data));

mchij = SharedArray{Float64}(num,reps,length(sigtauvec),length(tauvec));
mchia = SharedArray{Float64}(num,reps,length(sigtauvec),length(tauvec));
modchij = SharedArray{Float64}(num,reps,length(sigtauvec),length(tauvec));
modchia = SharedArray{Float64}(num,reps,length(sigtauvec),length(tauvec));
moddistchij = SharedArray{Float64}(num,reps,length(sigtauvec),length(tauvec));
moddistchia = SharedArray{Float64}(num,reps,length(sigtauvec),length(tauvec));

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

    for j=1:its
        pos = paramvec[j,:];
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
        toothlength = toothlength1[1,:];
        
        #plot it!
        simdensity_juv = toothdrop[:,1]./sum(toothdrop[:,1]);
        simdensity_adult = toothdrop[:,2]./sum(toothdrop[:,2]);
        U = kde(measures);
        maxy = maximum([simdensity_juv; simdensity_adult; U.density]);
        ply = lineplot(toothlength,simdensity_juv,xlim = [0, 40],ylim=[0,maxy],color=:green)
        lineplot!(ply,toothlength,simdensity_adult,color=:green)
        lineplot!(ply,U.x,U.density,color=:blue)

        meanchij, meanchia, modechij, modechia, modedistchij, modedistchia = empirical_sim_comparison(toothdrop,toothlength,measures)

        mchij[i,r,sigtau_pos,tau_pos] = meachij;
        mchia[i,r,sigtau_pos,tau_pos] = meachia;
        modchij[i,r,sigtau_pos,tau_pos] = modechij;
        modchia[i,r,sigtau_pos,tau_pos] = modechia;
        moddistchij[i,r,sigtau_pos,tau_pos] = modedistchij;
        moddistchia[i,r,sigtau_pos,tau_pos] = modedistchia;
    end

    

end