if homedir() == "/home/z840"
    loadfunc = include("$(homedir())/sharks_bodysize/src/loadfuncs.jl");
else
    loadfunc = include("$(homedir())/Dropbox/PostDoc/2018_sharks/src/loadfuncs.jl");
end



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


mintemp_j = [8,20,9,12];
maxtemp_j = [13,30,17,24];
mintemp_a = [8,20,12,9];
maxtemp_a = [13,30,24,17];

distvec = [200,400,650];
sigtauvec = collect(1:10);
tauvec = collect(5:25);
tauits = length(sigtauvec)*length(tauvec)*length(distvec);


#3 paramters: distance, sigtau, tau
distposvec = repeat(collect(1:3),inner = length(sigtauvec)*length(tauvec));
sigtauposvec = repeat(collect(1:length(sigtauvec)),inner = length(tauvec),outer=length(distvec));
tauposvec = repeat(collect(1:length(tauvec)),outer=length(distvec)*length(sigtauvec));
paramposvec_pre = [distposvec sigtauposvec tauposvec];
#temperature | distance | sigtauvec | tauposvec
paramposvec = [repeat(collect(1:4),inner = size(paramposvec_pre)[1]) repeat(paramposvec_pre,outer=4)];



its = size(paramposvec)[1];

@time @sync @distributed for i=1:its
    
    #set parameters
    pos = paramposvec[i,:];
    temp_pos = pos[1];
    dist_pos = pos[2];
    sigtau_pos = pos[3];
    tau_pos = pos[4];
    
    #Temperature range (high latitude: 8-13; Low latitude 22-30)
    #Juvenile site
    tempmin1 = mintemp_j[temp_pos]+273.15; tempmax1 = maxtemp_j[temp_pos]+273.15;
    #Adult site
    tempmin2 = mintemp_a[temp_pos]+273.15; tempmax2 = maxtemp_j[temp_pos]+273.15;

    tempvec1 = Array{Float64}(undef,0);
    tempvec2 = Array{Float64}(undef,0);
    if tempmin1 == tempmax1
        tempvec1 = repeat([tempmin1],inner=100);
    else
        tempvec1 = collect(tempmin1:((tempmax1-tempmin1)/(100-1)):tempmax1);
    end;
    if tempmin2 == tempmax2
        tempvec2 = repeat([tempmin2],inner=100);
    else
        tempvec2 = collect(tempmin2:((tempmax2-tempmin2)/(100-1)):tempmax2);
    end;

    #distance between site (km * 1000)
    distance = distvec[dist_pos]*1000; #3.779e6/10; #1500
    #Shark velocity (m/s)
    velocity = 1;
    D = 1;

    #Steepness of juvenile migration (function of mass)
    sigtau = sigtauvec[sigtau_pos];
    #Steepness of adult migration (function of time)
    tau = tauvec[tau_pos];

    #run sim
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
    
    #save data
    filename = "data/sharks_3oceans/simdata.jld";
    indices = [temp_pos,dist_pos,sigtau_pos,tau_pos];
    namespace = smartpath(filename,indices);
    
    @save namespace mass1 mass2 epsilonvec clock popstate toothdrop state toothlength1 toothlength2;
    
end

