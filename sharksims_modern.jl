if homedir() == "/home/z840"
    loadfunc = include("$(homedir())/sharks_bodysize/src/loadfuncs.jl");
else
    loadfunc = include("$(homedir())/Dropbox/PostDoc/2018_sharks/src/loadfuncs.jl");
end



#Size at birth (cm)
l0 = 90.0;
#Asymptotic size (cm)
# L = 295; #1500; #295.0;
# NOTE: L = 477 cm (15.65 ft) -> Mass 350 KG -> max tooth length 40 mm

L = 300;
# L = 477;

# # Mass from PRECAUDAL LENGTH (Schindler) 
# m0 = 0.00538776*l0^3.102025622731644;
# M = 0.00538776*L^3.102025622731644;

# Mass (kg) from TOTAL LENGTH (Goldman et al. 2006)
# NOTE: Mass of 350 KG (350000g) results in a tooth length of 40 mm (max emp. length)
# NOTE: Mass of 350 KG results in L = 477 cm = 15.65 ft
m0 = (0.00013*l0^2.4)*1000;
M = (0.00013*L^2.4)*1000;

#Sim params
n0=1000;
gen=1;


mintemp_j = 17;
maxtemp_j = 25;
mintemp_a = 13;
maxtemp_a = 23;

distvec = 700;
sigtauvec = collect(0.5:3:25); # collect(0.5:0.5:25)
tauvec = collect(1:5:50); #collect(1:1:50);
tauits = length(sigtauvec)*length(tauvec);


#3 paramters: distance, sigtau, tau
# distposvec = repeat(collect(1:3),inner = length(sigtauvec)*length(tauvec));
sigtauposvec = repeat(collect(1:length(sigtauvec)),inner = length(tauvec));
tauposvec = repeat(collect(1:length(tauvec)),outer=length(sigtauvec));
paramposvec_pre = [sigtauposvec tauposvec];
#temperature | distance | sigtauvec | tauposvec
# paramposvec = [repeat(collect(1:4),inner = size(paramposvec_pre)[1]) repeat(paramposvec_pre,outer=4)];

reps = 5;
paramvec_pre = repeat(paramposvec_pre,outer = reps);
paramvec = [repeat(collect(1:reps),inner=size(paramposvec_pre)[1]) paramvec_pre];

its = size(paramvec)[1];

@time @sync @distributed for i=1:its
    
    #set parameters
    pos = paramvec[i,:];
    r = pos[1];
    # temp_pos = pos[2];
    # dist_pos = pos[3];
    sigtau_pos = pos[2];
    tau_pos = pos[3];
    
    #Temperature range (high latitude: 8-13; Low latitude 22-30)
    #Juvenile site
    tempmin1 = mintemp_j+273.15; tempmax1 = maxtemp_j+273.15;
    #Adult site
    tempmin2 = mintemp_a+273.15; tempmax2 = maxtemp_a+273.15;

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
    distance = distvec*1000; #3.779e6/10; #1500
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
    filename = "data/sharks_modern/simdata.jld";
    indices = [r,sigtau_pos,tau_pos];
    namespace = smartpath(filename,indices);
    
    # @save namespace mass1 mass2 epsilonvec clock popstate toothdrop state toothlength1 toothlength2;
    
    @save namespace mass1 mass2 toothdrop toothlength1 toothlength2;
    
end

