function popgen_migrate_g(m0,M,tempvec1,tempvec2,n0,gen,distance,velocity,D,sigtau,tau)
    
    #Temperature is a vector (!!)
    
    # C = 18.47; #FISH
    # temp = 294.3; #324 to get B0 in natcomm; 291.5 for 65F; 294.3 for 70F
    # B0 = (exp(C))/(exp(0.63/((8.61733*10^(-5.0))*temp)));
    # Em = 5774; #Energy needed to synthesize a unit of mass
    # a = B0/Em;
    # # m0 = 136; #Birth size (grams)
    # # M = 500000; #Asymptotic size :: tiger shark: 380000 to 630000 grams
    eta = 3/4; #Scaling exponent
    C = 18.47; #FISH
    # Em = 5414; #J/g energy density of adult sharks from Schindler
    Em = 5774; #Energy needed to synthesize a unit of mass
    
    epsilonstep = 50;
    epsilonmax = 0.99;
    epsilonvec = collect(m0/M:(epsilonmax - (m0/M))/epsilonstep:epsilonmax);
    lspan = length(epsilonvec);

    #Growth function :: the time it takes to get to epsilon % body mass M from m0
    
    # function ts(epsilon1,epsilon2,a,eta)
    # 	# ts = log(((1-(m0/M)^(1-eta))/(1-epsilon^(1-eta))))*((M^(1-eta)/(a*(1-eta))));
    # 	# Epsilon 1: proportion of M to start
    # 	# Epsilon 2: proportion of M to end
    # 	ts = log(((1-(epsilon1)^(1-eta))/(1-epsilon2^(1-eta))))*((M^(1-eta)/(a*(1-eta))));
    # 	return ts
    # end

    #Gompertz parameters -- somewhat mysterious need to get Calder book
    #from Calder, W. A. Size, function, and life history (Harvard University Press, 1984)

    #Survivorship from Gompertz curve
    # function F(t)
    #     F0 = 1;
    #     a0 = 1.88*10^(-8.0);
    #     b0 = -0.56;
    #     a1 = 1.45*10^(-7.0);
    #     b1 = -0.27;
    #     c0 = a0*M^b0;
    #     c1 = a1*M^b1;
    # 	survship = minimum([1,F0*exp((c0/c1)*(1-exp(c1*t)))]);
    # 	return survship
    # end

    #The time it takes to get from m0 to epsilon M for increasing epsilons.
    #The time experienced by the individual is tvec[3] = tvec[2]-tvec[1]
    ltime = lspan - 1;
    ltemp = length(tempvec1);
    tint1 = Array{Float64}(undef,ltemp,ltime);
    mass1 = Array{Float64}(undef,ltemp,ltime);
    tint2 = Array{Float64}(undef,ltemp,ltime);
    mass2 = Array{Float64}(undef,ltemp,ltime);
    for k=1:ltemp
        temp1 = tempvec1[k]; #324 to get B0 in natcomm; 291.5 for 65F; 294.3 for 70F
        temp2 = tempvec2[k];
        B01 = (exp(C))/(exp(0.63/((8.61733*10^(-5.0))*temp1)));
        B02 = (exp(C))/(exp(0.63/((8.61733*10^(-5.0))*temp2)));
        a1 = B01/Em;
        a2 = B02/Em;
        # m0 = 200; #Birth size (grams)
        # M = 500000; #Asymptotic size :: tigerj shark: 380000 to 630000 grams
        for i=1:ltime
        	epsilon1 = epsilonvec[i];
        	epsilon2 = epsilonvec[i+1];
        	tint1[k,i] = ts(epsilon1,epsilon2,a1,eta);
            tint2[k,i] = ts(epsilon1,epsilon2,a2,eta);
        	mass1[k,i] = (epsilon1*M + epsilon2*M)*0.5;
            mass2[k,i] = (epsilon1*M + epsilon2*M)*0.5;
        end
    end
    
    #Cumsum across rows
    tvec1 = cumsum(tint1,dims=2);
    tvec2 = cumsum(tint2,dims=2);
    # R"plot($(tint1[1,:]))"    
    
    #COMPARE TO VON BUTTERFLY
    # VBlength = VB_Linf*(1 - exp(-VB_k*(VB_t - VB_t0)));
    
    # NOTE: this is currently not a function of size (it is static)
    # Survivorship: this is 1 - the probability of death at a given ageclass
    #mortality rate
    #Sand tiger shark mortality: Z = 0.18 to 0.097 yr-1
    #0.097/(60*60*24*365) = 3.07585*10^-9
    
    # m=(3.076*10^(-9.0))*1;
    m = (5.70776*10^(-9));
    
    #CURRENTLY NOT USING THIS
    # survship1 = Array{Float64}(undef,ltemp,ltime);
    # survship2 = Array{Float64}(undef,ltemp,ltime);
    # for k=1:ltemp
    #     for i=1:ltime
    #     	survship1[k,i] = F(tvec1[k,i],m);
    #         survship2[k,i] = F(tvec2[k,i],m);
    #     end
    # end
    # survship1[:,ltime] .= 0;
    # survship2[:,ltime] .= 0;
    
    
    # Reproduction
    
    #Fecundities
    
    #Litters every 2 years for blue sharks
    # littersize / 365*2 days
    # littersize_twoyrs = -3.349 .+ (0.179 .* (49.9394 .* ((mass ./ 1000) .^ 0.3223726627981947))); 
    # littersize_sec = littersize_twoyrs ./ (60*60*24*365*2);
    # offspring_persize = littersize_sec .* tint;
    
    #Offsping weight
    
    #Branstetter and Musick (1994) -- 1 litter per 2 years and 2 pups per litter = 3.17*10^-8
    #Female pups per female per second 0.5/(60*60*24*365) = 1.58549*10^-8
    #Number of female pups per female per second = 1.585*10^-8
    #NOTE: smaller size bins have 'more' offsprings because bins are larger
    
    #Branstetter and Musick (1994) -- 1 litter per 2 years and 2 pups per litter 
    # (1/2)/(60*60*24*365) = 3.17*10^-8
    # rmax_warm = (1.585*10^-8.0);
    
    #Cortes 1996 Bonnetthead shark 4.45 and 4.65 female pups/female
    #  4.65/(60*60*24*365)
    rmax_warm = (1.47451*10^(-7));
    
    #NOTE: could make this percent difference a function of Delta temp
    rmax_cold = (rmax_warm)*1;
    
    juvpos = Int64(floor(ltime/4)); #What defines 'juvenile state'
    
    #reproduction rate per state/location
    mintemp1 = minimum(tempvec1);
    maxtemp1 = maximum(tempvec1);
    mintemp2 = minimum(tempvec2);
    maxtemp2 = maximum(tempvec2);
    r_sizetemp1 = Array{Float64}(undef,ltemp,ltime);
    r_sizetemp2 = Array{Float64}(undef,ltemp,ltime);
    # offspring_sizetemp1 = Array{Float64}(undef,ltemp,ltime);
    # offspring_sizetemp2 = Array{Float64}(undef,ltemp,ltime);
    for i=1:ltemp
        for j=1:ltime
            
            #interpolated reproductive rate for a given temperature
            rtemp1 = rmax_cold + (rmax_warm - rmax_cold)*((tempvec1[i]-mintemp1)/(maxtemp1-mintemp1));
            
            #reproductive rate as a function of mass and temperature
            # r_sizetemp1[i,j] = (-4*(mass1[i,j] - mass1[i,ltime])*(mass1[i,j] - m0)*rtemp1)/((mass1[i,ltime] - m0)^2);
            
            #ALT sigmoid relationship
            # r_sizetemp1[i,j] = rtemp1 ./ (1 .+ exp.(- (1/1) .* (mass1[i,j] .- juvpos)));
            
            #ALT single value
            r_sizetemp1[i,j] = rtemp1;
            
            #offspring in time interval tint
            #offspring_sizetemp1[i,j] = r_sizetemp1[i,j]*tint1[i,j];
            
            rtemp2 = rmax_cold + (rmax_warm - rmax_cold)*((tempvec2[i]-mintemp2)/(maxtemp2-mintemp2));
            
            # r_sizetemp2[i,j] = (-4*(mass2[i,j] - mass2[i,ltime])*(mass2[i,j] - m0)*rtemp2)/((mass2[i,ltime] - m0)^2);
            
            #ALT sigmoid relationship
            # r_sizetemp2[i,j] = rtemp2 ./ (1 .+ exp.(- (1/1) .* (mass2[i,j] .- juvpos)));
            
            #ALT single value
            r_sizetemp2[i,j] = rtemp2;
            
            #offspring_sizetemp2[i,j] = r_sizetemp2[i,j]*tint2[i,j];
        end
    end
    
    # offspring_persize = r .* tint;
    
    #assume maturity at 6 years old #Otway et al. 2004
    #This gives the size bin at which reproduction is initiated for different temperature environments
    # sizeatmaturity = Array{Int64}(undef,ltemp)
    # for i=1:ltemp
    #     sizeatmaturity[i] = findall(x->x >= 2.0,tvec[i,:]/(60*60*24*365))[1];
    # end
    # 
    # #Assume minimum size at maturity
    # minmat = minimum(sizeatmaturity);
    # 
    # #Births per individual in each size class (starting at 2:100)
    # birthsperind = diff(floor.(cumsum(offspring_persize,dims=2)),dims=2);
    # 
    #Given offspring_persize, the rate of offspring production r is the same. The only thing that changes is the size of the time interval r*tint
    #For determining recruitment:
    #1) determine number of individuals > age at maturity B
    #2) R = Int64(floor(r*B))
    #3) N = N + R
    
    #If it is mass specific growth, just do this for each mass class
    
    # # Reproductive Output: 0.0066*mass kg per kg day (Schindler)
    # repout_gday = 0.0066 .* (mass); #g per g day
    # repout_ggs = repout_gday ./ (60*60*24); #grams per gram seconds
    # #length of time spent in each size class per gram
    # transitionpermass = (tint) .* mass; #grams/second
    # repout_g = (repout_ggs .* transitionpermass); #per gram
    
    # NOTE: Leslie matrix version
    # transitionprobs = survship[1,2:ltime];
    # Lmatrix = zeros(Float64,ltime,ltime);
    # Lmatrix[diagind(Lmatrix,-1)] = transitionprobs;
    # Lmatrix[1,sizeatmaturity:size(Lmatrix)[1]] = offspring_persize[1,sizeatmaturity:size(offspring_persize)[2]];
    
    #Check relationships

    #=====================#
    #Simulation
    #=====================#
    # tvec_yrs = tvec/60/60/24/365;
    #Reproduction parameters
    # alpha = 1;
    # beta = 0.1;
    # maturity = 10; #in years
    # repsilon = findmin((maturity-tvec_yrs).^2)[2]
    # repsilon = findmin(vec(findmin((maturity.-tvec_yrs).^2,dims=1)[1]))[2];

    # R"""
    # par(mfrow=c(2,1))
    # plot($(tvec_yrs),$(mass),xlab='Growth time',ylab='Mass')
    # lines(rep($(tvec_yrs[repsilon]),length.out=length(seq(1,10^6,by=1000))),seq(1,10^6,by=1000))
    # """



    # gen = 20;
    # tstep = minimum(tint); #The timesteps are set at the minimum time interval between age classes and temps
    tmax = mean([maximum(tvec1),maximum(tvec2)])*gen; #How much time to simulate?
    # totsteps = tmax/tstep;
    #How many individuals start
    # n0 = 1000;
    
    #Keep track of temperature changes by days
    daysinyear = 365;
    secondsinday = 24*60*60;
    secondsinyear = daysinyear*secondsinday;
    
    #temporal temperature vectors for sites 1 and 2
    temptime1 = mean(tempvec1) .+ (maximum(tempvec1) .- mean(tempvec1)).*sin.((pi/(daysinyear/2)).*collect(0:1:daysinyear));
    temptime2 = mean(tempvec2) .+ (maximum(tempvec2) .- mean(tempvec2)).*sin.((pi/(daysinyear/2)).*collect(0:1:daysinyear));
    
    
    
    
    
    #teeth lost at each state
    
    #mass vector for teeth in each timebin
    #NOTE: CHANGE TO THE MEASURE FROM Branstetter and MUSICK
    # bodylength = 5.38674.*mass.^(0.32237); #precaudal length in centimeters
    
    #Assume sharks lose 1 upper and 1 lower tooth per 40 days
    
    #From Shimata (A1 position)
    #Upper A1
    #Shimata 2004
    #bodylength (cm) = -26.665+12.499*crownheight (mm)
    # crownheight = 
    # toothlength = -0.0800064.*(-26.665 .- bodylength); #I think also in centimeters? NOTE check
    upperteethlostperday = 1/40;
    #Lower A1
    #Shimata 2004
    #bodylength (cm) = -24.722+10.305*crownheight (mm)
    #(teeth/s) :: number of teeth lost per second
    lowerteethlostperday = 1/40;
    
    toothlossrate = ((lowerteethlostperday + upperteethlostperday)/24/60/60); #*toothmass; 
    
    
    # N=100; #how many timebins to save to calculate steady state distribution
    # savestate = zeros(Int64,savebin,ltime,2);
    

    #state vector - number of individuals in given state
    state = zeros(Int64,ltime,2);
    # stateclock = zeros(Float64,ltime,2);
    state[1,1] = n0; #Nursury has n0 newborns
    state[juvpos,2] = n0; #Adult site has n0 juveniles
    
    toothdrop = zeros(Float64,ltime,2);
    
    pop1 = Array{Int64}(undef,0);
    pop2 = Array{Int64}(undef,0);
    # popstate[1,:] = sum(state,dims=1);
    
    clock = Array{Float64}(undef,0);
    
    staterate = [state;state;state;state];
    cistaterate = CartesianIndices(staterate);
    states = repeat(collect(1:ltime),outer=4);
    
    let n=0, tcum = 0, tictoc = 0, day = 1 #n=0; tcum = 0; tictoc = 0; day = 1
        
        while tcum < tmax #&& tictoc < 1000000
            
            tictoc += 1;
        	
            daytemp1 = temptime1[day];
            daytemp2 = temptime2[day];
            
            #Define rates for both sites 1 & 2
            
            #find index for closest match in tempvec
            k1 = findmin((tempvec1 .- daytemp1).^2)[2];
            k2 = findmin((tempvec2 .- daytemp2).^2)[2];
            
            #1) Calculate Rate
            N = sum(state,dims=1);
            N1 = state[:,1];
            N2 = state[:,2];
            
            if sum(state) == 0
                break
            end
            
            #reproduction rate per state/location
            #A nonlinear function, dependent on temp and size
            r1 = r_sizetemp1[k1,:];
            
            #A single value
            # r1 = repeat([rmax_warm],ltime);
            
            #turn off pre-juv reproduction
            r1[1:juvpos] .= 0;
            
            #No reproduction at adult site
            r2 = repeat([0],ltime);
    
            #growth rate per state/location
            g1 = 1 ./ tint1[k1,:]; #transpose(1/tint1[k1,:]);
            g2 = 1 ./ tint2[k2,:]; #transpose(1/tint2[k2,:]);
            
            #mortality rate per state/location
            d2 = repeat([m],ltime);
            #elevated mortality for last mass class
            d2[ltime] = m*10;
            d1 = copy(d2);
            
            #migration rate per state/location
            #peak travel day
            peakday = 180;
            
            # m1 = zeros(Float64,ltime);
            # m1[1:Int64(floor(juvpos/2))] .= 0;
            # m1[juvpos:ltime] .= D*velocity/distance;
            
            #For Juvenile site, migration is a function of size
            #sigmoidal increase over body mass to top migration above juvpos
            # sigtau = 5;
            m1 = D*(velocity/distance) ./ (1 .+ exp.(- (1/sigtau) .* (collect(1:ltime) .- juvpos)));

            # The corresponding masses along the sigtauvec are:
            # sigtauvec*(M/epsilonstep)

            # The last bit could be written as: (in terms of masses)
            # ((mass1[1,:] .- mass1[1,juvpos]) ./ maximum(mass1[1,:])) .* length(mass1[1,:])
            
            #For adult site, migration is a function of TIME
            # tau = 20;
            m2 = repeat([D*(velocity/distance)*exp(-((day .- peakday)^2)/(2*tau^2))],ltime);
            # m2 = zeros(Float64,ltime);
            # m2[1:juvpos] .= (velocity/2)/distance;
            # m2[juvpos:ltime] .= velocity/distance;
            
            Rate = dot(N1,r1 .+ g1 .+ d1 .+ m1) + dot(N2,r2 .+ g2 .+ d2 .+ m2);
            
            dt = 1/Rate;
            
            #Alternative NAline
            
            # Reproduce // Grow // Die // Move # Reproduce // Grow // Die // Move
            NAline_staterates = [state[:,1] .* r1; state[:,1] .* g1; state[:,1] .* d1; state[:,1] .* m1; state[:,2] .* r2; state[:,2] .* g2; state[:,2] .* d2; state[:,2] .* m2];
            
            NAline = cumsum(NAline_staterates/sum(NAline_staterates));
            
            NArand = rand();
            NApos = findall(isodd,NArand .<= NAline)[1];
            
            Nloc = cistaterate[NApos][2];
            Nstate = states[cistaterate[NApos][1]];
            Apos = findall(isodd,cistaterate[NApos][1] .<= (ltime .* [1,2,3,4]))[1];
            
            Altloc = setdiff([1,2],Nloc)[1];
            # 
            # 
            # #grab random individual from either site
            # Nline = cumsum(vec(reshape(state,length(state),1))/sum(state));
            # Nrand = rand();
            # Npos = findall(isodd,Nrand .<= Nline)[1];
            # 
            # cistate = CartesianIndices(state);
            # Nstate = cistate[Npos][1];
            # Nloc = cistate[Npos][2];
            # 
            # Altloc = setdiff([1,2],Nloc)[1];
            # 
            # #Now choose the activity to update
            # # Reproduce // Grow // Die // Move
            # 
            # if Nloc == 1
            #     #At site 1
            #     Aline = cumsum([r1[Nstate], g1[Nstate], d1[Nstate], m1[Nstate]]/sum([r1[Nstate], g1[Nstate], d1[Nstate], m1[Nstate]]));
            # else
            #     #At site 2
            #     Aline = cumsum([r2[Nstate], g2[Nstate], d2[Nstate], m2[Nstate]]/sum([r2[Nstate], g2[Nstate], d2[Nstate], m2[Nstate]]));
            # end
            # 
            # Arand = rand();
            # Apos = findall(isodd,Arand .<= Aline)[1];
            # 
            #Reproduce
            if Apos == 1
                #Grab a birth state at random
                rbirthstate = rand(collect(1:Int64(epsilonstep/10)));
                # rbirthstate = 1;
                #Add individual to first state
                state[rbirthstate,Nloc] += 1;
            end
            
            #Grow 
            if Apos == 2
                #Advance individual to next state
                state[Nstate,Nloc] -= 1;
                if Nstate + 1 < ltime
                    state[Nstate+1,Nloc] += 1;
                else
                    state[ltime,Nloc] += 1;
                end
            end
            
            #Die
            if Apos == 3
                #subtract individual from state
                state[Nstate,Nloc] -= 1;
                #drop teeth
                # toothdrop[Nstate,Nloc] += 8;
            end
            
            #Move (and grow while moving if distance is far enough)
            if Apos == 4
                #If there is growth to a new state over migration, advance state and location
                #Otherwise advance location
                if Nloc == 1
                    growthovermigration = 0; #(distance/velocity)/tint1[k1,Nstate];
                else
                    growthovermigration = 0; #(distance/velocity)/tint2[k2,Nstate];
                end
                if growthovermigration > 1 && Nstate < ltime
                    dstate = Int64(floor(growthovermigration));
                    state[Nstate,Nloc] -= 1;
                    if (Nstate + dstate) < ltime 
                        state[Nstate+dstate,Altloc] += 1;
                    else
                        state[ltime,Altloc] += 1;
                    end
                else
                    state[Nstate,Nloc] -= 1;
                    state[Nstate,Altloc] += 1;
                end
            end
            
            #Update time
            tcum += dt;
        	# stateclock .+= dt;
            push!(clock,tcum);
            
            #Advance temperature in both sites with time
            day = round(Int64,daysinyear*((tcum/secondsinyear) - floor(tcum/secondsinyear)));
            day = maximum([1,day]);
    
            #Update population
            pop = vec(sum(state,dims=1));
            push!(pop1,pop[1]);
            push!(pop2,pop[2]);
        	# push!(popstate,pop);
            
            #Tooth drop
            toothdrop .+= state .* (toothlossrate * dt);
            
            #Report
        	# if mod(tictoc,5000) == 0
        	# 	println("tcum ",tcum/tmax,"; Juv site = ",pop[1]," Adult site = ",pop[2])
        	# end
        	
            #test for negative values
            if any(state .< 0) == true
                println("RULES OF UNIVERSE VIOLATED")
                break
            end
        
        end #end while
    end #end let
    
    popstate = [pop1 pop2];
    
    
    
    # R"""
    # plot($clock/60/60/24/365,$(popstate[:,1]),pch='.',ylim=c(0,max($popstate)))
    # points($clock/60/60/24/365,$(popstate[:,2]),pch='.',col='blue')
    # """
    
    return(
    mass1,
    mass2,
    epsilonvec,
    clock,
    popstate,
    toothdrop,
    state
    )
    
    
    
end
