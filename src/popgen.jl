function popgen(m0,M,tempvec,n0,savebin,gen,poprep)
    
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
    
    epsilonstep = 100;
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
    ltemp = length(tempvec);
    tint = Array{Float64}(undef,ltemp,ltime);
    mass = Array{Float64}(undef,ltemp,ltime);
    for k=1:ltemp
        temp = tempvec[k]; #324 to get B0 in natcomm; 291.5 for 65F; 294.3 for 70F
        B0 = (exp(C))/(exp(0.63/((8.61733*10^(-5.0))*temp)));
        a = B0/Em;
        # m0 = 200; #Birth size (grams)
        # M = 500000; #Asymptotic size :: tigerj shark: 380000 to 630000 grams
        for i=1:ltime
        	epsilon1 = epsilonvec[i];
        	epsilon2 = epsilonvec[i+1];
        	tint[k,i] = ts(epsilon1,epsilon2,a,eta);
        	mass[k,i] = (epsilon1*M + epsilon2*M)*0.5;
        end
    end
    
    #Cumsum across rows
    tvec = cumsum(tint,dims=2);
    
    #COMPARE TO VON BUTTERFLY
    # VBlength = VB_Linf*(1 - exp(-VB_k*(VB_t - VB_t0)));
    
    # NOTE: this is currently not a function of size (it is static)
    # Survivorship: this is 1 - the probability of death at a given ageclass
    #mortality rate
    m=(3.076*10^(-9.0))*0.09;
    survship = Array{Float64}(undef,ltemp,ltime);
    for k=1:ltemp
        for i=1:ltime
        	survship[k,i] = F(tvec[k,i],m);
        end
    end
    survship[:,ltime] .= 0;
    
    
    # Reproduction
    
    #Fecundities
    
    #Litters every 2 years for blue sharks
    # littersize / 365*2 days
    # littersize_twoyrs = -3.349 .+ (0.179 .* (49.9394 .* ((mass ./ 1000) .^ 0.3223726627981947))); 
    # littersize_sec = littersize_twoyrs ./ (60*60*24*365*2);
    # offspring_persize = littersize_sec .* tint;
    
    #Offsping weight
    
    #Branstetter and Musick (1994) -- 1 litter per 2 years and 2 pups per litter = 3.17*10^-8
    
    #Number of female pups per female per second = 1.585*10^-8
    #NOTE: smaller size bins have 'more' offsprings because bins are larger
    rmax_warm = (1.585*10^-8.0);
    #NOTE: could make this percent difference a function of Delta temp
    rmax_cold = (1.585*10^-8.0)*0.5;
    
    mintemp = minimum(tempvec);
    maxtemp = maximum(tempvec);
    r_sizetemp = Array{Float64}(undef,ltemp,ltime);
    offspring_sizetemp = Array{Float64}(undef,ltemp,ltime);
    for i=1:ltemp
        for j=1:ltime
            rtemp = rmax_cold + (rmax_warm - rmax_cold)*((tempvec[i]-mintemp)/(maxtemp-mintemp));
            r_sizetemp[i,j] = (-4*(mass[i,j] - mass[i,ltime])*(mass[i,j] - m0)*rtemp)/((mass[i,ltime] - m0)^2);
            offspring_sizetemp[i,j] = r_sizetemp[i,j]*tint[i,j];
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
    tvec_yrs = tvec/60/60/24/365;
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
    tstep = minimum(tint); #The timesteps are set at the minimum time interval between age classes and temps
    tmax = maximum(tvec)*gen; #How much time to simulate?
    totsteps = tmax/tstep;
    #How many individuals start
    # n0 = 1000;
    
    #Keep track of temperature changes by days
    daysinyear = 365;
    secondsinday = 24*60*60;
    secondsinyear = daysinyear*secondsinday;
    temptime = mean(tempvec) .+ (maximum(tempvec) .- mean(tempvec)).*sin.((pi/(daysinyear/2)).*collect(0:1:daysinyear));

    #state vector - number of individuals in given state
    state = zeros(Int64,ltime);
    stateclock = zeros(Float64,ltime);
    state[1] = n0; #start out all individuals at birth size class

    #teeth lost at each state
    toothdrop = zeros(Float64,ltime);
    #mass vector for teeth in each timebin
    #NOTE: CHANGE TO THE MEASURE FROM Branstetter and MUSICK
    bodylength = 5.38674.*mass.^(0.32237); #precaudal length in centimeters
    
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
    savestate = zeros(Int64,savebin,ltime);
    
    
    juvpos = Int64(floor(ltime/4)); #What defines 'juvenile state'
    sharpness = 0.25; #high is sharp
    if poprep == "adult"
        probdrop = 1 ./ (1 .+ exp.(-sharpness .* (collect(1:ltime) .- juvpos)));
    end
    if poprep == "juv"
        probdrop = 1 .- (1 ./ (1 .+ exp.(-sharpness .* (collect(1:ltime) .- juvpos))));
    end
    if poprep == "both"
        probdrop = repeat([1.0],outer=ltime);
    end
    

    #Simulation
    popstate = Array{Int64}(undef,0);
    clock = Array{Float64}(undef,0);
    let n=0, tcum = 0, tictoc = 0
        while tcum < tmax
        	
        	tictoc += 1;
        	pop = sum(state);
        	push!(popstate,pop);
            
        	
        	if mod(tictoc,100) == 0
        		println("tictoc ",tictoc,"; Population = ",pop)
        	end
        	
        	tcum += tstep;
        	stateclock .+= tstep;
            push!(clock,tcum);
            
            #What is the temperature?
            day = round(Int64,daysinyear*((tcum/secondsinyear) - floor(tcum/secondsinyear)));
            day = maximum([1,day]);
            daytemp = temptime[day];
            #find index for closest match in tempvec
            k = findmin((tempvec .- daytemp).^2)[2];
            
            
            #drop teeth
            for i=1:ltime
                rdraw = rand(state[i]);
                toothdrop[i] += sum(rdraw .< probdrop[i]) * (toothlossrate * tstep);
            end
            
            newstate = zeros(Int64,ltime);
        	#Individuals grow
        	for i=1:(ltime-1)
        		#Grow
        		statetime = stateclock[i];
        		# if i < lspan
    			if state[i] > 0
                    #GROWTH
    				if statetime >= tint[k,i]
    					#Add growing individuals to new state
    					newstate[i+1] = state[i];
    					#reset clock for that bin
    					stateclock[i] = 0;
                    #NO GROWTH
    				else
                        newstate[i] = state[i];
                    end
    			end
        	end
            #save updated state as current state 
            state = copy(newstate);
            
            #Reproduction
            recruitmass = sum(r_sizetemp[k,:] .* tstep .* state);
            recruits = Int64(floor(recruitmass));
            state[1] += recruits;
            
            #Mortality
            #Remove indidivduals that die
            # prmort = (1 .- (survship[k,:]));
            #sets carrying capacity
            K = 100;
            prmort = 1 - exp(-tstep*m*(sum(state)/K));
            bdist = Binomial(1,prmort);
            die = zeros(Int64,ltime);
            for i=1:(ltime-1)
                if state[i] > 0
                    die[i] = sum(rand(bdist,state[i]));
                end
            end
            #all at last age-class die
            # die[ltime] = state[ltime];
            
            #Remove individuals who die
            state .-= die;
            
            #death tooth drop
            # toothdrop .+= die .* ((toothlossrate*8) * tstep);
            
            for i=1:ltime
                if die[i] > 0
                    rdraw = rand(die[i]);
                    toothdrop[i] += sum(rdraw .< probdrop[i]) * ((toothlossrate*8) * tstep);
                end
            end
            
        	
        	#save states from the last savebin timestep
        	if mod(tictoc,savebin) == 0
        		n=0;
        	end
        	n += 1;
        	savestate[n,:] = state;
        	
        end #end while
        
    end
    
    # R"plot($(popstate))"

    
    return(
    mass,
    epsilonvec,
    clock,
    popstate,
    savestate,
    toothdrop,
    state
    )
    
    
    
end
