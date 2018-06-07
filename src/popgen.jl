function popgen(m0,M,tempvec,aalpha,bbeta,maturity,n0,savebin,gen)
    
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
    ltime = lspan -1;
    ltemp = length(tempvec);
    tint = Array{Float64}(ltemp,ltime);
    mass = Array{Float64}(ltemp,ltime);
    for k=1:ltemp
        temp = tempvec[k]; #324 to get B0 in natcomm; 291.5 for 65F; 294.3 for 70F
        B0 = (exp(C))/(exp(0.63/((8.61733*10^(-5.0))*temp)));
        a = B0/Em;
        # m0 = 200; #Birth size (grams)
        # M = 500000; #Asymptotic size :: tiger shark: 380000 to 630000 grams
        for i=1:ltime
        	epsilon1 = epsilonvec[i];
        	epsilon2 = epsilonvec[i+1];
        	tint[k,i] = ts(epsilon1,epsilon2,a,eta);
        	mass[k,i] = (epsilon1*M + epsilon2*M)*0.5;
        end
    end
    
    #Cumsum across rows
    tvec = cumsum(tint,2);

    # Survivorship: this is 1 - the probability of death at a given ageclass
    survship = Array{Float64}(ltemp,ltime);
    for k=1:ltemp
        for i=1:ltime
        	survship[k,i] = F(tvec[k,i]);
        end
    end
    survship[:,ltime] = 0;
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
    repsilon = findmin(findmin((maturity-tvec_yrs).^2,1)[1])[2];

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
    temptime = mean(tempvec) + (maximum(tempvec) - mean(tempvec)).*sin.((pi/(daysinyear/2)).*collect(0:1:daysinyear));

    #state vector - number of individuals in given state
    state = zeros(Int64,lspan);
    stateclock = zeros(Float64,lspan);
    state[1] = n0; #start out all individuals at birth size class

    #teeth lost at each state
    toothdrop = zeros(Float64,lspan);
    #mass vector for teeth in each timebin
    bodylength = 5.38674.*mass.^(0.32237); #precaudal length in centimeters
    toothlength = -0.0800064.*(-26.665 - bodylength); #I think also in centimeters? NOTE check
    teethlostperday = 2/40;
    #(teeth/s) :: number of teeth lost per second
    toothlossrate = (teethlostperday/24/60/60); #*toothmass; 
    
    
    # N=100; #how many timebins to save to calculate steady state distribution
    N = savebin;
    savestate = zeros(Int64,N,lspan);
    n=0;

    #Simulation
    tcum = 0;
    tictoc = 0;
    popstate = Array{Int64}(0);
    while tcum < tmax
    	
    	tictoc += 1;
    	pop = sum(state);
    	push!(popstate,pop);
    	
    	if mod(tictoc,100) == 0
    		println("tictoc ",tictoc,"; Population = ",pop)
    	end
    	
    	tcum += tstep;
    	stateclock += tstep;
        
        #What is the temperature?
        day = convert(Int64,round(daysinyear*((tcum/secondsinyear) - floor(tcum/secondsinyear)),0));
        day = maximum([1,day]);
        daytemp = temptime[day];
        #find index for closest match in tempvec
        k = findmin((tempvec-daytemp).^2)[2];
    	
    	#Individuals grow
    	for i=1:ltime
    		
    		#Remove indidivduals that die
    		prmort = 1-survship[k,i];
    		bdist = Binomial(1,prmort);
    		die = sum(rand(bdist,state[i]));
    		#Remove individuals who die
    		state[i] = state[i] - die;
    		
    		#Grow
    		statetime = stateclock[i];
    		# if i < lspan
			if state[i] > 0
				if statetime >= tvec[k,i]
					#Add growing individuals to new state
					state[i+1] = state[i+1]+state[i];
					#Remove growing individuals from current state;
					state[i] = 0;
					#reset clock for that bin
					stateclock[i] = 0;
				end
			end
    		# else 
    		# 	state[i] = 0;
    		# end
    		
    		#Reproduce
    		if i >= repsilon
    			#Recruitment function
    			recruits = convert(Int64,floor((aalpha*state[i])/(1+bbeta*state[i])));
    			#Add the recruits to the initial state
    			state[1] = state[1] + recruits;
    		end
    		
    		#lose teeth
    		# number of ndividuals * tooth loss rate [toothloss/sec] * time interval [sec]
    		toothdrop[i] += state[i]*(toothlossrate * tstep);
    	end
    	
    	#save states from the last N timestep
    	
    	if mod(tictoc,N) == 0
    		n=0;
    	end
    	n+=1;
    	savestate[n,:] = state;
    	
    end


    
    return(
    popstate,
    savestate,
    toothdrop
    )
    
    
    
end
