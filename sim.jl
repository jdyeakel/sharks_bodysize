using Distributions
using JLD2
using RCall

#=====================#
#Setup
#=====================#

#Initialize variables


#Gillooly relationships (Nature, 2002)
# B0 = 0.047; #Metabolic normalization constant
# C = 19.50; #endotherms (Brown Ecology 2004)
C = 18.47; #FISH
temp = 324; #324 to get B0 in natcomm; 291.5 for 65F; 294.3 for 70F
B0 = (exp(C))/(exp(0.63/((8.61733*10^(-5.0))*temp)));
Em = 5774; #Energy needed to synthesize a unit of mass
a = B0/Em;
m0 = 200; #Birth size (grams)
M = 500000; #Asymptotic size :: tiger shark: 380000 to 630000 grams
eta = 3/4; #Scaling exponent

epsilonstep = 100;
epsilonmax = 0.99;
epsilonvec = collect(m0/M:(epsilonmax - (m0/M))/epsilonstep:epsilonmax);
lspan = length(epsilonvec);

#Growth function :: the time it takes to get to epsilon % body mass M from m0
function ts(epsilon1,epsilon2)
	# ts = log(((1-(m0/M)^(1-eta))/(1-epsilon^(1-eta))))*((M^(1-eta)/(a*(1-eta))));
	# Epsilon 1: proportion of M to start
	# Epsilon 2: proportion of M to end
	ts = log(((1-(epsilon1)^(1-eta))/(1-epsilon2^(1-eta))))*((M^(1-eta)/(a*(1-eta))));
	return ts
end

#Gompertz parameters -- somewhat mysterious need to get Calder book
#from Calder, W. A. Size, function, and life history (Harvard University Press, 1984)
F0 = 1;
a0 = 1.88*10^(-8.0);
b0 = -0.56;
a1 = 1.45*10^(-7.0);
b1 = -0.27;
c0 = a0*M^b0;
c1 = a1*M^b1;
#Survivorship from Gompertz curve
function F(t)
	survship = minimum([1,F0*exp((c0/c1)*(1-exp(c1*t)))]);
	return survship
end

#The time it takes to get from m0 to epsilon M for increasing epsilons.
#The time experienced by the individual is tvec[3] = tvec[2]-tvec[1]
tint = Array{Float64}(undef,lspan-1);
mass = Array{Float64}(undef,lspan-1);
for i=1:lspan-1
	epsilon1 = epsilonvec[i];
	epsilon2 = epsilonvec[i+1];
	tint[i] = ts(epsilon1,epsilon2);
	mass[i] = (epsilon1*M + epsilon2*M)*0.5;
end

tvec = cumsum(tint);

# Survivorship: this is 1 - the probability of death at a given ageclass
ltime = length(tvec);
survship = Array{Float64}(undef,ltime)
for i=1:ltime
	survship[i] = F(tvec[i]);
end
survship[ltime] = 0;
#Check relationships

#Body mass over time
# R"plot($(tvec/60/60/24/365),$(mass),xlab='Growth time',ylab='Mass')"


#=====================#
#Simulation
#=====================#
tvec_yrs = tvec/60/60/24/365;
#Reproduction parameters
alpha = 2;
beta = 0.1;
maturity = 10; #in years
repsilon = findmin((maturity-tvec_yrs).^2)[2]

# R"""
# par(mfrow=c(2,1))
# plot($(tvec_yrs),$(mass),xlab='Growth time',ylab='Mass')
# lines(rep($(tvec_yrs[repsilon]),length.out=length(seq(1,10^6,by=1000))),seq(1,10^6,by=1000))
# """



gen = 20;
tstep = minimum(tint); #The timesteps are set at the minimum time interval between age classes
tmax = maximum(tvec)*gen; #How much time to simulate?
totsteps = tmax/tstep;
#How many individuals start
n0 = 1000;
# initialstate = zeros(Int64,n0);

#state vector - number of individuals in given state
state = zeros(Int64,lspan);
stateclock = zeros(Float64,lspan);
state[1] = n0; #start out all individuals at birth size class

#teeth lost at each state
toothdrop = zeros(Float64,lspan);
#mass vector for teeth in each timebin
#NOTE this right now it's just linear, but should be allomtric and/or saturating
# mintoothmass = 5.0;
# maxtoothmass = 50.0;
# toothmass = collect(mintoothmass:((maxtoothmass-mintoothmass)/(lspan-1)):maxtoothmass);
bodylength = mass 
toothlength = -0.0800064.*(-26.665 - bodylength);
teethlostperday = 2/40;
#(teeth/s) :: number of teeth lost per second
toothlossrate = (teethlostperday/24/60/60); #*toothmass; 


N=1000; #how many timebins to save to calculate steady state distribution
savestate = zeros(Int64,N,lspan);
n=0;

#Simulation
tcum = 0;
tictoc = 0;
popstate = Array{Int64}(undef,0);
while tcum < tmax
	
	tictoc = tictoc + 1;
	pop = sum(state);
	push!(popstate,pop);
	
	if mod(tictoc,100) == 0
		println("tictoc ",tictoc,"; Population = ",pop)
	end
	
	tcum += tstep;
	stateclock += tstep;
	
	#Individuals grow
	for i=1:ltime
		
		#Remove indidivduals that die
		prmort = 1-survship[i];
		bdist = Binomial(1,prmort);
		die = sum(rand(bdist,state[i]));
		#Remove individuals who die
		state[i] = state[i] - die;
		
		#death tooth drop
		toothdrop[i] += die*((toothlossrate[i]*4) * tstep);
		
		#NOTE: Tooth pulse drop with death
		
		#Grow
		statetime = stateclock[i];
		# if i < lspan
			if state[i] > 0
				if statetime >= tvec[i]
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
			recruits = convert(Int64,floor((alpha*state[i])/(1+beta*state[i])));
			#Add the recruits to the initial state
			state[1] = state[1] + recruits;
		end
		
		#lose teeth
		# number of ndividuals * tooth loss rate [toothloss/sec] * time interval [sec]
		toothdrop[i] += state[i]*(toothlossrate[i] * tstep);
	end
	
	#save states from the last N timestep
	
	if mod(tictoc,N) == 0
		n=0;
	end
	n+=1;
	savestate[n,:] = state;
	
end

#Simulation time in years
maxlifetime = maximum(tvec)/60/60/24/365;
timesim = (collect(1:tictoc)*tstep)/60/60/24/365/maxlifetime;
# R"plot($(timesim),$popstate,type='l',xlab='Time (generations)',ylab='Population size',log='y')"




#Produce stationary distribution of mass NOTE OVER some time interval (rather than a single step)
# mvec = Array{Float64}(0);
# for j=1:N
# 	for i=1:ltime
# 		mvec = [mvec;vec(repeat([mass[i]],inner=savestate[j,i]))];
# 	end
# end
# R"hist($mvec,xlab='Steady state mass distribution (grams)',breaks=20,col='gray',main='',freq=FALSE)"


accumstate = sum(savestate,1);
laststate = findall(iszero,accumstate)[1];
R"""
par=mfrow=c(2,1)
barplot($(accumstate[1:laststate]),names.arg=$((epsilonvec*M)[1:laststate]));
""";



# #Tooth accumulation distribution
# R"""
# plot($(epsilonvec*M),$(toothdrop),xlab='Teeth by body mass (grams)',type='l',lwd=3)
# """

R"""
barplot($(toothdrop),names.arg=$(epsilonvec*M),xlab='Teeth by body mass (grams)');
""";
