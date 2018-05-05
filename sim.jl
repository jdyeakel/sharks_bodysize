using Distributions
using HDF5
using JLD
using RCall

#=====================#
#Setup
#=====================#

#Initialize variables
agestep = 100;
agevec = collect(1/agestep:1/agestep:0.99);
lspan = length(agevec);

#We will want to use the Gillooly relationships (Nature, 2002) rather than these (which are for mammals!)
# B0 = 0.047; #Metabolic normalization constant
# C = 19.50; #endotherms (Brown Ecology 2004)
C = 18.47; #FISH
temp = 294.3; #324 to get B0 in natcomm; 291.5 for 65F; 294.3 for 70F
B0 = (exp(C))/(exp(0.63/((8.61733*10^(-5.0))*temp)));
Em = 5774; #Energy needed to synthesize a unit of mass
a = B0/Em;
m0 = 100000; #Birth size (grams)
M = 500000; #Asymptotic size :: tiger shark: 380000 to 630000 grams
eta = 3/4; #Scaling exponent

#Gompertz parameters -- somewhat mysterious need to get Calder book
#from Calder, W. A. Size, function, and life history (Harvard University Press, 1984)
F0 = 1;
a0 = 1.88*10^(-8.0);
b0 = -0.56;
a1 = 1.45*10^(-7.0);
b1 = -0.27;
c0 = a0*M^b0;
c1 = a1*M^b1;

#Growth function :: the time it takes to get to epsilon % body mass M from m0
function ts(epsilon)
	ts = log(((1-(m0/M)^(1-eta))/(1-epsilon^(1-eta))))*((M^(1-eta)/(a*(1-eta))));
	return ts
end
#Survivorship from Gompertz curve
function F(t)
	survship = F0*exp((c0/c1)*(1-exp(c1*t)));
	return survship
end

#The time it takes to get from m0 to epsilon M for increasing epsilons.
#The time experienced by the individual is tvec[3] = tvec[2]-tvec[1]
tvec = Array{Float64}(lspan);
mass = Array{Float64}(lspan);
for i=1:length(agevec)
	epsilon = agevec[i];
	tvec[i] = ts(epsilon);
	mass[i] = epsilon*M;
end
tint = diff(tvec);

# Survivorship: this is 1 - the probability of death at a given ageclass
ltime = length(tvec);
survship = Array{Float64}(ltime)
tvecsum = cumsum(tvec);
for i=1:ltime
	survship[i] = F(tvecsum[i])
end
survship[ltime] = 0;
#Check relationships

#Body mass over time
R"plot($(tvec/60/60/24/365),$(mass),xlab='Growth time',ylab='Mass')"


#=====================#
#Simulation
#=====================#

#Reproduction parameters
alpha = 10;
beta = 0.001;
maturity = 
repepsilon = 50;

gen = 200;
tstep = minimum(tint);
tmax = maximum(tvec)*gen; #How much time to simulate?
totsteps = tmax/tstep;
#How many individuals start
n0 = 1000;
# initialstate = zeros(Int64,n0);

#state vector - number of individuals in given state
state = zeros(Int64,lspan);
stateclock = zeros(Float64,lspan);
state[1] = n0; #start out all individuals at birth size class


#Simulation
tcum = 0;
tictoc = 0;
popstate = Array{Int64}(0);
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
	for i=1:lspan
		
		#Remove indidivduals that die
		prmort = 1-survship[i];
		bdist = Binomial(1,prmort);
		die = sum(rand(bdist,state[i]));
		#Remove individuals who die
		state[i] = state[i] - die;
		
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
		if i >= repepsilon
			#Recruitment function
			recruits = convert(Int64,floor((alpha*state[i])/(1+beta*state[i])));
			#Add the recruits to the initial state
			state[1] = state[1] + recruits;
		end
		
		#Tooth stuff
	end
	
end

#Simulation time in years
timesim = (collect(1:tictoc)*tstep)/60/60/24/365;
R"plot($(timesim),$popstate,type='l',xlab='Time',ylab='Population size',log='y')"

#Produce stationary distribution of mass NOTE OVER some time interval (rather than a single step)
mvec = Array{Float64}(0);
for i=1:length(state)
	mvec = [mvec;vec(repeat([mass[i]],inner=state[i]))];
end
R"hist($mvec,xlab='Mass (grams)',breaks=20,col='gray',main='')"







