function gillespie(N0,r,m,tint)
    
    #Keep track of temperature changes by days
    daysinyear = 365;
    secondsinday = 24*60*60;
    secondsinyear = daysinyear*secondsinday;
    temptime = mean(tempvec) .+ (maximum(tempvec) .- mean(tempvec)).*sin.((pi/(daysinyear/2)).*collect(0:1:daysinyear));
    
    
    let n=0, tcum = 0,tictoc = 0
        while tcum < tmax
        
            #What is the temperature?
            day = round(Int64,daysinyear*((tcum/secondsinyear) - floor(tcum/secondsinyear)));
            day = maximum([1,day]);
            daytemp = temptime[day];
            #find index for closest match in tempvec
            k = findmin((tempvec .- daytemp).^2)[2];
            
            probline = cumsum()
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            #determine rate
            Rate = N[:,t] .* (r .+ m .+ (1/tint[temp,:]))