function ts(epsilon1,epsilon2,a,eta)
    # ts = log(((1-(m0/M)^(1-eta))/(1-epsilon^(1-eta))))*((M^(1-eta)/(a*(1-eta))));
    # Epsilon 1: proportion of M to start
    # Epsilon 2: proportion of M to end
    ts = log(((1-(epsilon1)^(1-eta))/(1-epsilon2^(1-eta))))*((M^(1-eta)/(a*(1-eta))));
    return ts
 end
