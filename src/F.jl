function F(t,m)
    # F0 = 1;
    # a0 = 1.88*10^(-8.0);
    # b0 = -0.56;
    # a1 = 1.45*10^(-7.0);
    # b1 = -0.27;
    # c0 = a0*M^b0;
    # c1 = a1*M^b1;
    # survship = minimum([1,F0*exp((c0/c1)*(1-exp(c1*t)))]);
    survship = exp(-t*m)
    return survship
  end
