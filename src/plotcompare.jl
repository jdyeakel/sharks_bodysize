function plotcompare(M,filename,measures)

    matchparams = findall(x->x==1,M)
    rmatch = rand(matchparams);
    
    datadensity = kde(measures);
    
    indices = [1,rmatch[1],rmatch[2]];
    namespace = smartpath(filename,indices);
    @load namespace mass1 mass2 toothdrop toothlength1 toothlength2;
    toothlength = toothlength1[1,:];
    
    simdensity = ((toothdrop[:,1])/sum(toothdrop[:,1])) #/maximum(toothdrop[:,1]);
    scaledsimdensity = simdensity ./ maximum(simdensity);
    

    return datadensity, toothlength, scaledsimdensity

end