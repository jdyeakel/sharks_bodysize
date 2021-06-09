function plotcompare(M,qM,filename,measures,r,site_type)
    datadensity = kde(measures);
    matchparams = findall(x->x==1,M)
    if length(matchparams) > 0
        rmatch = findmin(qM)[2];

        indices = [r,rmatch[1],rmatch[2]];
        namespace = smartpath(filename,indices);
        @load namespace mass1 mass2 toothdrop toothlength1 toothlength2;
        toothlength = toothlength1[1,:];
        
        if site_type == "juv"
            simdensity = ((toothdrop[:,1])/sum(toothdrop[:,1])) #/maximum(toothdrop[:,1]);
        else
            simdensity = ((toothdrop[:,1])/sum(toothdrop[:,1])) #/maximum(toothdrop[:,1]);
        end
        
        scaledsimdensity = simdensity ./ maximum(simdensity);
    else
        toothlength = collect(0:1:40);
        scaledsimdensity = repeat([0],inner=length(toothlength));
    end


    return datadensity, toothlength, scaledsimdensity

end