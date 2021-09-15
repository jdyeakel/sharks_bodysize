function plotdensityreturn(filename,indices,site_type)
    namespace = smartpath(filename,indices);
    @load namespace mass1 mass2 toothdrop toothlength1 toothlength2;
    toothlength = toothlength1[1,:];
    
    if site_type == "juv"
        simdensity = ((toothdrop[:,1])/sum(toothdrop[:,1])) #/maximum(toothdrop[:,1]);
    else
        simdensity = ((toothdrop[:,2])/sum(toothdrop[:,2])) #/maximum(toothdrop[:,1]);
    end

    scaledsimdensity = simdensity ./ maximum(simdensity);
    # else
    #     toothlength = collect(0:1:40);
    #     scaledsimdensity = repeat([0],inner=length(toothlength));
    # end


    return toothlength, scaledsimdensity

end