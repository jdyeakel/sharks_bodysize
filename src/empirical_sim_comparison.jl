function empirical_sim_comparison(toothdrop,toothlength,measures)


    mvj, mva, varj, vara, pj, pa, peakjuvbin, peakadultbin, peakjuvdist, peakadultdist, modesj, modesa = toothdist_analysis(toothdrop,toothlength);


    #derive empirical toothdrop and toothlength vectors
    U = kde(measures);
    emp_toothlength = U.x;
    emp_toothdrop = U.density;
    mve, vare, pe, peakebin, peakedist, modes = toothdist_emp_analysis(emp_toothdrop,emp_toothlength);

    #Error function
    function errorfunc(pred,obs)
        err = sum(((pred .- obs) ./ obs));
        return err
    end

    #Measure the distance
    meanchij = errorfunc(mvj,mve);
    meanchia = errorfunc(mva,mve);

    #HA! THIS IS CHI SQUARE STAT
    # meandist_juv = 1 - (mve - sqrt((mvj - mve)^2))/mve; 
    # meandist_adult = 1 - (mve - sqrt((mva - mve)^2))/mve;

    if peakebin == 0
        # modechij = (maximum(modesj) - maximum(modes))/maximum(modes);
        modechij = errorfunc(maximum(modesj),maximum(modes));

        # modechia = (maximum(modesa) - maximum(modes))/maximum(modes);
        modechia = errorfunc(maximum(modesa),maximum(modes));

        #Binary matching
        if peakjuvbin == 0
            modedistchij = 0.;
        else
            modedistchij = 1.;
        end
        if peakadultbin == 0
            modedistchia = 0.;
        else
            modedistchia = 1.;
        end
    else
        # modechij = ((minimum(modesj) - minimum(modes))/minimum(modes)) + ((maximum(modesj) - maximum(modes))/maximum(modes));
        modechij = errorfunc(sort(modesj),sort(modes));

        # modechia = ((minimum(modesa) - minimum(modes))/minimum(modes)) + ((maximum(modesj) - maximum(modes))/maximum(modes));
        modechia = errorfunc(sort(modesa),sort(modes));


        # modedistchij = (peakjuvdist - peakedist)/peakedist;
        modedistchij = errorfunc(peakjuvdist,peakedist);
        # modedistchia = (peakadultdist - peakedist)/peakedist;
        modedistchia = errorfunc(peakadultdist,peakedist);

    end

    return meanchij, meanchia, modechij, modechia, modedistchij, modedistchia

end





