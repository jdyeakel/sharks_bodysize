function empirical_sim_comparison(toothdrop,toothlength1,measures)


    mvj, mva, varj, vara, pj, pa, peakjuvbin, peakadultbin, peakjuvdist, peakadultdist = toothdist_analysis(toothdrop,toothlength1);


    #derive empirical toothdrop and toothlength vectors
    U = kde(measures);
    emp_toothlength = U.x;
    emp_toothdrop = U.density;
    mve, vare, pe, peakebin, peakedist = toothdist_emp_analysis(emp_toothdrop,emp_toothlength);

    #Measure the distance
    meandist

    modedist1

    modedist2


