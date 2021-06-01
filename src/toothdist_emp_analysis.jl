function toothdist_emp_analysis(emp_toothdrop,emp_toothlength)

    #MEASURE SOMETHING
    mve = dot(emp_toothlength,(emp_toothdrop/sum(emp_toothdrop)));
    
    # meanjuv[r,temp_pos,dist_pos,sigtau_pos,tau_pos] = mvj;
    # meanadult[r,temp_pos,dist_pos,sigtau_pos,tau_pos] = mva;
    
    #NOTE: collect later
    # meanjuv[r,sigtau_pos,tau_pos] = mvj;
    # meanadult[r,sigtau_pos,tau_pos] = mva;
    
    #calculate variance
    meanofxsquared = dot(emp_toothlength.^2,(emp_toothdrop/sum(emp_toothdrop)));
    squareofmeanofx = mve^2;
    vare = meanofxsquared - squareofmeanofx;
    # varjuv[r,temp_pos,dist_pos,sigtau_pos,tau_pos] = varj;
    #NOTE: collect later
    # varjuv[r,sigtau_pos,tau_pos] = varj;
    
   
    pe = (emp_toothdrop/sum(emp_toothdrop));

    #NOTE: place modality analysis as function
    secondpeake, secondtoothe, troughtoothe, troughmine, maxtoothe = modality_analysis(pe,collect(emp_toothlength));

    
    #Record presence of multiple modes
    if secondpeake == 0.0
        # peakjuv[r,temp_pos,dist_pos,sigtau_pos,tau_pos] = 0;
        # peakjuvquant[r,temp_pos,dist_pos,sigtau_pos,tau_pos] = 0;
        #NOTE: collect later
        # peakjuv[r,sigtau_pos,tau_pos] = 0;
        # peakjuvquant[r,sigtau_pos,tau_pos] = 0;

        peakebin = 0;
        peakedist = 0;

    else
        # peakjuv[r,temp_pos,dist_pos,sigtau_pos,tau_pos] = 1;
        # peakjuvquant[r,temp_pos,dist_pos,sigtau_pos,tau_pos] = abs(secondtoothj-maxtoothj);
        #NOTE: collect later
        # peakjuv[r,sigtau_pos,tau_pos] = 1;
        # peakjuvquant[r,sigtau_pos,tau_pos] = abs(secondtoothj-maxtoothj);

        peakebin = 1;
        peakedist = abs(secondtoothe-maxtoothe);
    end

    modes = sort([secondtoothe,maxtoothe]);

    return mve, vare, pe, peakevbin, peakedist, modes





end
