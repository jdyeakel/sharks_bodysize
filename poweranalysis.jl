if homedir() == "/home/z840"
    loadfunc = include("$(homedir())/sharks_bodysize/src/loadfuncs.jl");
    filename = "SandTiger_all.csv";
    namespace = string("$(homedir())/sharks_bodysize/",filename);
    data = CSV.read(namespace,header=true,DataFrame)
else
    loadfunc = include("$(homedir())/Dropbox/PostDoc/2018_sharks/src/loadfuncs.jl");
    filename = "SandTiger_all.csv";
    namespace = string("$(homedir())/Dropbox/PostDoc/2018_sharks/",filename);
    data = CSV.read(namespace,header=true,DataFrame)
end

ndata = names(data);
num = length(ndata);
reps = 1000;
samplesizevec = collect(10:5:1000);
lsamples = length(samplesizevec);

meandiff = SharedArray{Float64}(num,lsamples,reps);
mode1diff = SharedArray{Float64}(num,lsamples,reps);
mode2diff = SharedArray{Float64}(num,lsamples,reps);


for i = 1:num
    measures = Array{Float64}(data[!,i][findall(!ismissing,data[!,i])]);

    U = kde(measures);
    emp_toothlength = U.x;
    emp_toothdrop = U.density;
    mve, vare, pe, peakebin, peakedist, modes, mediane, quartile25e, quartile75e = toothdist_emp_analysis(emp_toothdrop,emp_toothlength);

    eden = (emp_toothdrop/sum(emp_toothdrop));

    #sampling scheme

    
    for j=1:lsamples
        #draw samples from empirical distirbution
        numsamples = samplesizevec[j];
        
        @sync @distributed for r=1:reps
            #sample the empirical distirbution

            

            sampden = rand(numsamples);

            smeasures = Array{Float64}(undef,numsamples);

            for z = 1:numsamples

                pos = findall(x->x>sampden[z], cumsum(eden));
                # if length(pos) == 0
                #     new_sampden = rand();
                #     pos = findall(x->x>(new_sampden[1]), cumsum(eden));
                # end
                smeasures[z] = emp_toothlength[pos[1]];
            end

            sU = kde(smeasures);
            semp_toothlength = sU.x;
            semp_toothdrop = sU.density;
            smve, svare, spe, speakebin, speakedist, smodes, smediane, squartile25e, squartile75e = toothdist_emp_analysis(semp_toothdrop,semp_toothlength);
        
            #difference from means and modes
            
            meandiff[i,j,r] = sqrt((mve - smve)^2);
            if modes[1] == 0.
                mode1diff[i,j,r] = 0.;
            else
                mode1diff[i,j,r] = sqrt((modes[1] - smodes[1])^2);
            end
            mode2diff[i,j,r] = sqrt((modes[2] - smodes[2])^2);

        end
    end
end

avg_meandiff = mean(meandiff,dims=3)[:,:,1];
avg_mode1diff = mean(mode1diff,dims=3)[:,:,1];
avg_mode2diff = mean(mode2diff,dims=3)[:,:,1];

filename = "figures/fig_poweranalysis.pdf";
namespace = smartpath(filename);
R"""
pdf($namespace,width=12,height=6)
par(mfrow=c(2,3))
plot($samplesizevec,$(avg_meandiff[1,:]),type='l',ylim=c(0.1,10),main=$(ndata[1]),xlab='Sample size',ylab='Difference',log='y')
lines($samplesizevec,$(avg_mode1diff[1,:]),col='red',lty=2)
lines($samplesizevec,$(avg_mode2diff[1,:]),col='red',lty=1)
"""
for i=2:num
    R"""
    plot($samplesizevec,$(avg_meandiff[i,:]),type='l',ylim=c(0.1,max(c(10,max($(avg_mode1diff[i,:]))))),main=$(ndata[i]),xlab='Sample size',ylab='Difference',log='y')
    lines($samplesizevec,$(avg_mode1diff[i,:]),col='red',lty=2)
    lines($samplesizevec,$(avg_mode2diff[i,:]),col='red',lty=1)
    """
end
R"""
dev.off()
"""