function findlocalminima(signal::Vector)
   inds = Int[]
   if length(signal)>1
       if signal[1]<signal[2]
           push!(inds,1)
       end
       for i=2:length(signal)-1
           if signal[i-1]>signal[i]<signal[i+1]
               push!(inds,i)
           end
       end
       if signal[end]>signal[end-1]
           push!(inds,length(signal))
       end
   end
   inds
 end
 