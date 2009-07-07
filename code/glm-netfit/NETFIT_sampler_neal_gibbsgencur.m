%SCRIPT TO DRAW SAMPLE OF SPIKE HISTORIES FOR ESTIMATING W
% --- neal's MCMC + Gibbs
%GENERATE CURRENTS FOR ALL NEURONS
for k=1:N
  J{k}=zeros(1,T);
  for t=tmin:tmax
    offset=t-tmin+1;
    
    for l=1:N
      J{k}=J{k}+Weights(k,l,offset)*[zeros(1,t),H_loc{l}(1:end-t)];
    end
  end
end

%clean the mess
clear k l
