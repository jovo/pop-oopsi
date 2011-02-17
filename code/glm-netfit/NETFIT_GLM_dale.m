%SCRIPT TO ESTIMATE W FROM SAMPLE OF SPIKE HISTORIES USING DALE PRIOR
% AN ATTEMPT USING HARD WAY - DO IN A LOOP AND DECIDE ON SIGN ASSIGNMENT BY
% COMPARING FULL MATRIX FIT LIKELIHOODS WITH EITHER SIGN FOR j
par=[];
par.tmin=tmin;                      % min coupling time-depth
par.tmax=tmax;                      % max coupling time-depth
par.dt=netSim.FR;                   % time-step
par.samples=spkM;                   % samples in the data
par.lambda=0;
if(exist('lambda')==1) par.lambda=lambda; end % L1-regularizer strength


%%RECOVER WEIGHTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dalepr=zeros(N,1);                  % initial guess

flg=1;
while(flg)
  dalepr_old=dalepr;
  for ik=1:N
    dalepr(ik)=1;                    % check with pos sign
    par.signs=dalepr;
    NETFIT_GLM_solve
    lik_pos=lik_tot;
    abs_pos=max(Weights(:,ik));
    clear weights_idal
    
    dalepr(ik)=-1;                   % check with neg sign
    par.signs=dalepr;    
    NETFIT_GLM_solve
    lik_neg=lik_tot;
    abs_neg=min(Weights(:,ik));
    
                                    % make decision
    if(lik_pos<lik_neg) dalepr(ik)=1; 
    elseif(lik_pos>lik_neg) dalepr(ik)=-1;
    else dalepr(ik)=0;
    end
    fprintf('===============\n+%.3g vs -%.3g == %.3g/%i (%.3g vs %.3g)\n',...
      lik_pos,lik_neg,(lik_pos-lik_neg)*N^2,N^2,abs_pos,abs_neg);
  end
  
  flg=max(abs(dalepr-dalepr_old));  % continue if changed
end

%clean the mess
save(fname,'J','Rate','Omega','Weights');
clear par flg dalepr_old k lik_pos lik_neg lik_tot ik
    