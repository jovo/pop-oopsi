%SCRIPT TO DRAW SAMPLE OF SPIKE HISTORIES FOR ESTIMATING W
% --- neal's MCMC + Gibbs
% ---- this doesn't support skipped frame-rates for now

par=[];
par.Grid_points = 25;                   % points in Grid in neal's MCMC
par.Grid_var    = [0.1 0.1];            % min Grid variance
par.MCMC_steps  = 5;                    % MCMC steps before sampling
par.Gibbs_steps = 1*spkM;               % Gibbs steps
par.Gibbs_freq  = 1;                    % Gibbs sampling freq

T=length(F{1});                         % recording time

if(cnt==0)                              % in the first round everything
  par.Gibbs_steps=spkM;                 % is independent, no need to relax
  par.Gibbs_freq=1;
end

if(exist('H_loc')~=1)                   % internal Gibbs state
  H_loc=cell(N,1);                      % H is spike history for GLM
  for k=1:N H_loc{k}=zeros(1,T); end
end


NETFIT_sampler_neal_gibbsgencur;        % prepare all currents


gcnt=0;
H_new=cell(N,1);
syncdata={'H_loc'};                     % equalize proc-lists !MUST DO!
xrange=nrange;
for k=length(nrange{end})+1:length(nrange{1}) nrange{end}(k)=0; end
while(gcnt<par.Gibbs_steps)               % OUTER GIBBS LOOP
  for k=nrange{id_proc}
    if(k==0)                              % need do nothing, dummy
      NETFIT_sync;
      continue;
    end      
                                          % inner MCMC:
    NETFIT_sampler_neal_mcmcpart;         % one neuron spike train
    H_loc{k}=H_new{k};                    % replace spike train    
    
    NETFIT_sync                           % SYNCHRONIZE
    
    NETFIT_sampler_neal_gibbsgencur;      % update currents
  end
  
  gcnt=gcnt+1;                            % accumulate samples/Gibbs
  if( mod(gcnt,par.Gibbs_freq) == 0 )       
    icnt=round(gcnt/par.Gibbs_freq);
    for kk=1:N
      n{kk}(icnt,:)=H_loc{kk};
      H{kk}(icnt,:)=H_loc{kk};
    end
  end
  
end
nrange=xrange;

save(fname,'n','H','H_loc','J');

%clean the mess
clear par T gcnt icnt H_new
