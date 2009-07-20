%SCRIPT TO ESTIMATE W FROM SAMPLE OF SPIKE HISTORIES
par=[];
par.tmin=tmin;                        % min coupling time-depth
par.tmax=tmax;                        % max coupling time-depth
par.dt=netSim.FR;                     % time-step
par.samples=spkM;                     % samples in the data
par.Tp=Tp;                            % coupling temporal depth
if(exist('bbox')) par.box=bbox; else par.box=[]; end % bounding box for solver
if(exist('rate_b')==1) par.rate=rate_b; else par.rate=5; end
if(exist('lambda')==1) par.lambda=lambda; end % L1-regularizer strength
if(exist('dalepr')==1) par.signs=dalepr; end  % dale prior signs

%%RECOVER WEIGHTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Weights,Omega,Rate,J,ivals]=NETFIT_GLM_solve(n,H,nrange{id_proc},par,ivals);

%clean the mess
save(fname,'J','Rate','Omega','Weights','ivals');
clear par lik_tot