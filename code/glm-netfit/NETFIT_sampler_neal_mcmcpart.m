%SCRIPT TO DRAW SAMPLE OF SPIKE HISTORIES FOR ESTIMATING W
% --- neal's MCMC + Gibbs
% SAMPLE ONE NEURON USING MCMC, FOR GIBBS LOOP

%%PF-pass
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sim=[];
Sim.N       = 25;                       % # of particles
Sim.dt      = netSim.FR;                % time step size
Sim.freq    = 1;                        % # of time steps between observations
Sim.M       = 1;                        % number of spike history dimensions
Sim.pf      = 1;                        % use conditional sampler (not prior) when possible

                                        % don't infer anything here
Sim.Mstep   = 0;                        % do M-step
Sim.C_params= 0;                        % whether to estimate calcium parameters {tau,A,C_0,sig}
Sim.n_params= 0;                        % whether to estimate rate governing parameters {b,k}
Sim.h_params= 0;                        % whether to estimate spike history parameters {h}
Sim.F_params= 0;                        % whether to estimate observation parameters {alpha,beta}
Sim.G_params= 0;                        % whether to estimate observation parameters {gamma}
Sim.MaxIter = 1;                        % max # of EM iterartions

Sim.Scan    = 1;                        % end-of-frame data

Sim.T       = length(F{k});             % # of samples
Sim.Nsec    = Sim.T*Sim.dt;             % # of actual seconds
Sim.T_o     = Sim.T;                    % # of observations
Sim.tvec    = Sim.dt:Sim.dt:Sim.Nsec;   % vector of times

Sim.n=repmat(NaN,[1,Sim.T]);            % for plotting purposes in PF
if(exist('n_GT')==1) Sim.n(min(Sim.T,ceil(n_GT{k}/netSim.K)))=1; end

%% 2) initialize parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P=M{k};
if(isfield(P,'omega')) Sim.M=1; end

Sim.StimDim=1;  
Sim.x=ones(1,Sim.T);
if(~isempty(J{k})) Sim.x=J{k}(1:Sim.T); end % past-estimated incoming currents


%% 3) PF pass get priors for stochastic-grids
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
warning off
tic;[I.M I.P]   = GOOPSI_main_v1_0(F{k},P,Sim);toc
fprintf('\n');
warning on


%% 4) Corrected sample with nealMCMC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% par structure for mcmc
P1=I.P;

P1.grid=par.Grid_points;            % size of stochastic grid
P1.iter=par.MCMC_steps;             % number of mcmc steps

P1.wf=I.M.w;                        % copy particle swarms
for t=1:Sim.T                       % collect centers for KDE
  P1.nf{t}=[I.M.n(:,t) I.M.C(:,t)]; % swarm of particles
  v1=sqrt(I.M.nvar(t));             % variance in nbar
  A=repmat(I.M.C(:,t),[1,Sim.N]);   % variance in C
  A=(A-A');
  v2=std(A(:));
  P1.vf(:,t)=max(par.Grid_var,[v1 v2])';%for use with MCMC
end

P1.dt=Sim.dt;                       % some auxiliary stuff for speed
P1.sig2=P1.sigma_c^2*P1.dt;
P1.rate=max(1e-3,exp(-exp(P1.k'*Sim.x)*P1.dt));
P1.pf=cell(1,Sim.T);                % dummy, for using spPG sampler

Y1=cell(1,Sim.T);
for t=1:Sim.T Y1{t}=F{k}(t); end    % BEWARE, available observations
                                    % this doesn't do right now
                                    % missing observations, BEWARE

fprintf('nealMCMC for %i\n###################################\n',k);

tic                                 % sample spike train
Xr=sampleMCMC(P1,Y1,@spPY_GOOPSI,@spPXX_GOOPSI,@spPG_GOOPSI);
toc

H_new{k}=zeros(1,Sim.T);            % new spike train for k
for t=1:Sim.T  H_new{k}(t)=Xr{t}(1); end;
fprintf('GOT %i vs %i\n',sum(H_new{k}),length(n_GT{k}));

clear Sim P I P1 v1 A v2 Y1 t Xr
