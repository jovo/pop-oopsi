% this script pulls set of spike trains from fluorescence
clear, close, clc, fprintf('\nnetSim.SPIKE\n')

netsim_name='netSim0315N50S5678mod';
load([netsim_name,'.mat'],'netSim','n','F');

%DOWNSAMPLE
F=F(:,2:2:end);
netSim.FR=2*netSim.FR;

T_M=15;                                 % time to estimate model
T_I=300;                                % time to obtain spike

Sim.dt      = netSim.FR;                % time step size
Sim.freq    = 1;                        % # of time steps between observations
Sim.N       = 50;                       % # of particles
Sim.M       = 0;                        % number of spike history dimensions
Sim.pf      = 1;                        % use conditional sampler (not prior) when possible
Sim.StimDim = 1;                        % # of stimulus dimensions

Sim.Mstep   = 1;                        % do M-step
Sim.C_params= 1;                        % whether to estimate calcium parameters {tau,A,C_0,sig}
Sim.n_params= 1;                        % whether to estimate rate governing parameters {b,k}
Sim.h_params= 0;                        % whether to estimate spike history parameters {h}
Sim.F_params= 1;                        % whether to estimate observation parameters {alpha,beta}
Sim.G_params= 1;                        % whether to estimate observation parameters {gamma}
Sim.MaxIter = 15;                       % max # of EM iterartions

Sim.Scan    = 1;                        % end-of-frame data
Sim.minVar  = [0.1 0.1];                % min variances on nbar & [Ca] for mcmc
                                        % this should be made bigger, eg [0.1 1]

if(exist('Js')~=1) Js=ones(size(F)); end% input currents, supply if known

% 1) infer parameters & connectivity from resampleddata
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Np=size(F,1);                           % number of neurons
P_infer=cell(1,Np);                     % infered models
n_infer=zeros(Np,ceil(T_I/netSim.FR));  % infered spikes
flg_start=1;
for k=31:Np                              % loop over neurons
  if(exist('tmp.mat') && flg_start)     % load previous data, if any
    flg_start=0;
    load tmp.mat;
    k=31;                              % get to next cell
  end

  fprintf('Neuron %i\n############################################\n',k);

  %ESTIMATE CELL MODEL
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  Sim.T       = ceil(T_M/netSim.FR);
  Sim.Nsec    = Sim.T*Sim.dt;             % # of actual seconds
  Sim.T_o     = Sim.T;                    % # of observations
  Sim.tvec    = Sim.dt:Sim.dt:Sim.Nsec;   % vector of times
  Sim.M=0;                                % number of spike history dimensions
  Sim.x       = Js(k,1:Sim.T);            % input current
  Sim.Mstep   = 1;                        % do M-step  
  Sim.MaxIter = 15;                       % max # of EM iterartions  

  Sim.n=zeros(1,Sim.T);                   % for plotting purposes in PF
  Sim.n(max(1,min(Sim.T,round(n{k}*netSim.dt/Sim.dt))))=1;
  Sim.n(Sim.n==0)=NaN;              

  %% 2) initialize parameters
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  P=[];
  P.rate  = 5;                        % expected rate
  % initialize particle filter parameters
  rf          = log(-log(1-P.rate/Sim.T)/Sim.dt);
  P.k         = [rf,zeros(1,Sim.StimDim-1)]';% linear filter
  P.tau       = 0.5;                  % calcium decay time constant (sec)
  P.A         = 50;                   % jump size ($\mu$M)
  P.sig       = 25;                   % standard deviation of noise (\mu M)
  P.C_0       = 25;                   % baseline [Ca++]
  P.n         = 1;                    % hill equation exponent
  P.k_d       = 100^P.n;              % hill coefficient
  P.alpha     = 1;                    % F_max
  P.beta      = 0;                    % F_min
  P.gamma     = 1e-5;                 % scaled variance
  P.zeta      = 4*P.gamma;            % constant variance

  P.C_init    = P.C_0;                % initial [Ca++]
  P.tau_c     = P.tau;
  P.sigma_c   = P.sig;
  P.a         = Sim.dt/P.tau_c;
  P.lam       = Sim.T/(P.rate*P.A)*Sim.dt;% expected jump size ber time bin

  %% 3a) estimate parameters
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  warning off

  [I.M I.P]   = GOOPSI_main_v1_0(F(k,1:Sim.T),P,Sim);
  P1=I.P;
  fprintf('\n');

  % 3b) infer nbar
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  fprintf('Particle Filter sample %i\n############################################\n',k);        
  %get swarms with PF
  Sim.T       = ceil(T_I/netSim.FR);
  Sim.Nsec    = Sim.T*Sim.dt;             % # of actual seconds
  Sim.T_o     = Sim.T;                    % # of observations
  Sim.tvec    = Sim.dt:Sim.dt:Sim.Nsec;   % vector of times
  Sim.x       = Js(k,1:Sim.T);            % input current
  Sim.Mstep   = 0;                        % do M-step  
  Sim.MaxIter = 1;                        % max # of EM iterartions 
  
  Sim.n=zeros(1,Sim.T);
  Sim.n(max(1,min(Sim.T,round(n{k}*netSim.dt/Sim.dt))))=1;
  Sim.n(Sim.n==0)=NaN;                    % for plotting purposes in PF

  P=I.P;
  [I.M I.P]   = GOOPSI_main_v1_0(F(k,1:Sim.T),P,Sim);
  P1.nbar=I.M.nbar;
  P_infer{k}=P1;                      % SAVE CELL MODEL    
  fprintf('\n');  
  warning on
  
  

  %% 4) draw bias-corrected sample with MCMC
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % par structure for mcmc
  P1.grid=25;                         % size of stochastic grid
  P1.iter=10;                         % number of mcmc steps

  P1.wf=I.M.w;                        % copy particle swarms
  for t=1:Sim.T                       % collect centers for KDE
    P1.nf{t}=[I.M.n(:,t) I.M.C(:,t)]; % swarm of particles
    v1=sqrt(I.M.nvar(t));             % variance in nbar
    A=repmat(I.M.C(:,t),1,Sim.N);     % variance in C
    A=(A-A');
    v2=max(P1.sigma_c,std(A(:)));
    P1.vf(:,t)=max(Sim.minVar/sqrt(Sim.N),[v1 v2])';% for use with SMC
  end

  P1.dt=Sim.dt;                       % some auxiliary stuff for speed
  P1.sig2=P1.sigma_c^2*P1.dt;
  P1.rate=max(1e-3,exp(-exp(P1.k'*Sim.x)*P1.dt));
  P1.pf=cell(1,Sim.T);                % dummy, for using spPG sampler

  Y1=cell(1,Sim.T); for t=1:Sim.T Y1{t}=NaN; end;% missing observations
  for t=1:Sim.freq:Sim.T Y1{t}=F(k,t); end% available observations
                                      % this probably gonna suck with
                                      % missing observations

  fprintf('Bias-corrected sample %i\n############################################\n',k);    
  Xr=sampleMCMC(P1,Y1,@spPY_GOOPSI,@spPXX_GOOPSI,@spPG_GOOPSI);% sample spike train
  for t=1:Sim.T  n_infer(k,t)=Xr{t}(1); end;

  fprintf('sample[n|F,]: h:%i|m:%i vs %i spk\n',...
    sum(I.M.nbar>0.5),sum(n_infer(k,:)),sum(n{k}*netSim.dt<T_I));

  save tmp.mat n_infer P_infer k Sim
end

save([netsim_name,'INF33Hz.mat'],'Sim','P_infer','n_infer');