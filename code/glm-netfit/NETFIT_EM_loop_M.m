%SCRIPT TO PERFORM INNER EM-LOOP FOR ESTIMATION OF 
%CALCIUM MODELS M ON A SUBSET OF DATA

if(exist('setFQ')~=1) setFQ=1; end      % frame skipping
if(exist('setH')~=1) setH=1; end        % spike-history dimensions

N_spk=100;                              % spikes to use for estimating [Ca] model

Sim=[];
Sim.N       = 50;                       % # of particles
Sim.dt      = netSim.FR/setFQ;          % time step size
Sim.freq    = setFQ;                    % # of time steps between observations
Sim.M       = setH;                     % number of spike history dimensions
Sim.pf      = 1;                        % use conditional sampler (not prior) when possible

Sim.Mstep   = 1;                        % do M-step
Sim.C_params= 1;                        % whether to estimate calcium parameters {tau,A,C_0,sig}
Sim.n_params= 1;                        % whether to estimate rate governing parameters {b,k}
Sim.h_params= 1;                        % whether to estimate spike history parameters {h}
Sim.F_params= 1;                        % whether to estimate observation parameters {alpha,beta}
Sim.G_params= 1;                        % whether to estimate observation parameters {gamma}
Sim.MaxIter = 15;                       % max # of EM iterartions

Sim.Scan    = 1;                        % end-of-frame data?

P=[];
if Sim.M==1                                         % if there are spike history terms
    P.omega = -1;                                   % weight
    P.tau_h = 0.01;                                 % time constant
    P.tau_h=Sim.dt/(1-exp(-Sim.dt/P.tau_h));        % correction -- large dt
    P.sigma_h = 0.01;                               % stan dev of noise
end

K=Sim.dt/netSim.dt;                     %sim/inf resampling constant
for k=nrange{id_proc}
  rate=5;
  if(~isempty(M{k})) rate=exp(M{k}.k(1)); end %expected neuron rate
  rate=max(5,rate);
  T=N_spk/rate;                           % sub-time to use for models estimation
  
  Sim.T       = ceil(T/Sim.dt);           % # of samples
  Sim.Nsec    = Sim.T*Sim.dt;             % # of actual seconds
  Sim.T_o     = Sim.T;                    % # of observations
  Sim.tvec    = Sim.dt:Sim.dt:Sim.Nsec;   % vector of times

  Sim.n=repmat(NaN,[1,Sim.T]);            % for plotting purposes in PF
  if(exist('n_GT')==1) Sim.n(min(Sim.T,ceil(n_GT{k}/K)))=1; end

  Sim.StimDim=1;  
  Sim.x=ones(1,Sim.T);
  if(~isempty(J{k})) Sim.x=J{k}(1:Sim.T); end% estimated incoming currents
  
  %% 2) initialize parameters
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(isempty(M{k}))
    % initialize particle filter parameters    
    P.rate      = rate;                 % expected rate
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
  else                                  % use previous iteration to start off
    P=M{k};
  end

  %% 3) estimate parameters
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  warning off
  Ft=repmat(F{k}(:)',[Sim.freq 1]); Ft=Ft(:);
  [I.M I.P]   = GOOPSI_main_v1_0(Ft(1:Sim.T),P,Sim);
  I.P.FQ=Sim.freq;
  M{k}=I.P;
  fprintf('\n');
  warning on
end

save(fname,'M');

%clean the mess
clear N_spk rate T rf P Sim I
