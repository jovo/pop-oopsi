%SCRIPT TO PERFORM INNER EM-LOOP FOR ESTIMATION OF
%CALCIUM MODELS M ON A SUBSET OF DATA

if( exist('setFQ') ~=1 ) setFQ=1; end      % frame skipping
if( exist('setH')  ~=1 ) setH=1; end        % spike-history dimensions
if( exist('rate_b') == 1 ) rate_b=5; end

N_spk=100;                              % spikes to use for estimating [Ca] model
rate=rate_b;                            %expected spike rate

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

if(exist('holdTau') == 1) Sim.holdTau=holdTau; end % hold tau in M-step

Sim.Scan    = 1;                        % end-of-frame data?

P=[];
if Sim.M>0                                          % if there are spike history terms
    P.omega = -1;                                   % weight
    P.tau_h = 0.01;                                 % time constant
    P.tau_h=Sim.dt/(1-exp(-Sim.dt/P.tau_h));        % correction -- large dt
    P.sigma_h = 0.01;                               % stan dev of noise
end

K=Sim.dt/netSim.dt;                     %sim/inf resampling constant
for k=nrange{id_proc}
    rate1=rate;
    if(exist('cP')) rate1=cP.rate; end    %expected spike rate from init
    if(~isempty(M{k})) rate1=exp(M{k}.k(1)); end %expected spike rate from before

    T=N_spk/rate1;                           % sub-time to use for models estimation

    Sim.T       = min(length(F{k}),ceil(T/Sim.dt));% # of samples
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
    if(isempty(M{k}) && exist('cP')==1)   %load custom initialization
        P=cP;
    elseif(isempty(M{k}))

        % initialize particle filter parameters
        P.rate      = rate;                 % expected rate
        P.k         = log(P.rate);          % linear filter
        P.A         = 50;                   % jump size ($\mu$M)
        P.C_0       = 25;                   % baseline [Ca++]
        P.n         = 1;                    % hill equation exponent
        P.k_d       = 100^P.n;              % hill coefficient
        P.alpha     = 1;                    % F_max
        P.beta      = 0;                    % F_min
        P.gamma     = 1e-5;                 % scaled variance
        P.zeta      = 4*P.gamma;            % constant variance

        P.C_init    = P.C_0;                % initial [Ca++]
        P.tau_c     = 0.5;
        P.sigma_c   = 25;
        P.a         = Sim.dt/P.tau_c;

        % initialize fast-oopsi parameters
        %     P.sig       = 25;                   % standard deviation of noise (\mu M)
        %     P.lam       = Sim.T/(rate*P.A)*Sim.dt;% expected jump size ber time bin
        %     P.tau       = 0.5;                  % calcium decay time constant (sec)

    else                                  % use previous iteration to start off
        P=M{k};
        if(Sim.n_params == 0) P.k=1; end    %if k is not est, set at  P.k=1 since J is exact
    end

    %% 3) estimate parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    warning off
    Ft=repmat(F{k}(:)',[Sim.freq 1]); Ft=Ft(:);
    %   Sim.FastInit=0;
    Sim.b_est=1;
    [I.M I.P]   = GOOPSI_main_v1_0(Ft(1:Sim.T),P,Sim);
    I.P.FQ=Sim.freq;
    M{k}=I.P;
    fprintf('\n');
    warning on
end

save(fname,'M');

%clean the mess
clear N_spk rate T rf P Sim I
