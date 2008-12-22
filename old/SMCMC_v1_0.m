% this script implements the particle-filter-gibbs-metropolis sampler
% 
% 1) set simulation metadata
% 2) initialize parameters
% 3) simulate data
% 4) SMCMC
% 5) plot stuff
% 
% Remarks:
% v1_0: only has simulation for a single cell

clear all; clc

%% 1) set simulation metadata
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Sim.T       = 500;                                  % # of time steps
Sim.dt      = 1/200;                                % time step size
Sim.freq    = 1;                                    % # of time steps between observations
Sim.T_o     = Sim.T/Sim.freq;                       % # of observation time steps
Sim.N       = 20;                                   % # of particles
Sim.M       = 1;                                    % number of spike history dimensions
Sim.pf      = 1;                                    % use conditional sampler (not prior) when possible
Sim.StimDim = 1;                                    % # of stimulus dimensions
Sim.x       = 1*randn(1,Sim.T);                     % stimulus

Sim.Mstep   = false;                                % do M-step
Sim.MaxIter = 0;                                    % max # of EM iterartions
Sim.C_params = false;                               % whether one must estimate calcium parameters

%% 2) initialize parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rate        = 10;                                   % expected spike rate
P.k         = log(-log(1-rate*Sim.dt)/Sim.dt);      % linear filter
P.tau_c     = 0.5;                                  % calcium decay time constant (sec)
P.A         = 20;                                   % jumps size (\mu M)
P.C_0       = 20;                                   % baseline [Ca++] (\mu M)
P.C_init    = P.C_0;                                % initial [Ca++] (\mu M)
P.sigma_c   = 1;                                    % noise on 
P.n         = 1.0;                                  % hill equation exponent
P.k_d       = 200;                                  % hill coefficient
P.alpha     = 1;                                    % F_max
P.beta      = 0;                                    % F_min
P.gamma     = 4e-6;                                 % scaled variance
P.zeta      = 4*P.gamma;                            % constant variance
P.a         = Sim.dt/P.tau_c;

if Sim.M==1                                         % if there are spike history terms
    P.omega = -0.3;                                 % weight
    P.tau_h = 0.2;                                 % time constant
    P.sigma_h = 0.01;                               % stan dev of noise
end

%% 3) simulate data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kx        = P.k'*Sim.x;                             % external input to neuron
epsilon_c = P.sigma_c*sqrt(Sim.dt)*randn(1,Sim.T);  % generate noise on calcium
U_sampl   = rand(1,Sim.T);                          % generate random number to use for sampling
if Sim.M==1                                         % if spike history terms, recursively
    p         = zeros(1,Sim.T);                     % extize p_t because it must be updated iteratively
    n         = zeros(1,Sim.T);                     % extize n_t because it must be updated iteratively
    h         = zeros(1,Sim.T);                     % extize spike history because it must be updated iteratively
    epsilon_h = repmat(P.sigma_h*sqrt(Sim.dt),1,Sim.T).*randn(1,Sim.T); % generate noise on spike history
    for t=2:Sim.T                                   % update states
        h(:,t)= (1-Sim.dt./P.tau_h).*h(:,t-1)+n(t-1) + epsilon_h(:,t); % update h terms
        y_t=kx(t)+P.omega'*h(:,t);                  % generate operand for rate function
        p(t)=1-exp(-exp(y_t)*Sim.dt);               % generate rate
        n(t)  = U_sampl(t)<p(t);                    % sample from bernoulli with prob p_t
    end %time loop
else
    n=rand(Sim.T,1)<1-exp(-exp(P.k)*Sim.dt);
end
C=zeros(1,Sim.T);
epsilon_c = P.sigma_c*sqrt(Sim.dt)*randn(1,Sim.T);  % generate noise on calcium
for t=2:Sim.T                                       % recursively update calcium
    C(t)  = (1-Sim.dt/P.tau_c)*C(t-1) + P.a*P.C_0 + P.A*n(t) + epsilon_c(t);
end
S=Hill_v1(P,C);
F=(P.alpha*S+P.beta+sqrt(P.gamma*S+P.zeta).*randn(1,Sim.T));
F(F<=0)=eps;
figure(1), clf
subplot(211), plot(F+1); hold on, stem(n),
if Sim.M==1, plot(P.omega*h); end
subplot(212),
if Sim.M==1, plot(p), hold on, end
plot(Sim.x,'k');
Sim.n = double(n); Sim.n(Sim.n==0)=NaN;          % for plotting purposes in ParticleFiltD

disp([{'expected rate=', num2str(rate)}; {'actual rate=',num2str(sum(n)/(Sim.T*Sim.dt))}])              % # spikes

%% 4) generate samples

profile on
[I.M I.P]   = GOOPSI_main_v1_0(F,P,Sim);
profile off; profile viewer