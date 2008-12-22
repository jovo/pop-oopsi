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
% v1_1: simulates for 2 cells (they both have identical parameters)

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
Sim.x       = ones(1,Sim.T);                        % stimulus
Sim.Ncells  = 2;                                    % # of cells

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

S.h     = zeros(Sim.Ncells,Sim.T);
S.n     = zeros(Sim.Ncells,Sim.T);
S.C     = zeros(Sim.Ncells,Sim.T);
S.F     = zeros(Sim.Ncells,Sim.T);

kx      = P.k'*Sim.x;                             % external input to neuron
eps_c   = P.sigma_c*sqrt(Sim.dt)*randn(Sim.Ncells,Sim.T);  % generate noise on calcium
U_sampl = rand(Sim.Ncells,Sim.T);                          % generate random number to use for sampling

eps_h   = repmat(P.sigma_h*sqrt(Sim.dt),Sim.Ncells,Sim.T).*randn(Sim.Ncells,Sim.T); % generate noise on spike history
eps_c   = P.sigma_c*sqrt(Sim.dt)*randn(Sim.Ncells,Sim.T);  % generate noise on calcium
for t=2:Sim.T                                   % update states
    S.h(:,t)= (1-Sim.dt./P.tau_h).*S.h(:,t-1)+sum(S.n(:,t-1)) + eps_h(:,t); % update h terms
    y_t     = kx(t)+P.omega'*S.h(:,t);                  % generate operand for rate function
    p(:,t)  = 1-exp(-exp(y_t)*Sim.dt);               % generate rate
    S.n(:,t)= U_sampl(:,t)<p(:,t);                    % sample from bernoulli with prob p_t
    S.C(:,t)= (1-Sim.dt/P.tau_c)*S.C(:,t-1) + P.a*P.C_0 + P.A*S.n(:,t) + eps_c(:,t);
end

s   = Hill_v1(P,S.C);
S.F = (P.alpha*s+P.beta+sqrt(P.gamma*s+P.zeta).*randn(Sim.Ncells,Sim.T));
S.F(S.F<=0) = eps;

figure(1), clf
subplot(311), plot(z1(S.F)+1); hold on, stem(S.n'), axis('tight')
subplot(312), stem(S.n'), hold on, plot(P.omega*S.h'), axis('tight')
subplot(313), plot(p'), %hold on, plot(Sim.x,'k'), axis('tight')

Sim.n = double(S.n); Sim.n(Sim.n==0)=NaN;          % for plotting purposes in ParticleFiltD

disp([{'expected rate=', num2str(rate)}; {'actual rate=',num2str(sum(S.n(:))/(Sim.T*Sim.dt)/Sim.Ncells)}])              % # spikes

%% 4) generate samples

% [I.M I.P]   = GOOPSI_main_v1_0(F,P,Sim);
