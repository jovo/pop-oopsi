% this script implements the particle-filter-gibbs-metropolis sampler
%
% 1) set simulation metadata
% 2) initialize parameters
% 3) simulate data
% 4) plot simulation (ie, truth)
% 4) SMCMC
% 5) plot results (ie, inference)
%
% Remarks:
% v1_0: only has simulation for a single cell
% v1_1: simulates for 2 cells (they both have identical parameters)
% v1_2: 2 cell simulation working with stupid code

clear; clc

%% 1) set simulation metadata
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Sim.T       = 500;                                  % # of time steps
Sim.dt      = 1/100;                                % time step size
Sim.Np      = 20;                                   % # of particles
Sim.x       = ones(1,Sim.T);                        % stimulus
Sim.Nc      = 2;                                    % # of cells

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
P.tau_h     = 0.2;                                  % time constant
P.sigma_h   = 0.01;                                 % stan dev of noise
P.omega     = [.2 -.3; .2 -1.5];                     % weights

%% 3) simulate data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S(1).h     = zeros(1,Sim.T);
S(1).n     = zeros(1,Sim.T);
S(1).C     = zeros(1,Sim.T);
S(1).F     = zeros(1,Sim.T);
for ii=2:Sim.Nc; S(ii)=S(1); end

kx      = P.k'*Sim.x;                                   % external input to neuron
eps_c   = P.sigma_c*sqrt(Sim.dt)*randn(Sim.Nc,Sim.T);   % generate noise on calcium
U_sampl = rand(Sim.Nc,Sim.T);                           % generate random number to use for sampling
eps_h   = repmat(P.sigma_h*sqrt(Sim.dt),Sim.Nc,Sim.T).*randn(Sim.Nc,Sim.T); % generate noise on spike history
eps_F   = randn(Sim.Nc,Sim.T);                          % generate noise on fluorescence
p       = zeros(Sim.Nc,Sim.T);                          % prob of spiking for each cell at each tim
y       = zeros(Sim.Nc,Sim.T);                          % input to each cell at each time

for t=2:Sim.T                                           % update states
    for i=1:Sim.Nc                                      % loop over presynaptic cells
        S(i).h(t)   = (1-Sim.dt./P.tau_h).*S(i).h(t-1) + S(i).n(t-1) + eps_h(i,t); % update h terms
        y(i,t)      = kx(t);                            % initialize input to cell
        for j=1:Sim.Nc                                  % loop of post-synaptic cells
            y(i,t)  = y(i,t)+P.omega(i,j)*S(j).h(t);    % generate operand for rate function
        end
        p(i,t)      = 1-exp(-exp(y(i,t))*Sim.dt);       % generate prob of spiking
        S(i).n(t)   = U_sampl(i,t)<p(i,t);              % sample from bernoulli with prob p_t
        S(i).C(t)   = (1-Sim.dt/P.tau_c)*S(i).C(t-1) +...
            (Sim.dt/P.tau_c)*P.C_0 + P.A*S(i).n(t) + eps_c(i,t); %update calcium
        s           = Hill_v1(P,S(i).C(t));             % compute saturated calcium
        S(i).F(t)   = (P.alpha*s+P.beta)+sqrt(P.gamma*s+P.zeta).*eps_F(i,t); % update fluorescence
        if S(i).F(t)<=0; S(i).F(t) = eps; end           % keep fluorescence positive
    end
end

%% 4) plot simulation results
col = [1 0 0; 0 0 1];          % define colors for mean
figure(1), clf
for i=1:Sim.Nc
    subplot(311), plot(z1(S(i).F)+1,'Color',col(i,:)); hold on, stem(S(i).n,'Color',col(i,:)), axis('tight'), ylabel('F')
    subplot(312), stem(S(i).n,'Color',col(i,:)), hold on, plot(y(i,:),'Color',col(i,:)), axis('tight'), ylabel('y')
    subplot(313), plot(p(i,2:end),'Color',col(i,:)), hold on, axis('tight'), ylabel('p')
end

%% 5) pop pf

F = zeros(Sim.Nc,Sim.T);
for i=1:Sim.Nc,F(i,:)=S(i).F;end

[I.M I.P]   = Pop_PF_v1_0(F,P,Sim);
