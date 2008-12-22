% this script simulates Sim.Nc cells, and then infers spikes for each,
% assuming they are independent.
%
% 1) set simulation metadata
% 2) initialize parameters
% 3) simulate data
% 4) plot simulation (ie, truth)
% 5) SMCMC
% 6) plot results (ie, spike inference)
% 7) estimate connection matrix from calcium
% 8) estimate connection matrix from spikes
% 9) plot omega and estimate
% 
% Remarks:
% a) if Nc > 2, code automatically sets all coupling terms for neurons j=3,4,... to zero
% b) some of the code is general for Sim.M spike history terms per neurons, but not all (eg, time constants)
% c) when we infer spikes, we assume no spike history terms (not even self, coupling terms).  
% d) thus, we add them in appropriately for each particle before estimating omega 
% e) inference assumes all correct parameters (except those governing GLM)
% f) you'll want to play with section 7).  not sure whether it is working. it does something.
% g) # of external stimulus dimensions is currently restricted to 1
% h) Enew.k and Enew.omega together comprise the estimates for k and omega

clear; clc

%% 1) set simulation metadata

% metaparameters to simulate data
Sim.T       = 500;                                  % # of time steps
Sim.dt      = 1/100;                                % time step size
Sim.D       = 1;                                    % # dimensions of external stimulus
Sim.x       = ones(Sim.D,Sim.T);                    % stimulus
Sim.Nc      = 4;                                    % # of cells

% metaparameters necessary to run smc-em code
Sim.N       = 20;                                   % # of particles
Sim.Mstep   = 0;                                    % whether to estimate parameters
Sim.M       = 0;                                    % # of spike history terms per neuron
Sim.pf      = 1;                                    % if 1, then conditional sampler, if 0, then prior sampler
Sim.freq    = 1;                                    % # time steps per observation (d in BJ08)
Sim.T_o     = Sim.T/Sim.freq;                       % # of observable time steps
Sim.C_params= 0;                                    % whether to compute suff stats for calcium parameters

%% 2) initialize parameters

rate        = 10;                                   % expected spike rate (assuming no spike history terms and Sim.x=1)
P.k         = log(-log(1-rate*Sim.dt)/Sim.dt);      % linear filter
P.k         = P.k*ones(Sim.D,1);                    % initialize k to the right number of dimensions
P.tau_c     = 0.5;                                  % calcium decay time constant (sec)
P.A         = 20;                                   % jumps size (\mu M)
P.C_0       = 20;                                   % baseline [Ca++] (\mu M)
P.C_init    = P.C_0;                                % initial [Ca++] (\mu M)
P.sigma_c   = 1;                                    % noise on
P.n         = 1.0;                                  % hill equation exponent
P.k_d       = 200;                                  % hill coefficient
P.alpha     = 1;                                    % F_max
P.beta      = 0;                                    % F_min
P.gamma     = 4e-5;                                 % scaled variance
P.zeta      = 4*P.gamma;                            % constant variance
P.tau_h     = 0.2;                                  % time constant
P.sigma_h   = 0.01;                                 % stan dev of noise
P.omega     = [.1 -.3; -.1 -1.0];                    % weights
if Sim.Nc>2
    omega = zeros(Sim.Nc);
    omega(1:2,1:2)=P.omega;
    P.omega = omega;
end

%% 3) simulate data

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

col = [1 0 0; 0 0 1; 0 0.5 0; 1 0.5 0];          % define colors for mean
figure(1), clf
for i=1:Sim.Nc
    subplot(311), plot(z1(S(i).F)+1,'Color',col(i,:)); hold on, stem(S(i).n,'Color',col(i,:)), axis('tight'), ylabel('F')
    subplot(312), stem(S(i).n,'Color',col(i,:)), hold on, plot(y(i,:),'Color',col(i,:)), axis('tight'), ylabel('y')
    subplot(313), plot(p(i,2:end),'Color',col(i,:)), hold on, axis('tight'), ylabel('p')
end

 %% 5) pop pf
% 
% % infer spikes for each neuron
% F = zeros(1,Sim.T);
% for i=1:Sim.Nc,
%     F=S(i).F;
%     [I{i}.S I{i}.M I{i}.P] = GOOPSI_main_v2_0(F,P,Sim);
% end
% 
% %% 6) plot inference
% 
% figure(2), clf, nrows=Sim.Nc;
% for i=1:Sim.Nc
%     subplot(nrows,1,i), hold on,
%     stem(S(i).n,'LineStyle','none','Color','k'),
%     stem(I{i}.M.nbar,'Marker','none','Color',col(1,:))
% %     axis('tight'), 
%     ylabel(num2str(i))
% end
% 
% %% 7) estimate GLM parameters
% 
% Sim.M       = 1;                                    % # spike history terms per neuron (fixed at one for this version of code)
% Sim.n_params= 1;                                    % if 1, estimate k
% Sim.h_params= 1;                                    % if 1, estimate omega (self-coupling)
% Sim.F_params= 0;                                    % if 1, estimate observation parameters
% Sim.StimDim = Sim.Nc;                               % set external stim dimesions to # cells
% Tim         = Sim;                                  % copy Sim structure for input to Mstep function
% P.g         = 1-Sim.dt/P.tau_h;                     % for brevity
% 
% for i=1:Sim.Nc;
%     I{i}.S.h = zeros(Sim.N,Sim.T,Sim.M);            % extize spike history terms
%     h = zeros(Sim.Nc-1,Sim.T);                      % we append this to x to generate input into neuron from other neurons
%     Pre=1:Sim.Nc;                                   % generate list of presynaptic neurons
%     Pre(Pre==i)=[];                                 % remove self                                
%     k=0;                                            % counter of dimension
%     for j=Pre                                       % loop thru all presynaptic neurons
%         k=k+1;                                      % generate input to neuron based on posterior mean spike train from neuron j
%         h(k,:) = filter(1,[1 -(1-Sim.dt/P.tau_h)],I{j}.M.nbar);
%     end
%     Tim.x = [Sim.x; h];                             % append input from other neurons onto external stimulus
%     for m=1:Sim.M                                   % loop thru number of spike history terms per neuron (note that code only works for Sim.M=1)                                   
%         for t=2:Sim.T                               % obtain spike history terms (this is necessary because above we inferred spike trains assuming no spike history terms, even self-coupling)
%             I{i}.S.h(:,t,m)=P.g(m)*I{i}.S.h(:,t-1,m)+I{i}.S.n(:,t-1);
%         end
%     end
%     if i==1
%         figure(2),
%         subplot(nrows,1,1), hold on, plot(S(1).h+1,'k'), plot(I{1}.S.h(1,:,1)+1)
%         for j=2:Sim.Nc
%             subplot(nrows,1,j), hold on, plot(S(j).h+1,'k'), plot(Tim.x(j,:)+1)
%         end
%     end
%     
%     E = I{i}.P;
%     E.omega = E.omega(i,i);                         % initialize self-coupling term
%     E.k     = E.k*ones(Sim.StimDim,1);              % initialize external stim and cross-coupling terms
%     Enew{i}  = GOOPSI_Mstep_v1_0(Tim,I{i}.S,I{i}.M,E,F);
% end

%% 8) estimate connection matrix directly from spikes

% J = I;
for i=1:Sim.Nc;
    J{i}.S.h = zeros(Sim.N,Sim.T,Sim.M);            % extize spike history terms
    h = zeros(Sim.Nc-1,Sim.T);                      % we append this to x to generate input into neuron from other neurons
    Pre=1:Sim.Nc;                                   % generate list of presynaptic neurons
    Pre(Pre==i)=[];                                 % remove self                                
    k=0;                                            % counter of dimension
    for j=Pre                                       % loop thru all presynaptic neurons
        k=k+1;                                      % generate input to neuron based on posterior mean spike train from neuron j
        h(k,:) = filter(1,[1 -(1-Sim.dt/P.tau_h)],S(j).n);
    end
    Tim.x = [Sim.x; h];                             % append input from other neurons onto external stimulus
    for m=1:Sim.M                                   % loop thru number of spike history terms per neuron (note that code only works for Sim.M=1)                                   
        for t=2:Sim.T                               % obtain spike history terms (this is necessary because above we inferred spike trains assuming no spike history terms, even self-coupling)
            J{i}.S.h(:,t,m)=P.g(m)*J{i}.S.h(:,t-1,m)+S(i).n(:,t-1);
        end
    end
    J{i}.S.w_b=1/Sim.N*ones(Sim.N,Sim.T);
    
    E = P;
    E.omega = E.omega(i,i);                         % initialize self-coupling term
    E.k     = E.k*ones(Sim.Nc,1);              % initialize external stim and cross-coupling terms
    Enew2{i}  = GOOPSI_Mstep_v1_0(Tim,J{i}.S,J{i}.M,E,F);
end


%% 9) plot omega
% figure(4), clf,
% Phat.omega=zeros(Sim.Nc);
% for i=1:Sim.Nc
%     Phat.omega(i,i)=Enew{i}.omega;
%     Pre=1:Sim.Nc;                                   % generate list of presynaptic neurons
%     Pre(Pre==i)=[];                                 % remove self                                
%     k=0;                                            % counter of dimension
%     for j=Pre
%         k=k+1;                                      % generate input to neuron based on posterior mean spike train from neuron j
%         Phat.omega(i,j)=Enew{i}.k(k);
%     end
% end
% 
% Phat2.omega=zeros(Sim.Nc);
% for i=1:Sim.Nc
%     Phat2.omega(i,i)=Enew2{i}.omega;
%     Pre=1:Sim.Nc;                                   % generate list of presynaptic neurons
%     Pre(Pre==i)=[];                                 % remove self                                
%     k=0;                                            % counter of dimension
%     for j=Pre
%         k=k+1;                                      % generate input to neuron based on posterior mean spike train from neuron j
%         Phat2.omega(i,j)=Enew2{i}.k(k);
%     end
% end
% 
% clims(1)=min(min(P.omega(:)),min(Phat.omega(:)));
% clims(2)=max(max(P.omega(:)),max(Phat.omega(:)));
% subplot(131), imagesc(P.omega,clims), colormap(gray), %colorbar
% subplot(132), imagesc(Phat.omega,clims), %colorbar
% subplot(133), imagesc(Phat2.omega,clims), %colorbar