function PF_Pop_Inf_v2_5(T,pl)
clear all; close all; T = 10000; pl = 1;
% this script simulates Sim.Nc cells, and then infers spikes for each,
% assuming they are independent.
%
% 1) set simulation metadata
% 2) initialize parameters
% 3) simulate data
% 4) plot simulation (ie, truth)
% 5) estimate connection matrix from spikes
% 6) prep for spike inference and matrix learning
% 7) omega inference using fast filter
% 8) omega inference using particle filter
% 9) make pretty fig
%
% Remarks:
% a) some of the code is general for Sim.M spike history terms per neurons, but not all (eg, time constants)
% b) inference assumes all correct parameters (except those governing GLM)
% c) # of external stimulus dimensions is currently restricted to 1
%
% Version Updates:
% v2_1: put EM in loop, assume that each neuron has 1 spike history term in
% E-step, don't plot inference between each EM iteration
%
% v2_2: code is general enough for arbitrary # neurons.  also started
% writing code to estimate matrix using fast method
%
% v2_3: is now a function that takes input T for time. also only saves
% stuff needed for plotting, not other crap.
%
% v2_4: now takes another input, 'pl', which, if 1, plots
% 
% v2_5: renamed simulation structure 'D' from 'S', changed some stuff to
% reduce memory requirements

%% 1) set simulation metadata

% metaparameters to simulate data
Sim.T       = T;                                    % # of time steps
Sim.dt      = 1/60;                                 % time step size
Sim.D       = 1;                                    % # dimensions of external stimulus
Sim.x       = ones(Sim.D,Sim.T);                    % stimulus
Sim.Nc      = 5;                                    % # of cells

%% 2) initialize parameters

rate        = 5;                                    % expected spike rate (assuming no spike history terms and Sim.x=1)
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
P.gamma     = 1e-5;                                 % scaled variance
P.zeta      = 4*P.gamma;                            % constant variance
P.tau_h     = 0.05;                                  % time constant
P.sigma_h   = 0.01;                                 % stan dev of noise
% P.omega     = [.1 -3; -3 .1];                   % weights
% if Sim.Nc>2
%     omega = zeros(Sim.Nc);
%     omega(1:2,1:2)=P.omega;
%     P.omega = omega;
% end
w=1;
P.omega=diag(-2*ones(Sim.Nc,1));
P.omega(1,2)=w;
P.omega(end,end-1)=w/2;
for i=2:2:Sim.Nc-1
    P.omega(i,i-1)=w/2;
    P.omega(i,i+1)=w;
end
for i=3:2:Sim.Nc-1
    P.omega(i,i-1)=-w;
    P.omega(i,i+1)=-w*2;
end
% figure(3), clf, imagesc(P.omega), colormap(gray)

%% 3) simulate data

D(1).h     = zeros(1,Sim.T);
D(1).n     = zeros(1,Sim.T);
D(1).C     = zeros(1,Sim.T);
D(1).F     = zeros(1,Sim.T);
for ii=2:Sim.Nc; D(ii)=D(1); end

kx      = P.k'*Sim.x;                                   % external input to neuron
eps_c   = P.sigma_c*sqrt(Sim.dt)*randn(Sim.Nc,Sim.T);   % generate noise on calcium
U_sampl = rand(Sim.Nc,Sim.T);                           % generate random number to use for sampling
eps_h   = repmat(P.sigma_h*sqrt(Sim.dt),Sim.Nc,Sim.T).*randn(Sim.Nc,Sim.T); % generate noise on spike history
eps_F   = randn(Sim.Nc,Sim.T);                          % generate noise on fluorescence
p       = zeros(Sim.Nc,Sim.T);                          % prob of spiking for each cell at each tim
y       = zeros(Sim.Nc,Sim.T);                          % input to each cell at each time

for t=2:Sim.T                                           % update states
    for i=1:Sim.Nc                                      % loop over presynaptic cells
        D(i).h(t)   = (1-Sim.dt./P.tau_h).*D(i).h(t-1) + D(i).n(t-1) + eps_h(i,t); % update h terms
    end

    for i=1:Sim.Nc                                      % loop over presynaptic cells
        y(i,t)      = kx(t);                            % initialize input to cell
        for j=1:Sim.Nc                                  % loop of post-synaptic cells
            y(i,t)  = y(i,t)+P.omega(i,j)*D(j).h(t);    % generate operand for rate function
        end
        p(i,t)      = 1-exp(-exp(y(i,t))*Sim.dt);       % generate prob of spiking
        D(i).n(t)   = U_sampl(i,t)<p(i,t);              % sample from bernoulli with prob p_t
        D(i).C(t)   = (1-Sim.dt/P.tau_c)*D(i).C(t-1) +...
            (Sim.dt/P.tau_c)*P.C_0 + P.A*D(i).n(t) + eps_c(i,t); %update calcium
        s           = Hill_v1(P,D(i).C(t));             % compute saturated calcium
        D(i).F(t)   = (P.alpha*s+P.beta)+sqrt(P.gamma*s+P.zeta).*eps_F(i,t); % update fluorescence
        if D(i).F(t)<=0; D(i).F(t) = eps; end           % keep fluorescence positive
    end
end

%% 4) plot simulation results

if pl==1
    tt  = 500;
    col = [1 0 0; 0 0 1; 0 0.5 0; 1 0.5 0; 1 0 1];          % define colors for mean
    figure(1), clf
    for i=1%:Sim.Nc
        h1=subplot(311); hold on, plot(((D(i).F(tt:2*tt))./max(D(i).F(tt:2*tt)))+1,'Color',col(i,:)); stem(D(i).n(tt:2*tt),'Color',col(i,:)), axis('tight'), ylabel('F')
        h2=subplot(312); hold on, plot(y(i,(tt:2*tt)),'Color',col(i,:)), axis('tight'), ylabel('y') %stem(D(i).n,'Color',col(i,:)),
        h3=subplot(313); hold on, plot(p(i,(tt:2*tt)),'Color',col(i,:)), axis('tight'), ylabel('p')
    end
    linkaxes([h1 h2 h3],'x')
    legend('1','2')
    for i=1:Sim.Nc, disp(sum(D(i).n)); end
    %keyboard
end
%% 5) estimate connection matrix directly from spikes

Sim.M       = 1;                                    % # spike history terms per neuron (fixed at one for this version of code)
Sim.n_params= 1;                                    % if 1, estimate k
Sim.h_params= 1;                                    % if 1, estimate omega (self-coupling)
Sim.F_params= 0;                                    % if 1, estimate observation parameters
Sim.C_params= 0;                                    % whether to compute
Sim.StimDim = Sim.Nc;                               % set external stim dimesions to # cells
Tim         = Sim;                                  % Tim is Sim for this estimation
Tim.N       = 1;                                    % # of particles
for i=1:Sim.Nc
    h=zeros(Sim.Nc-1,Sim.T);
    Pre=1:Sim.Nc;                                   % generate list of presynaptic neurons
    Pre(Pre==i)=[];                                 % remove self
    k=0;                                            % counter of dimension
    for j=Pre                                       % loop thru all presynaptic neurons
        k=k+1;                                      % generate input to neuron based on posterior mean spike train from neuron j
        h(k,:) = D(j).h;
    end
    Tim.x       = [Sim.x; h];                       % append input from other neurons onto external stimulus
    E           = P;
    E.omega     = E.omega(i,i);                     % initialize self-coupling term
    E.k         = E.k*ones(Sim.Nc,1);               % initialize external stim and cross-coupling terms
    D(i).w_b    = ones(1,Sim.T);
    Enew2{i}    = GOOPSI_Mstep_v1_0(Tim,D(i),0,E,D(i).F);
end

Phat{1}.omega = GetMatrix(Sim.Nc,Enew2);

% Phat{1}.omega=zeros(Sim.Nc);
% for i=1:Sim.Nc
%     Phat{1}.omega(i,i)=Enew2{i}.omega;
%     Pre=1:Sim.Nc;                                   % generate list of presynaptic neurons
%     Pre(Pre==i)=[];                                 % remove self
%     k=0;                                            % counter of dimension
%     for j=Pre
%         k=k+1;                                      % generate input to neuron based on posterior mean spike train from neuron j
%         Phat{1}.omega(i,j)=Enew2{i}.k(k+1);
%     end
% end
% [P.omega; round(Phat{1}.omega*10)/10]
%
if pl==1
    figure(3), clf,
    clims(1)=min(min(P.omega(:)),min(Phat{1}.omega(:)));
    clims(2)=max(max(P.omega(:)),max(Phat{1}.omega(:)));
    subplot(121), imagesc(P.omega,clims), colormap(gray), %colorbar
    subplot(122), imagesc(Phat{1}.omega,clims), %colorbar
    Fs=1024; ts=0:1/Fs:1; sound(sin(2*pi*ts*200)),
    %keyboard
end

%% 6) pop pf preparation

% metaparameters necessary to run smc-em code
Sim.N       = 99;                                   % # of particles
Sim.Mstep   = 1;                                    % whether to estimate parameters
Sim.pf      = 1;                                    % if 1, then conditional sampler, if 0, then prior sampler
Sim.freq    = 1;                                    % # time steps per observation (d in BJ08)
Sim.T_o     = Sim.T/Sim.freq;                       % # of observable time steps
Sim.M       = 1;                                    % # of spike history terms per neuron
Sim.Mstep   = 0;                                    % whether to estimate parameters
Sim.ptiles  = 1;                                    % generate percentiles

% initialize stuff for each cell
Tim         = Sim;                                  % copy Sim structure for input to Mstep function
E           = P;
E.omega     = E.omega(i,i);                         % initialize self-coupling term
E.k         = E.k*ones(Sim.StimDim,1);              % initialize external stim and cross-coupling terms
for i=1:Sim.Nc,
    I{i}.M.nbar = zeros(1,Sim.T);                   % initialize spike trains
    I{i}.P      = E;                                % initialize parameters
end

%% 7) Pop inference using fast filter

% E.tau       = E.tau_c;                                  % calcium decay time constant (sec)
% E.sig       = E.sigma_c;                                    % noise on
% E.lam       = Sim.T/(rate*E.A)*Sim.dt;              % expected jump size ber time bin
% Sim.Plot    = 1;
% Sim.MaxIter = 10;                                   % max number of EM iterations
%
% for i=1:Sim.Nc,                                     % initialize estimate of spike trains
%     [Fn{i} FP{i}]       = FOOPSI_v1_9(D(i).F',E,Sim);
% end
%
% for i=1:Sim.Nc
%
%     % append external stimulus for neuron 'i' with spike histories from other cells
%     h = zeros(Sim.Nc-1,Sim.T);                  % we append this to x to generate input into neuron from other neurons
%     Pre=1:Sim.Nc;                               % generate list of presynaptic neurons
%     Pre(Pre==i)=[];                             % remove self
%     k=0;                                        % counter of dimension
%     for j=Pre                                   % loop thru all presynaptic neurons
%         k=k+1;                                  % generate input to neuron based on posterior mean spike train from neuron j
%         h(k,:) = filter(1,[1 -(1-Sim.dt/P.tau_h)],Fn{j}/max(Fn{j}));
%     end
%     Tim.x   = [Sim.x; h];                         % append input from other neurons onto external stimulus
%     Tim.N   = 1;
%     T.w_b   = ones(1,Sim.T);
%     T.n     = Fn{i};
%     T.C     = filter(1,[1 -(1-Sim.dt/FP{i}.tau)],Fn{j}/max(Fn{j}));
%     T.h     = filter(1,[1 -(1-Sim.dt/E.tau_h)],Fn{j}/max(Fn{j}));;
%     FP{i}   = GOOPSI_Mstep_v1_0(Tim,T,0,II{i}.P,D(i).F);
% end
% PlotPop_v1_0(Sim,FE,P.omega,Phat{1}.omega)

%% 8)  Pop EM iteration

for tr=1:1
    fprintf('\nTrial # %g\n',tr)

    % for each neuron, infer spike train conditioned on previous EM
    % iterations spike history terms
    for i=1:Sim.Nc,                                 % infer spikes for each neuron

        % append external stimulus for neuron 'i' with spike histories from other cells
        h = zeros(Sim.Nc-1,Sim.T);                  % we append this to x to generate input into neuron from other neurons
        Pre=1:Sim.Nc;                               % generate list of presynaptic neurons
        Pre(Pre==i)=[];                             % remove self
        k=0;                                        % counter of dimension
        for j=Pre                                   % loop thru all presynaptic neurons
            k=k+1;                                  % generate input to neuron based on posterior mean spike train from neuron j
            h(k,:) = filter(1,[1 -(1-Sim.dt/P.tau_h)],I{j}.M.nbar);
        end
        Tim.x = [Sim.x; h];                         % append input from other neurons onto external stimulus

        % infer spike train for neuron 'i'
        fprintf('\nNeuron # %g\n',i)
        [S M II{i}.P]   = GOOPSI_main_v3_0(D(i).F,I{i}.P,Tim);
        II{i}.S.n       = S.n; 
        II{i}.S.h       = S.h; 
        II{i}.S.w_b     = S.w_b; 
        II{i}.M.nbar    = M.nbar;
    end

    % set inference for each neuron to the newly updated inference
    for i=1:Sim.Nc, I{i}.S = II{i}.S; I{i}.M = II{i}.M; end

    % given new inference for each neuron, update parameters
    for i=1:Sim.Nc

        % append external stimulus for neuron 'i' with spike histories from other cells
        h = zeros(Sim.Nc-1,Sim.T);                  % we append this to x to generate input into neuron from other neurons
        Pre=1:Sim.Nc;                               % generate list of presynaptic neurons
        Pre(Pre==i)=[];                             % remove self
        k=0;                                        % counter of dimension
        for j=Pre                                   % loop thru all presynaptic neurons
            k=k+1;                                  % generate input to neuron based on posterior mean spike train from neuron j
            h(k,:) = filter(1,[1 -(1-Sim.dt/P.tau_h)],I{j}.M.nbar);
        end
        Tim.x = [Sim.x; h];                         % append input from other neurons onto external stimulus

        fprintf('\nNeuron # %g\n',i)
        I{i}.P  = GOOPSI_Mstep_v1_0(Tim,I{i}.S,0,II{i}.P,D(i).F);
        EE{i}    = I{i}.P;
    end
    Phat{2}.omega = GetMatrix(Sim.Nc,EE);
    PlotMatrix(D,Sim,P,Phat)
    save(['SimConnector', num2str(Sim.T), '_', num2str(Sim.Nc), '_', num2str(tr)],'Phat','D','Sim','P')
    Fs=tr*1024; ts=0:1/Fs:1; sound(sin(2*pi*ts*200)),
end

end

function omega = GetMatrix(Nc,P)

omega=zeros(Nc);
for i=1:Nc
    omega(i,i)=P{i}.omega;
    Pre=1:Nc;                                       % generate list of presynaptic neurons
    Pre(Pre==i)=[];                                 % remove self
    k=0;                                            % counter of dimension
    for j=Pre
        k=k+1;                                      % generate input to neuron based on posterior mean spike train from neuron j
        omega(i,j)=P{i}.k(k+1);
    end
end

end

function PlotMatrix(D,Sim,P,Phat)

%%
fig=figure(5); clf,
col = [1 0 0; 0 0 1; 0 0.5 0; 1 0.5 0; 1 0 1];          % define colors for mean

% Phat.omega=zeros(Sim.Nc);
% for i=1:Sim.Nc
%     Phat.omega(i,i)=I{i}.P.omega;
%     Pre=1:Sim.Nc;                                   % generate list of presynaptic neurons
%     Pre(Pre==i)=[];                                 % remove self
%     k=0;                                            % counter of dimension
%     for j=Pre
%         k=k+1;                                      % generate input to neuron based on posterior mean spike train from neuron j
%         Phat.omega(i,j)=I{i}.P.k(k+1);
%     end
% end

tstart=2/Sim.dt;
tdt=1/Sim.dt;
tend = min(tstart+8*tdt,Sim.T);
subplot(2,3,[1 2 3]), hold on
for i=1:2
    plot(((D(i).F(tstart:tend))./(D(i).F(tstart:tend)))+1,'Color',col(i,:));
end
axis('tight'),
xticks  = tstart:tdt:tend;               % XTick positions
set(gca,'YTick',[],'YTickLabel',[])
set(gca,'XTick',xticks,'XTickLabel',xticks*Sim.dt)
ylabel('Fluorescence')
xlabel('Time (sec)')

clims(1)=min([P.omega(:)' Phat{1}.omega(:)' Phat{2}.omega(:)']);
clims(2)=max([P.omega(:)' Phat{1}.omega(:)' Phat{2}.omega(:)']);
subplot(234), imagesc(P.omega,clims), colormap(gray),
set(gca,'XTick',[1:Sim.Nc],'YTick',[1:Sim.Nc]) %colorbar
title('True matrix')
ylabel('Presynaptic'), xlabel('Postsynaptic')

subplot(235), imagesc(Phat{1}.omega,clims), %colorbar
set(gca,'XTick',[1:Sim.Nc],'YTick',[1:Sim.Nc]) %colorbar
title('Matrix from spikes')
% ylabel('Presynaptic'), xlabel('Postsynaptic')

subplot(236), imagesc(Phat{2}.omega,clims), %colorbar
set(gca,'XTick',[1:Sim.Nc],'YTick',[1:Sim.Nc]) %colorbar
title('Matrix from fluorescence')
% ylabel('Presynaptic'), xlabel('Postsynaptic')

% print fig
wh=[7 5];   %width and height
set(fig,'PaperPosition',[0 11-wh(2) wh]);
print('-depsc','SimConnector')
print('-dpdf','SimConnector')

end