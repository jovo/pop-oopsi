% function weight_inference_example(T,pl)
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
% Input:
%   T:  # of time steps (min of 1000, when plotting)
%   pl: whether to plot stuff between iterations
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
Sim.T       = 10000;                                    % # of time steps
Sim.dt      = 1/60;                                 % time step size
Sim.D       = 1;                                    % # dimensions of external stimulus
Sim.x       = ones(Sim.D,Sim.T);                    % stimulus
Sim.Nc      = 10;                                    % # of cells
pl          = 1;

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

%%
% w=1;
% P.omega=diag(-2*ones(Sim.Nc,1));
% P.omega(1,2)=w;
% P.omega(end,end-1)=w/2;
% for i=2:2:Sim.Nc-1
%     P.omega(i,i-1)=w/2;
%     P.omega(i,i+1)=w;
% end
% for i=3:2:Sim.Nc-1
%     P.omega(i,i-1)=-w;
%     P.omega(i,i+1)=-w*2;
% end
%%
netSim.N=Sim.Nc;                    %number of cells
netSim.SP=0.20;                  %sparseness of connections
netSim.scaleE=0.025;            %excitatiry abs scale
netSim.scaleI=0.1;              %inhibitory abs scale
netSim.balance=0.2;             %fraction of inhibitory neurons
netSim.T=6e5;                   %simulation duration, steps
netSim.dt=0.001;                %simulation time-step
netSim.FR=0.015;                %sampling time step
netSim.transf=inline('exp(x)','x');   %nonlinearity
netSim.invtransf=inline('log(x)','x');%nonlinear inverse

%FEEDBACK PROFILES/should be same length
t=(100:-1:0);                   %profile length             
h=exp(-t/10)-exp(-t/1); 
netSim.hepsp=h/max(h);          %EPSP profile
h=exp(-t/20)-exp(-t/1); 
netSim.hipsp=h/max(h);          %IPSP profile
h=exp(-t/10);
netSim.hrefr=h/max(h);          %REFR profile

netSim.A=[20 40 12];            %Ca-spike intensity, min, mean & STD
netSim.sigma=[0.05 0.35 0.1];   %Ca-noise/spike ratio, min, mean & STD
netSim.tau=[0.15 0.25 0.05];    %Ca-decay constant, min, mean & STD
netSim.C_0=[0.05 0.30 0.1];     %Ca-baseline/spike ratio, min, mean & STD
netSim.n=1;                     %Hill-profile exponent
netSim.k_d=100;                 %Hill-profile scale

netSim.alpha=[0 1 0];           %fluorescence scaling, min, mean & STD
netSim.beta =[0 0 0];           %fluorescence offset, min mean & STD
netSim.gamma=4e4*[0.3 1 0.5];   %fluorescence generation, shot noise, min, mean & STD
netSim.zeta=[0 4 0];            %fluorescence generation, background noise
netSim.hill=inline('C.^(P.n)./(C.^(P.n)+P.k_d)','P','C');%nonlinearity

netSim.rate =[0 5 0];           %spontaneous ignition rate, min, mean & STD
netSim.irate=[0 5 0];           %spontaneous ignition rate for inh neurons, same
netSim.omega=[0 100 0];         %refraction strength, min, mean & STD

netSim.K=netSim.FR/netSim.dt;   %sampling/sim time-steps ratio

rate=netSim.rate(2);                %mean firing rate
identie=false(netSim.N,1);          %I/E identity
Ninh=round(netSim.N*netSim.balance);%number of inhibitory neurons

netSim.weights=full(sprand(netSim.N,netSim.N*(1-netSim.SP),netSim.SP));
netSim.weights=[netSim.weights rand(netSim.N,netSim.N*(netSim.SP))];

%excitatory neurons -- log-rates distribution, VogMish setup
Zmin=1e-6;                          %truncation of weights
idx=1:netSim.N-Ninh;
identie(idx)=1;
W=netSim.weights(:,idx); ind=find(W>0);
z=netSim.scaleE*exprnd(1,length(ind),1);
z=max(Zmin,exp(-rate*netSim.dt)-z);
z=log(-log(z)/netSim.dt/rate);
W(ind)=z; netSim.weights(:,idx)=W;

%inhibitory neurons -- log-rates distribution
Zmin=1-1e-6;                       %truncation of weights
idx=netSim.N-Ninh+1:netSim.N;
identie(idx)=0;
W=netSim.weights(:,idx); ind=find(W>0);
z=netSim.scaleI*exprnd(1,length(ind),1);
z=min(Zmin,exp(-rate*netSim.dt)-z);
% z=max(1-Zmin,exp(-rate*netSim.dt)+z);%this is correct form, but OK
z=-log(-log(z)/netSim.dt/rate);
W(ind)=z; netSim.weights(:,idx)=W; 

%remove self-couplings
for k=1:netSim.N netSim.weights(k,k)=0; end
P.omega=netSim.weights;
% P.omega(:,idx)=P.omega(:,idx)/max(max(abs(P.omega(:,idx))))*max(P.omega(:));
%%
cmap=zeros(64,3);
zero_point=floor((max(P.omega(:))-min(P.omega(:)))/(max(P.omega(:))-min(P.omega(:))*2)*64);
cmap(:,1)=[zeros(zero_point,1)' linspace(0,1,64-zero_point)];
cmap(:,3)=[linspace(1,0,64-zero_point) zeros(1,zero_point)];

figure(2), clf, imagesc(P.omega), colormap(cmap), colorbar

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
%         D(i).C(t)   = (1-Sim.dt/P.tau_c)*D(i).C(t-1) +...
%             (Sim.dt/P.tau_c)*P.C_0 + P.A*D(i).n(t) + eps_c(i,t); %update calcium
%         s           = Hill_v1(P,D(i).C(t));             % compute saturated calcium
%         D(i).F(t)   = (P.alpha*s+P.beta)+sqrt(P.gamma*s+P.zeta).*eps_F(i,t); % update fluorescence
%         if D(i).F(t)<=0; D(i).F(t) = eps; end           % keep fluorescence positive
    end
end

%% 4) plot simulation results

if pl==1
    tt  = 500;
    col = [1 0 0; 0 0 1; 0 0.5 0; 1 0.5 0; 1 0 1];          % define colors for mean
    figure(1), clf
    nnum=zeros(Sim.Nc,1); for i=1:Sim.Nc, nnum(i)=sum(D(i).n); end; [foo ind]=max(nnum);
    h1=subplot(311); hold on, plot(z1(D(ind).F(tt:2*tt))+1,'Color',col(1,:)); stem(D(ind).n(tt:2*tt),'Color',col(1,:)), axis('tight'), ylabel('F')
    h2=subplot(312); hold on, plot(y(ind,(tt:2*tt)),'Color',col(1,:)), axis('tight'), ylabel('y') %stem(D(i).n,'Color',col(i,:)),
    h3=subplot(313); hold on, plot(p(ind,(tt:2*tt)),'Color',col(1,:)), axis('tight'), ylabel('p')
    linkaxes([h1 h2 h3],'x')
    legend('1','2')
    for i=1:Sim.Nc, disp(sum(D(i).n)); end
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
    display(i)
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

%%
omega=zeros(Sim.Nc);
for i=1:Sim.Nc
    omega(i,i)=Enew2{i}.omega;
    Pre=1:Sim.Nc;                                       % generate list of presynaptic neurons
    Pre(Pre==i)=[];                                 % remove self
    k=0;                                            % counter of dimension
    for j=Pre
        k=k+1;                                      % generate input to neuron based on posterior mean spike train from neuron j
        omega(i,j)=Enew2{i}.k(k+1);
    end
end


Phat{1}.omega = omega;

save('~/Research/oopsi/pop-oopsi/figs/weight_inference_example','Phat','P')
%%

fs=12;

% make true W matrix fig
figure(3), clf,
cmap=zeros(64,3);
zero_point=ceil((max(P.omega(:))-min(P.omega(:)))/(max(P.omega(:))-min(P.omega(:))*2)*64);
cmap(:,1)=[zeros(zero_point,1)' linspace(0,1,64-zero_point)];
cmap(:,3)=[linspace(1,0,64-zero_point) zeros(1,zero_point)];

imagesc(P.omega), 
set(gca,'FontSize',fs)
colormap(cmap),
colorbar(...
    'location','southoutside',...
    'XTick',linspace(min(P.omega(:)),max(P.omega(:)),3),...
    'XTickLabel',[-1 0 1])
set(gca,'XAxisLocation','top')
title('True W')

wh = [4 4];
set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
FigName = ['~/Research/oopsi/pop-oopsi/figs/W_star'];
print('-depsc',FigName)
print('-dpdf',FigName)

% make estimated W matrix fig
cmap2=zeros(64,3);
zero_point2=ceil((max(Phat{1}.omega(:))-min(Phat{1}.omega(:)))/(max(Phat{1}.omega(:))-min(Phat{1}.omega(:))*2)*64);
cmap2(:,1)=[zeros(zero_point2,1)' linspace(0,1,64-zero_point2)];
cmap2(:,3)=[linspace(1,0,64-zero_point2) zeros(1,zero_point2)];

figure(4), 
imagesc(Phat{1}.omega), colormap(cmap2)
colorbar(...
    'location','southoutside',...
    'XTick',linspace(min(Phat{1}.omega(:)),max(Phat{1}.omega(:)),3),...
    'XTickLabel',[-1 0 1])
set(gca,'XAxisLocation','top')
set(gca,'YAxisLocation','right')
title('Estimated W')

wh = [4 4];
set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
FigName = ['~/Research/oopsi/pop-oopsi/figs/W_hat'];
print('-depsc',FigName)
print('-dpdf',FigName)