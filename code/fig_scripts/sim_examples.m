% this script does makes a figure showing how the smc-em approach
% outperforms other approaches even for just noisy data, by doing:
% 
% 1) set simulation metadata (eg, dt, T, # particles, etc.)
% 2) initialize parameters
% 3) generate fake data
% 4) infers spikes using a variety of approaches
% 5) plots results

clear, clc, fprintf('\nNoisy Simulation Fig\n')

%% 1) set simulation metadata

Sim.T       = 405;                                  % # of time steps
Sim.dt      = 1/60;                                % time step size
Sim.freq    = 1;                                    % # of time steps between observations
Sim.Nsec    = Sim.T*Sim.dt;                         % # of actual seconds
Sim.T_o     = Sim.T;                                % # of observations
Sim.tvec    = Sim.dt:Sim.dt:Sim.Nsec;               % vector of times
Sim.N       = 200;                                  % # of particles
Sim.M       = 0;                                    % number of spike history dimensions
Sim.pf      = 1;                                    % use conditional sampler (not prior) when possible
Sim.StimDim = 1;                                    % # of stimulus dimensions
Sim.x       = ones(1,Sim.T);                        % stimulus

Sim.Mstep   = false;                                % do M-step
Sim.C_params = true;                                 % whether to estimate calcium parameters {tau,A,C_0,sig}
Sim.n_params = true;                                 % whether to estimate rate governing parameters {b,k}
Sim.h_params = false;                                % whether to estimate spike history parameters {h}
Sim.F_params = false;                                % whether to estimate observation parameters {alpha,beta,gamma,zeta}
Sim.MaxIter = 0;                                    % max # of EM iterartions

%% 2) initialize parameters

% initialize barrier and wiener filter parameters

spt     = [47   105   117   162   267   332   376];
spt     = [spt 88 200 290 300 360 400];
P.rate  = numel(spt)/(Sim.T*Sim.dt);                % expected spike rate
P.A     = 1;                                        % jump size ($\mu$M)
P.tau   = 0.5;                                        % calcium decay time constant (sec)
P.lam   = 1/sqrt(P.rate*P.A);                           % expected jump size ber time bin
P.sig   = 1;                                        % standard deviation of noise (\mu M)

% initialize particle filter parameters
P.k         = log(-log(1-4*Sim.dt)/Sim.dt);    % linear filter
P.tau_c     = P.tau;
P.A         = 50;
P.C_0       = 1;                                   % baseline [Ca++]
P.C_init    = P.A;                                % initial [Ca++]
P.sigma_c   = 0.1;
P.n         = 1.0;                                  % hill equation exponent
P.k_d       = 200;                                  % hill coefficient
P.alpha     = 2;                                    % F_max
P.beta      = 0;                                    % F_min
P.gamma     = .002;                                 % scaled variance
P.zeta      = 0;                            % constant variance
P.a         = Sim.dt/P.tau_c;

%% 3) simulate data

% n=double(rand(Sim.T,1)<1-exp(-exp(P.k)*Sim.dt));
n=zeros(Sim.T,1); n(spt)=1;
C=zeros(Sim.T,1);
C(1)=P.C_init;
epsilon_c = P.sigma_c*sqrt(Sim.dt)*randn(1,Sim.T);      %generate noise on calcium
for t=2:Sim.T                                           %recursively update calcium
    C(t)  = (1-P.a)*C(t-1) + P.a*P.C_0 + P.A*n(t) + epsilon_c(t);   
end
R.C=C;
S=Hill_v1(P,R.C);
eps_t=randn(Sim.T,1);
F=P.alpha*S+P.beta+sqrt(P.gamma*S+P.zeta).*eps_t;

spt=find(n==1);

%%
diffF=diff(F);
dF1=diffF(spt-1);

nospt=1:Sim.T;
nospt(spt)=[];
dF0=diffF(nospt(2:end)-1);

mse1=sum(dF1.^1)/numel(dF1);
mse0=sum(dF0.^2)/numel(dF0);

eSNR=mse1/sqrt(mse0);

%% 4) plot black and white results
noisyF=P.alpha*S+P.beta+sqrt(P.gamma*8*S+P.zeta).*eps_t;
diffF=diff(noisyF);
dF1=diffF(spt-1);

nospt=1:Sim.T;
nospt(spt)=[];
dF0=diffF(nospt(2:end)-1);

mse1=sum(dF1.^1)/numel(dF1);
mse0=sum(dF0.^2)/numel(dF0);

eSNR2=mse1/sqrt(mse0);

% eSNR2=round(eSNR2*100)/100;
%%
% load('../data/NoisySim.mat')

fig=figure(1); clf,
nrows = 1;
ncols = 2;
gray  = [.5 .5 .5];                 % define gray color
col   = [1 0 0; 0.2 0.2 1];         % define colors for mean
ccol  = col+.4; ccol(ccol>1)=1;     % define colors for std
inter = 'none';                     % interpreter for axis labels
xlims = [Sim.T-200 Sim.T]-50;                % xmin and xmax for current plot
fs=12;                              % font size
ms=4;                               % marker size for real spike
sw=1.5;                             % spike width
lw=1;                               % line width
% make xticks
Nsec = floor(Sim.T*Sim.dt);
secs = zeros(1,Nsec-1);
for i=1:Nsec
    secs(i) = find(Sim.tvec>=i,1);
end
xticks=xlims(1)+[0:60:Sim.T];
I{2}.name=[{'Wiener Filter'}];
i=0;

col   = [1 0 0; 0.2 0.2 1];         % define colors for mean
ccol  = col+.4; ccol(ccol>1)=1;     % define colors for std
I{7}.name=[{'Linear Observation'}; {'PFS Spike Inference'}];

i=0;
% plot fluorescence data
i=i+1; subplot(nrows,ncols,i), hold on
plot(z1(F)+1/4,'k','LineWidth',lw);
stem(n/4,'Marker','none','MarkerSize',ms,'LineWidth',sw,'Color','k','MarkerFaceColor',gray,'MarkerEdgeColor',gray)

% stem(Sim.n,'Marker','none','LineWidth',sw,'Color','k')
ylab=ylabel([{'F'}; {''}; {''}; {''};  {'n'}],'Interpreter',inter,'FontSize',fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','top','Position',[140 .95 17])
set(gca,'YTick',1,'YTickLabel',[])
set(gca,'XTick',secs,'XTickLabel',[],'FontSize',fs)
axis([xlims 0 1.25])
% title(['eSNR=', num2str(eSNR)],'fontsize',fs)
title('Low Noise','fontsize',fs)

% % plot spike train
% subplot(nrows,ncols,i+2), hold on
% stem(n,'Marker','none','MarkerSize',ms,'LineWidth',sw,'Color','k','MarkerFaceColor',gray,'MarkerEdgeColor',gray)
% ylab=ylabel([{'n'}],'Interpreter',inter,'FontSize',fs);
% set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
% set(gca,'YTick',1,'YTickLabel',[])
% set(gca,'XTick',secs,'XTickLabel',[],'FontSize',fs)
% axis([xlims 0 1])

i=i+1; subplot(nrows,ncols,i), hold on
plot(z1(noisyF)+1/4,'k','LineWidth',lw);
stem(n/4,'Marker','none','MarkerSize',ms,'LineWidth',sw,'Color','k','MarkerFaceColor',gray,'MarkerEdgeColor',gray)

% stem(Sim.n,'Marker','none','LineWidth',sw,'Color','k')
% ylab=ylabel([{'Simulated'}; {'Fluorescence'}],'Interpreter',inter,'FontSize',fs);
% set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'YTick',1,'YTickLabel',[])
set(gca,'XTick',secs,'XTickLabel',[],'FontSize',fs)
axis([xlims 0 1.25])
% title(['eSNR=', num2str(eSNR2)],'fontsize',fs)
title('High Noise','fontsize',fs)

% % plot spike train
% subplot(nrows,ncols,i+2), hold on
% stem(n,'Marker','none','MarkerSize',ms,'LineWidth',sw,'Color','k','MarkerFaceColor',gray,'MarkerEdgeColor',gray)
% % ylab=ylabel([{'Simulated'}; {'Spike Train'}],'Interpreter',inter,'FontSize',fs);
% % set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
% set(gca,'YTick',1,'YTickLabel',[])
% set(gca,'XTick',secs,'XTickLabel',[],'FontSize',fs)
% axis([xlims 0 1])

minsec=0;%floor(xlims*Sim.dt);
subplot(nrows,ncols,1)
set(gca,'XTick',xticks,'XTickLabel',(xticks-xlims(1))*Sim.dt,'FontSize',fs)
xlabel('Time (sec)','FontSize',fs)

subplot(nrows,ncols,nrows+1)
set(gca,'XTick',xticks,'XTickLabel',(xticks-xlims(1))*Sim.dt,'FontSize',fs)
xlabel('Time (sec)','FontSize',fs)

% print fig
wh=[7 2];   %width and height
% set(fig,'PaperPosition',[0 11-wh(2) wh]);
set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
print('-depsc','../figs/sim_examples')
print('-dpdf','../figs/sim_examples')

[eSNR2 eSNR]
