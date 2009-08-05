% clear, clc
% load('~/Research/oopsi/meta-oopsi/data/tom/20081126_13_05_43_orientation_Bruno_reg1_ori1_135umdepth.mat')
% %%
% nrows=2;
% ncols=3;
% fs=12;
% figure(1), clf
% 
% %%
% a=Im.MeanFrame;
% MeanFrame = (a-min(a(:)))/(max(a(:))-min(a(:)));
% roi_edges = Im.roi_edges;
% roi_edges(roi_edges>1)=1;
% 
% seg_frame = MeanFrame+roi_edges;
% subplot(nrows,ncols,[1 4])
% subimage(seg_frame)
% set(gca,'XTick',[0 30 60 90 120])
% colormap(gray)
% title('Mean Frame','FontSize',fs)
% 
% tvec=1000:2800;
% xtick=[0:900:tvec(end)];
% fs=12;
% xticklabel=xtick/30;
% subplot(nrows,ncols,2)
% plot(z1(detrend(double(Im.F(3,tvec))))+1,'k')
% axis([0 max(tvec)-min(tvec) 0 2])
% title('High SNR','FontSize',fs)
% set(gca,'XTick',xtick,'XTickLabel',[])
% set(gca,'YTick',[])
% ylab=ylabel([{'$\mathbf{F}$'}; {''}; {''}; {'$\widehat{\mathbf{n}}$'}],'FontSize',fs,'Interpreter','latex');
% set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
% 
% subplot(nrows,ncols,5)
% plot(z1(detrend(double(Im.F(4,tvec))))+1,'k')
% axis([0 max(tvec)-min(tvec) 0 2])
% title('Low SNR','FontSize',fs)
% set(gca,'XTick',xtick,'XTickLabel',xticklabel,'FontSize',fs)
% set(gca,'YTick',[])
% ylab=ylabel([{'$\mathbf{F}$'}; {''}; {''}; {'$\widehat{\mathbf{n}}$'}],'FontSize',fs,'Interpreter','latex');
% set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')


%%
netSim.N=100;                    %number of cells
netSim.SP=0.10;                  %sparseness of connections
netSim.scaleE=0.025;            %excitatiry abs scale
netSim.scaleI=0.1;              %inhibitory abs scale
netSim.balance=0.1;             %fraction of inhibitory neurons
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
P.omega=60*(P.omega-min(P.omega(:)))/(max(P.omega(:))-min(P.omega(:)));
% P.omega=z1(P.omega);


%%
figure(1), clf
cmap=zeros(64,3);
zero_point=33; 
cmap(:,1)=[zeros(zero_point,1)' linspace(0,1,64-zero_point)];
cmap(:,3)=[linspace(1,0,64-zero_point) zeros(1,zero_point)];

% subplot(nrows,ncols,[3 6])
subimage(P.omega,cmap), 
set(gca,'FontSize',fs)
colormap(cmap),
colorbar(...
    'location','eastoutside',...
    'XTick',linspace(1,60,3),...
    'XTickLabel',[-1 0 1])
% set(gca,'XTick',[])
title('Weight Matrix')


%%
wh = [7 4];
set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
FigName = ['~/Research/oopsi/pop-oopsi/figs/W_matrix'];
print('-depsc',FigName)
print('-dpdf',FigName)
