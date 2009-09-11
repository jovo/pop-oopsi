%THIS FILE GENERATES figure;S FOR POP-INF PAPER WITH Liam AND Josh. 
%THIS FILE SHOULD BE PLACED IN THE SAME PLACE WITH THE DATA FILES.
close all

%% Some settings -- also make sure to fix throughout below after 'clear's
font_siz=14;                                %font size
R=(1-exp(-0.015/0.01))/(0.015/0.01);        %scale-correction thing

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% figure; 3) MCMC vs IID weights inferrence, N=25, T=600s
clear
font_siz=14;                                %font size
R=(1-exp(-0.015/0.01))/(0.015/0.01);        %scale-correction thing

load('fluor-0610-25-3-glm_iid0proc1.mat')
W_iid=Weights;
load('fluor-0610-25-3-glm_neal2proc1.mat')
W_mcmc=Weights;
load('fluor-0610-25-4-glm_base0proc1.mat')
W_base=Weights;
load('fluor-0610-25-3-result_base.mat')
W_GT=full(result_base.GT_W);

h=figure;                              % figure; 3a), scatter plot
plot(W_GT(:),W_iid(:),'.','Color',[0 0 0],'MarkerSize',12)
hh=get(h,'Children');set(hh,'FontSize',font_siz); set(hh,'FontName','Arial');
hold on
plot(W_GT(:),W_mcmc(:),'+','Color',[1 0 0],'MarkerSize',6)
plot(W_GT(:),W_base(:),'*','Color',[0 0 1])
plot([-2 2],[-2 2],'k--','LineWidth',2,'Color',[0.5 0.5 0.5])
axis([-2 1.5 -1.25 0.75])
xlabel('Actual connection weights')
ylabel('Inferred connection weights')
title('No scale correction')
legend({'Factorized approx.','Embedded-chain-Gibbs','Downsampled trains'},'Location','NorthWest')

h=figure;                              % figure; 3b), scatter plot bias corrected
plot(W_GT(:),W_iid(:)/R,'.','Color',[0 0 0],'MarkerSize',12)
hh=get(h,'Children');set(hh,'FontSize',font_siz); set(hh,'FontName','Arial');
hold on
plot(W_GT(:),W_mcmc(:)/R,'+','Color',[1 0 0],'MarkerSize',6)
r_base=regress(W_base(:),W_GT(:));  % find scaling bias, base
plot(W_GT(:),W_base(:)/R,'*','Color',[0 0 1])
plot([-2 2],[-2 2],'k--','LineWidth',2,'Color',[0.5 0.5 0.5])
axis([-2 2 -2 2])
xlabel('Actual connection weights')
ylabel('Inferred connection weights')
title('Scale correction')
legend({'Factorized approx.','Embedded-chain-Gibbs','Downsampled trains'},'Location','NorthWest')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% figure; 4) scaling bias as function of FR & r^2 as function of FR
fr=[15,33,66,100,200,1000];
dt=1./fr; theor=(1-exp(-dt/0.01))./(dt/0.01);
r=[]; c=[];
for k=fr
  fname=sprintf('frame-0717-%iFR-result_25_600.mat',k);
  load(fname);  
  idx=result.GT_W>0;
  [r1,b1]=regress(result.Weights_glm(idx),full(result.GT_W(idx)));
  r=[r,b1'];
  r1=corrcoef(result.Weights_glm(:),result.GT_W(:));
  c=[c,r1(3)^2];
end
h=figure;      %scaling bias
plot([dt*1000,0],[theor,1],'k-','LineWidth',2)
hh=get(h,'Children');set(hh,'FontSize',font_siz); set(hh,'FontName','Arial');
hold on
errorbar(dt*1000,mean(r,1),diff(r,1),'ko')
legend('theor.','actual')
xlabel('Time discretization, ms'),ylabel('Scaling factor')
%title('Scale bias for different time-bin size');
axis([0 70 -0.02 1.2])

h=figure;      %r^2
plot(fr(1:end-1),c(1:end-1),'k.-','LineWidth',2,'MarkerSize',font_siz);
hh=get(h,'Children');set(hh,'FontSize',font_siz); set(hh,'FontName','Arial');
xlabel('Frame rate, Hz'),ylabel('r^2')
%title('Impact of imaging frame-rates on inferrence');
axis([0 210 0 0.8])



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% figure; 5) r^2 vs SNR for N=25, 50, 100, THIS WILL BE MOD BY OTHER FR
gamma=[1e3,5e3,10e3,20e3,40e3,80e3];
gamma=[2.2,4.8,6.4,8.3,10.3,12]*1e3;
load('fluor-0610-25-1-result_base.mat')
W_GT=result_base.GT_W;
r=zeros(1,6);
for k=1:6
  load(sprintf('fluor-0610-25-%i-glm_iid0proc1.mat',k))
  W_iid=Weights;
  c=corrcoef(W_GT(:),W_iid(:));
  r(k)=c(3)^2;
end
h=figure;;
plot(gamma/1000,r,'k.-','LineWidth',2,'MarkerSize',16);
hh=get(h,'Children');set(hh,'FontSize',font_siz); set(hh,'FontName','Arial');
% Upper SNR bar is made manually, don't know how to do it in matlab
% commands; translation of photon budget to SNR:
% SNR: 4.8 5.9 6.7 7.5 7.8 8.1 8.4 8.7 dimless
% ph.: 10  20  30  40  50  60  70  80  Kph/frame/neuron

r=zeros(1,6);
for k=1:6
  load(sprintf('fluor-0717-25-1%i-result_bf.mat',k))
  W_GT=result_iid.GT_W;
  W_iid=result_iid.Weights_glm;
  c=corrcoef(W_GT(:),W_iid(:));
  r(k)=c(3)^2;
end
hold on
plot(gamma/1000,r,'k.--','LineWidth',2,'MarkerSize',16);

r=zeros(1,6);
for k=1:6
  load(sprintf('fluor-0717-25-2%i-result_bf.mat',k))
  W_GT=result_iid.GT_W;
  W_iid=result_iid.Weights_glm;
  c=corrcoef(W_GT(:),W_iid(:));
  r(k)=c(3)^2;
end
hold on
plot(gamma/1000,r,'k.:','LineWidth',2,'MarkerSize',16);

%bars
plot([3 3], [0 1],'k-','LineWidth',2)
plot([9 9],[0 1],'k-','LineWidth',2)

xlabel('eSNR')
ylabel('r^2')
legend({'66Hz','33Hz','15Hz'},'Location','NorthWest')
%title('Impact of imaging noise on inferrence')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% figure; 6) r^2 vs T for sparse prior for N=100, 200
%PREPROCESSING PART ++++++++++++++++++++++++++++++++++++++++++++
clear
font_siz=14;                                %font size
R=(1-exp(-0.015/0.01))/(0.015/0.01);        %scale-correction thing

N=[10,25,50,75,100,150,200,250];      %regular
rreg=zeros(4,8);
rspa=zeros(4,8);
rdal=zeros(4,8);
for k=[5,7]
  clear result_300 result_600 result_1800 result_3600
  load(sprintf('probe-0528-result_%i.mat',N(k)),'result_300');
  c=corrcoef(result_300.GT_W(:),result_300.Weights_glm(:));
  rreg(1,k)=c(3)^2;
  c=corrcoef(result_300.GT_W(:),result_300.Weights_spa(:));
  rspa(1,k)=c(3)^2;  
  c=corrcoef(result_300.GT_W(:),result_300.Weights_dal(:));
  rdal(1,k)=c(3)^2;    

  clear result_300
  load(sprintf('probe-0528-result_%i.mat',N(k)),'result_600');  
  c=corrcoef(result_600.GT_W(:),result_600.Weights_glm(:));
  rreg(2,k)=c(3)^2;
  c=corrcoef(result_600.GT_W(:),result_600.Weights_spa(:));
  rspa(2,k)=c(3)^2;  
  c=corrcoef(result_600.GT_W(:),result_600.Weights_dal(:));
  rdal(2,k)=c(3)^2;    
  
  clear result_600
  load(sprintf('probe-0528-result_%i.mat',N(k)),'result_1800');    
  c=corrcoef(result_1800.GT_W(:),result_1800.Weights_glm(:));
  rreg(3,k)=c(3)^2;
  c=corrcoef(result_1800.GT_W(:),result_1800.Weights_spa(:));
  rspa(3,k)=c(3)^2;  
  c=corrcoef(result_1800.GT_W(:),result_1800.Weights_dal(:));
  rdal(3,k)=c(3)^2;    
  
  clear result_1800
  load(sprintf('probe-0528-result_%i.mat',N(k)),'result_3600');    
  c=corrcoef(result_3600.GT_W(:),result_3600.Weights_glm(:));
  rreg(4,k)=c(3)^2;
  c=corrcoef(result_3600.GT_W(:),result_3600.Weights_spa(:));
  rspa(4,k)=c(3)^2;  
  c=corrcoef(result_3600.GT_W(:),result_3600.Weights_dal(:));
  rdal(4,k)=c(3)^2;      
end
clear result_300 result_600 result_1800 result_3600
%PREPROCESSING PART ++++++++++++++++++++++++++++++++++++++++++++

load('scan-0528-probe.mat')
T=zeros(size(result_50));
r_base=T;
r_spa=r_base;
r_dal=r_base;
h=figure;
axis([0 3600 0 1])
hh=get(h,'Children');set(hh,'FontSize',font_siz); set(hh,'FontName','Arial');

%N=100 neurons overlay
T=[];
r_base=[];
r_spa=[];
r_dal=[];
for k=300:300:2100
  sname=sprintf('scan-0528-result_%i_100.mat',k);
  load(sname);
  eval(sprintf('result_100=result_%i_100;',k));
  T(end+1)=result_100.T;
  c=corrcoef(result_100.GT_W(:),result_100.Weights_glm(:));
  r_base(end+1)=c(3)^2;
  c=corrcoef(result_100.GT_W(:),result_100.Weights_spa(:));
  r_spa(end+1)=c(3)^2;
  c=corrcoef(result_100.GT_W(:),result_100.Weights_dal(:));
  r_dal(end+1)=c(3)^2;
  clear(sprintf('result_%i_100',k));
end
%add missing T=3600 point from a different run
T=[T,3600]; 
r_base=[r_base,rreg(4,5)]; 
r_spa=[r_spa,rspa(4,5)];
r_dal=[r_dal,rdal(4,5)];
plot(T,r_base,'k-','LineWidth',2,'Color',[0.3,0.3,0.3])
hold on
plot(T,r_spa,'k--','LineWidth',2,'Color',[0.3,0.3,0.3])
plot(T,r_dal,'k:','LineWidth',2,'Color',[0.3,0.3,0.3])
axis([0 3600 0 1])

%N=200 neurons overlay
plot([300 600 1800 3600],rreg(:,7),'k-','LineWidth',2,'Color',[0.6,0.6,0.6]);
plot([300 600 1800 3600],rspa(:,7),'k--','LineWidth',2,'Color',[0.6,0.6,0.6]);
plot([300 600],rdal(1:2,7),'k:','LineWidth',2,'Color',[0.6,0.6,0.6]);
xlabel('Data available, sec')
ylabel('r^2')
legend({'Regular GLM','Sparse GLM','Sparse & Dale GLM'},'Location','SouthEast')
%title('Impact of imaging time on inferrence')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% figure; -) Histograms
load('probe-0528-result_200_600.mat')
W_GT=result_600.GT_W;
W_glm=result_600.Weights_glm;
W_spa=result_600.Weights_spa;

r=regress(W_glm(:),full(W_GT(:)));
hh=figure;;                           % figure; 3d), inferred distributions
plot([0 0],[0 1],'k--','LineWidth',3);        % zero branch
hh=get(hh,'Children');set(hh,'FontSize',font_siz); set(hh,'FontName','Arial');
hold on
[h,n]=hist(W_glm(W_GT==0)/r,13);
h=[0,h,0]; n=[2*n(1)-n(2),n,2*n(end)-n(end-1)];
plot(n,h/sum(h),'k','LineWidth',2);
legend({'True','Inferred'},'Location','NorthWest')

[h,n]=hist(W_GT(W_GT>0),13);         % original distributions of weights
h=[0,h,0]; n=[2*n(1)-n(2),n,2*n(end)-n(end-1)];
plot(n,h/sum(h),'r--','LineWidth',3);% positive branch
hold on
[h,n]=hist(W_GT(W_GT<0),13);
h=[0,h,0]; n=[2*n(1)-n(2),n,2*n(end)-n(end-1)];
plot(n,h/sum(h),'b--','LineWidth',3); % negative branch
[h,n]=hist(W_glm(W_GT>0)/r,13);       % inferred distributions mcmc
h=[0,h,0]; n=[2*n(1)-n(2),n,2*n(end)-n(end-1)];
plot(n,h/sum(h),'r','LineWidth',2);
hold on
[h,n]=hist(W_glm(W_GT<0)/r,13);
h=[0,h,0]; n=[2*n(1)-n(2),n,2*n(end)-n(end-1)];
plot(n,h/sum(h),'b','LineWidth',2);
xlabel('Connection weights')
ylabel('Histogram')
title('No prior')
axis([-6 2 0 1])

r=regress(W_spa(:),full(W_GT(:)));
hh=figure;;                           % figure; 3e), inferred distributions
[h,n]=hist(W_spa(W_GT>0)/r,13);       % inferred distributions baseline
h=[0,h,0]; n=[2*n(1)-n(2),n,2*n(end)-n(end-1)];
plot(n,h/sum(h),'r','LineWidth',2);
hh=get(hh,'Children');set(hh,'FontSize',font_siz); set(hh,'FontName','Arial');
hold on
[h,n]=hist(W_spa(W_GT<0)/r,13);
h=[0,h,0]; n=[2*n(1)-n(2),n,2*n(end)-n(end-1)];
plot(n,h/sum(h),'b','LineWidth',2);
hold on
[h,n]=hist(W_spa(W_GT==0)/r,99);
h=[0,h,0]; n=[2*n(1)-n(2),n,2*n(end)-n(end-1)];
plot(n,h/sum(h),'k','LineWidth',2);
legend({'Positive weights','Negative weights','Zero weights'},'Location','NorthWest')

plot([0 0],[0 1],'k--','LineWidth',3);        % zero branch
[h,n]=hist(W_GT(W_GT>0),13);          % original distributions of weights
h=[0,h,0]; n=[2*n(1)-n(2),n,2*n(end)-n(end-1)];
plot(n,h/sum(h),'r:','LineWidth',3); % positive branch
hold on
[h,n]=hist(W_GT(W_GT<0),13);
h=[0,h,0]; n=[2*n(1)-n(2),n,2*n(end)-n(end-1)];
plot(n,h/sum(h),'b--','LineWidth',3); % negative branch
xlabel('Connection weights')
ylabel('Histogram')
title('Sparse Prior')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% figure; 7) Strong/weak solution
h=figure;                          %figure; 10c) weak correlations
load fluor-0610-50-1-data.mat netSim
load fluor-0610-50-1-glm_base0proc1.mat
plot(netSim.weights(:),Weights(:)/R,'k.')
hh=get(h,'Children');set(hh,'FontSize',font_siz); set(hh,'FontName','Arial');
hold on
plot([-10 10],[-10 10],'k--','LineWidth',2,'Color',[0.5 0.5 0.5])
axis([-2 2 -2 2])
xlabel('Actual connection weights')
ylabel('Inferred connection weights')
% title('Weak correlations')

h=figure;                            %figure; 10d) strong correlations
load fluor-0706-50-9-result_bf.mat
plot(result_base.GT_W(:),result_base.Weights_glm(:)/R,'k.')
hh=get(h,'Children');set(hh,'FontSize',font_siz); set(hh,'FontName','Arial');
hold on
plot([-10 10],[-10 10],'k--','LineWidth',2,'Color',[0.5 0.5 0.5])
axis([-4 4 -4 4])
xlabel('Actual connection weights')
ylabel('Inferred connection weights')
% title('Strong correlations')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% figure; 8) variation in EPSP time constants, tau
clear
font_siz=14;                                %font size
R=(1-exp(-0.015/0.01))/(0.015/0.01);        %scale-correction thing

load vartau-0801-25-5-result_bf.mat
W_GT=full(result_iid.GT_W);
W_iid=result_iid.Weights_glm;

h=figure;;                             % figure; 8a), scatter plot
plot(W_GT(:),W_iid(:)/R,'k.')
hh=get(h,'Children');set(hh,'FontSize',font_siz); set(hh,'FontName','Arial');
hold on
plot([-2 2],[-2 2],'k--','LineWidth',2,'Color',[0.5 0.5 0.5])
axis([-2 2 -2 2])
xlabel('Actual connection weights')
ylabel('Inferred connection weights')
title('Constant  \tau_w')


load vartauv-0801-25-5-result_bf.mat
W_GT=full(result_iid.GT_W);
W_iid=result_iid.Weights_glm;

h=figure;;
plot(result_iid.GT_W(:),result_iid.Weights_glm(:)/R,'k.')
hh=get(h,'Children');set(hh,'FontSize',font_siz); set(hh,'FontName','Arial');
hold on
plot([-10 10],[-10 10],'k--','LineWidth',2,'Color',[0.5 0.5 0.5])
axis([-2 2 -2 2])
xlabel('Actual connection weights')
ylabel('Inferred connection weights')
title('Variable  \tau_w')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% figure; 9) regular glm & sparse overlaid for N=50, T=800
h=figure;;                    % figure; 3b, scatter plot bias corrected
load fluor-0610-50-5-data.mat netSim
Wt=full(netSim.weights);
load fluor-0610-50-5-glm_base0proc1.mat Weights
plot(Wt(:),Weights(:),'k.')
hh=get(h,'Children');set(hh,'FontSize',font_siz); set(hh,'FontName','Arial');
hold on
plot([-10 10],[-10 10],'k--','LineWidth',2,'Color',[0.5 0.5 0.5])
r=regress(Weights(:),Wt(:))
plot([-10 10],[-10 10]*r,'k-.','LineWidth',1,'Color',[0.5 0.5 0.5])
axis([-2 2 -2 2])
xlabel('Actual connection weights')
ylabel('Inferred connection weights')
title('No prior')
axis([-2 1.5 -1.2 0.75])

h=figure;;                     % figure; 3b, scatter plot bias corrected
load fluor-0610-50-5-data.mat netSim
Wt=full(netSim.weights);
hold on
load fluor-0610-50-5-spa_base0proc1.mat Weights
plot(Wt(:),Weights(:),'k.')
hh=get(h,'Children');set(hh,'FontSize',font_siz); set(hh,'FontName','Arial');
hold on
plot([-10 10],[-10 10],'k--','LineWidth',2,'Color',[0.5 0.5 0.5])
r=regress(Weights(:),Wt(:))
plot([-10 10],[-10 10]*r,'k-.','LineWidth',1,'Color',[0.5 0.5 0.5])
axis([-2 2 -2 2])
xlabel('Actual connection weights')
ylabel('Inferred connection weights')
title('Sparse Prior')
axis([-2 1.5 -1.2 0.75])
box on


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %figure; 10) Real data
% load invivo_cut.mat
% load invivo-0716-23-1-result_iid.mat
% 
% figure;                      % figure; 10a) traces
% t=(1:1500)/30;
% subplot(3,1,1)
% plot(t,im2double(F{1}(1:1500)))
% subplot(3,1,2)
% plot(t,im2double(F{15}(1:1500)))
% subplot(3,1,3)
% plot(t,im2double(F{23}(1:1500)))
% xlabel('Time, sec')
% 
% figure;                      % figure; 10b) spike raster
% n=zeros(23,500);
% for k=1:23 n(k,:)=mean(result_iid.n{k}(:,1:500),1); end
% imagesc(n),colormap gray
% 
% figure;                      % figure; 10c) X-correlation
% n=zeros(23,length(F{1}));
% for k=1:23 n(k,:)=mean(result_iid.n{k},1); end
% C=zeros(23,23);
% for k=1:23
%   for l=1:23
%     if(k==l) continue; end
%     c=corrcoef(n(k,2:end),n(l,1:end-1));
%     C(k,l)=c(3)^2;
%   end
% end
% imagesc(max(-2,min(2,C))),colorbar
% 
% figure;                      % figure; 10d) glm-solution
% imagesc(max(-2,min(2,result_iid.Weights_glm))),colorbar
% 
% figure;                      % figure; 10e) sparse-solution
% imagesc(max(-2,min(2,result_iid.Weights_spa))),colorbar