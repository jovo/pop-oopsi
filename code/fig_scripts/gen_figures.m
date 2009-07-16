%THIS FILE GENERATES FIGURES FOR POP-INF PAPER WITH Liam AND Josh. 
%THIS FILE SHOULD BE PLACED IN THE SAME PLACE WITH THE DATA FILES.
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure 2) MCMC vs IID weights inferrence, N=25, T=600s
clear
load('fluor-0610-25-3-glm_iid0proc1.mat')
W_iid=Weights;
load('fluor-0610-25-3-glm_neal2proc1.mat')
W_mcmc=Weights;
load('fluor-0610-25-3-result_base.mat')
W_GT=full(result_base.GT_W);

figure                      % Figure 1a), scatter plot
plot(W_GT(:),W_iid(:),'r*')
hold on
plot(W_GT(:),W_mcmc(:),'bd','Color',[0.5 0 0.25],'MarkerSize',6)
plot([-2 2],[-2 2],'k--','LineWidth',2,'Color',[0.5 0.5 0.5])
axis([-2 1.5 -1.25 0.75])
xlabel('Actual connection weights')
ylabel('Inferred connection weights')
title('Figure 2a: MCMC-Gibbs vs Indep. approx. scatter plot')
legend({'Indep. approx.','MCMC-Gibbs'},'Location','NorthWest')

figure                      % Figure 1b), scatter plot bias corrected
r=regress(W_iid(:),W_GT(:));% find scaling bias, iid
plot(W_GT(:),W_iid(:)/r,'r*')
hold on
r=regress(W_iid(:),W_GT(:));% find scaling bias, mcmc
plot(W_GT(:),W_mcmc(:)/r,'bd','Color',[0.5 0 0.25],'MarkerSize',6)
plot([-2 2],[-2 2],'k--','LineWidth',2,'Color',[0.5 0.5 0.5])
axis([-2 2 -2 2])
xlabel('Actual connection weights')
ylabel('Inferred connection weights')
title('Figure 2b: MCMC-Gibbs vs Indep. approx. scatter plot bias corrected')
legend({'Indep. approx.','MCMC-Gibbs'},'Location','NorthWest')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure 3) MCMC vs IID weights inferrence, N=25 T=600s
clear
load('fluor-0610-25-4-glm_iid0proc1.mat')
W_iid=Weights;
load('fluor-0610-25-4-glm_base0proc1.mat')
W_base=Weights;
load('fluor-0610-25-4-result_base.mat')
W_GT=full(result_base.GT_W);

figure                  % Figure 2a, scatter plot
plot(W_GT(:),W_iid(:),'r*')
hold on
plot(W_GT(:),W_base(:),'.')
plot([-2 2],[-2 2],'k--','LineWidth',2,'Color',[0.5 0.5 0.5])
axis([-2 1.5 -1.25 0.75])
xlabel('Actual connection weights')
ylabel('Inferred connection weights')
title('Figure 3a: Indep. approx. vs Orig. spikes, scatter plot')
legend({'Indep. approx.','Orig. spikes'},'Location','NorthWest')

figure                  % Figure 2a, scatter plot bias corrected
r=regress(W_iid(:),W_GT(:));% find scaling bias, iid
plot(W_GT(:),W_iid(:)/r,'r*')
hold on
r=regress(W_base(:),W_GT(:));% find scaling bias, base
plot(W_GT(:),W_base(:)/r,'.')
plot([-2 2],[-2 2],'k--','LineWidth',2,'Color',[0.5 0.5 0.5])
axis([-2 2 -2 2])
xlabel('Actual connection weights')
ylabel('Inferred connection weights')
title('Figure 3b: Indep. approx. vs Orig. spikes, scatter plot bias corrected')
legend({'Indep. approx.','Orig. spikes'},'Location','NorthWest')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure 4) regular glm & sparse overlaid for N=50, T=800
figure                        % Figure 3a, scatter plot
load fluor-0610-50-5-data.mat netSim
Wt=full(netSim.weights);
load fluor-0610-50-5-glm_base0proc1.mat Weights
plot(Wt(:),Weights(:),'.')
hold on
load fluor-0610-50-5-spa_base0proc1.mat Weights
plot(Wt(:),Weights(:),'ro')
hold on
plot([-10 10],[-10 10],'k--','LineWidth',2,'Color',[0.5 0.5 0.5])
axis([-2 1.5 -1.25 0.75])
xlabel('Actual connection weights')
ylabel('Inferred connection weights')
title('Figure 4a: Regular GLM vs. sparse GLM for orig. spikes')
legend({'Regular GLM','Sparse GLM'},'Location','NorthWest')

figure                        % Figure 3b, scatter plot bias corrected
load fluor-0610-50-5-data.mat netSim
Wt=full(netSim.weights);
load fluor-0610-50-5-glm_base0proc1.mat Weights
r=regress(Weights(:),Wt(:));  % find scaling bias, base
plot(Wt(:),Weights(:)/r,'b.')
hold on
load fluor-0610-50-5-spa_base0proc1.mat Weights
r=regress(Weights(:),Wt(:));  % find scaling bias, sparse
plot(Wt(:),Weights(:)/r,'ro')
hold on
plot([-10 10],[-10 10],'k--','LineWidth',2,'Color',[0.5 0.5 0.5])
axis([-2 2 -2 2])
xlabel('Actual connection weights')
ylabel('Inferred connection weights')
title('Figure 4b: Regular GLM vs. sparse GLM for orig. spikes, bias corrected')
legend({'Regular GLM','Sparse GLM'},'Location','NorthWest')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure 5) Residuals for N=50, T=800s
clear
load('fluor-0610-50-5-glm_iid0proc1.mat')
W_glm=Weights;
load('fluor-0610-50-5-spa_base0proc1.mat')
W_spa=Weights;
load('fluor-0610-50-5-glm_base0proc1.mat')
W_base=Weights;
load('fluor-0610-50-1-data.mat','netSim')
W_GT=full(netSim.weights);

figure                              %figure of residuals iid
r=regress(W_glm(:),W_GT(:));
res=W_GT(:)-W_glm(:)/r;
[h,n]=hist(res(W_GT>0),9);         %original distributions of weights
plot(n,h/max(h),'r','LineWidth',2); %positive branch
hold on
[h,n]=hist(res(W_GT<0),9);
plot(n,h/max(h),'b','LineWidth',2); %negative branch
[h,n]=hist(res(W_GT==0),9);
plot(n,h/max(h),'k','LineWidth',2);
xlabel('Residual')
ylabel('Normalized histogram')
axis([-1 1 0 1])
title('Figure 5a: Residuals for indep. approx. inference')

figure                              %figure of residuals orig. spikes
r=regress(W_base(:),W_GT(:));
res=W_GT(:)-W_base(:)/r;
[h,n]=hist(res(W_GT>0),9);         %original distributions of weights
plot(n,h/max(h),'r','LineWidth',2); %positive branch
hold on
[h,n]=hist(res(W_GT<0),9);
plot(n,h/max(h),'b','LineWidth',2); %negative branch
[h,n]=hist(res(W_GT==0),9);
plot(n,h/max(h),'k','LineWidth',2);
xlabel('Residual')
ylabel('Normalized histogram')
axis([-1 1 0 1])
title('Figure 5b: Residuals for orig. spikes inference')

figure                              %figure of residuals orig. spikes, sparse pr
r=regress(W_spa(:),W_GT(:));
res=W_GT(:)-W_spa(:)/r;
[h,n]=hist(res(W_GT>0),9);         %original distributions of weights
plot(n,h/max(h),'r','LineWidth',2); %positive branch
hold on
[h,n]=hist(res(W_GT<0),9);
plot(n,h/max(h),'b','LineWidth',2); %negative branch
[h,n]=hist(res(W_GT==0),27);
plot(n,h/max(h),'k','LineWidth',2);
xlabel('Residual')
ylabel('Normalized histogram')
axis([-1 1 0 1])
title('Figure 5c: Residuals for orig. spikes inference, sparse prior')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure 6) histogram of weights
clear
load('fluor-0610-50-5-glm_iid0proc1.mat')
W_glm=Weights;
load('fluor-0610-50-5-spa_base0proc1.mat')
W_spa=Weights;
load('fluor-0610-50-1-data.mat','netSim')
W_GT=full(netSim.weights);

figure                              %figure of inferred distributions
r=regress(W_glm(:),W_GT(:));        %get scaling bias
[h,n]=hist(W_GT(W_GT>0),9);         %original distributions of weights
plot(n,h/sum(h),'r','LineWidth',2); %positive branch
hold on
[h,n]=hist(W_GT(W_GT<0),9);
plot(n,h/sum(h),'b','LineWidth',2); %negative branch
stem(0,1,'k','LineWidth',2);        %zero branch

[h,n]=hist(W_glm(W_GT>0)/r,9);        %inferred distributions regular GLM
plot(n,h/sum(h),'r--','LineWidth',2);
hold on
[h,n]=hist(W_glm(W_GT<0)/r,9);
plot(n,h/sum(h),'b--','LineWidth',2);
[h,n]=hist(W_glm(W_GT==0)/r,9);
plot(n,h/sum(h),'k--','LineWidth',2);
                                    %inferred distributions sparse GLM
r=regress(W_spa(:),W_GT(:));        %get scaling bias
[h,n]=hist(W_spa(W_GT>0)/r,9);
plot(n,h/sum(h),'r:','LineWidth',2);
[h,n]=hist(W_spa(W_GT<0)/r,9);
plot(n,h/sum(h),'b:','LineWidth',2);
[h,n]=hist(W_spa(W_GT==0)/r,27);
plot(n,h/sum(h),'k:','LineWidth',2);
xlabel('Connection weights')
ylabel('Histogram')
title('Figure 6: True and inferred distribution of connectivity weights')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure 7) r^2 vs SNR for N=25, 50, 100, THIS WILL BE MOD BY OTHER FR
gamma=[1e3,5e3,10e3,20e3,40e3,80e3];
load('fluor-0610-25-1-result_base.mat')
W_GT=result_base.GT_W;
r=zeros(1,6);
for k=1:6
  load(sprintf('fluor-0610-25-%i-glm_iid0proc1.mat',k))
  W_iid=Weights;
  c=corrcoef(W_GT(:),W_iid(:));
  r(k)=c(3)^2;
end
h=figure;
plot(gamma/1000,r,'k.-','LineWidth',2,'MarkerSize',16);
% Upper SNR bar is made manually, don't know how to do it in matlab
% commands; translation of photon budget to SNR:
% SNR: 4.8 5.9 6.7 7.5 7.8 8.1 8.4 8.7 dimless
% ph.: 10  20  30  40  50  60  70  80  Kph/frame/neuron


% % This part is doing N=50 neurons, drop in light of the mod above
% gamma=[1e3,5e3,10e3,20e3,40e3,80e3];
% load('fluor-0610-50-1-data.mat','netSim')
% W_GT=netSim.weights;
% r=zeros(1,6);
% for k=1:6
%   load(sprintf('fluor-0610-50-%i-glm_iid0proc1.mat',k))
%   W_iid=Weights;
%   c=corrcoef(W_GT(:),W_iid(:));
%   r(k)=c(3)^2;
% end
% figure(h)
% hold on
% plot(gamma/1000,r,'k.--','LineWidth',2);

xlabel('Photon budget, Kph/neuron/frame')
ylabel('Recovered variance')
legend({'N=25, T=600s'},'Location','NorthWest')
title('Figure 7: Impact of noise in calcium imaging data')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure 8) r^2 vs T for sparse prior for N=50,100, 200
%PREPROCESSING PART ++++++++++++++++++++++++++++++++++++++++++++
clear
N=[10,25,50,75,100,150,200,250];      %regular
r=zeros(4,8);
for k=[5,7]
  clear result_300 result_600 result_1800 result_3600
  load(sprintf('probe-0528-result_%i.mat',N(k)),'result_300');
  c=corrcoef(result_300.GT_W(:),result_300.Weights_glm(:));
  r(1,k)=c(3)^2;
  clear result_300
  load(sprintf('probe-0528-result_%i.mat',N(k)),'result_600');  
  c=corrcoef(result_600.GT_W(:),result_600.Weights_glm(:));
  r(2,k)=c(3)^2;
  clear result_600
  load(sprintf('probe-0528-result_%i.mat',N(k)),'result_1800');    
  c=corrcoef(result_1800.GT_W(:),result_1800.Weights_glm(:));
  r(3,k)=c(3)^2;
  clear result_1800
  load(sprintf('probe-0528-result_%i.mat',N(k)),'result_3600');    
  c=corrcoef(result_3600.GT_W(:),result_3600.Weights_glm(:));
  r(4,k)=c(3)^2;
end
rreg=r;
clear result_300 result_600 result_1800 result_3600

N=[10,25,50,75,100,150,200,250];        %sparse
r=zeros(4,8);
for k=[5,7]
  clear result_300 result_600 result_1800 result_3600
  load(sprintf('probe-0528-result_%i.mat',N(k)),'result_300');
  c=corrcoef(result_300.GT_W(:),result_300.Weights_spa(:));
  r(1,k)=c(3)^2;
  clear result_300
  load(sprintf('probe-0528-result_%i.mat',N(k)),'result_600');
  c=corrcoef(result_600.GT_W(:),result_600.Weights_spa(:));
  r(2,k)=c(3)^2;
  clear result_600
  load(sprintf('probe-0528-result_%i.mat',N(k)),'result_1800');
  c=corrcoef(result_1800.GT_W(:),result_1800.Weights_spa(:));
  r(3,k)=c(3)^2;
  clear result_1800
  load(sprintf('probe-0528-result_%i.mat',N(k)),'result_3600');
  c=corrcoef(result_3600.GT_W(:),result_3600.Weights_spa(:));
  r(4,k)=c(3)^2;
end
rspa=r;
clear result_300 result_600 result_1800 result_3600
%PREPROCESSING PART ++++++++++++++++++++++++++++++++++++++++++++

load('scan-0528-probe.mat')
T=zeros(size(result_50));
r_base=T;
r_spa=r_base;
r_dal=r_base;
for k=1:length(result_50)
  T(k)=result_50(k).T;
  c=corrcoef(result_50(k).GT_W(:),result_50(k).Weights_glm(:));
  r_base(k)=c(3)^2;
  c=corrcoef(result_50(k).GT_W(:),result_50(k).Weights_spa(:));
  r_spa(k)=c(3)^2;
  c=corrcoef(result_50(k).GT_W(:),result_50(k).Weights_dal(:));
  r_dal(k)=c(3)^2;
end
figure
plot(T,r_base,'k.-','LineWidth',2)
hold on
plot(T,r_spa,'k.--','LineWidth',2)
plot(T,r_dal,'k.:','LineWidth',2)
axis([0 3600 0 1])
xlabel('Data available, sec')
ylabel('Recovered variance')
legend({'Regular GLM','Sparse GLM','Sparse & Dale GLM'},'Location','SouthEast')
title('Figure 8: Performance limit vs T for N=50, 100 and 200, orig. spikes')

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
T=[T,3600]; r_base=[r_base,rreg(4,5)]; r_spa=[r_spa,rspa(4,5)];
plot(T,r_base,'k.-','LineWidth',2,'Color',[0.3,0.3,0.3])
hold on
plot(T,r_spa,'k.--','LineWidth',2,'Color',[0.3,0.3,0.3])
plot(T(1:end-1),r_dal,'k.:','LineWidth',2,'Color',[0.3,0.3,0.3])
axis([0 3600 0 1])

%N=200 neurons overlay
plot([300 600 1800 3600],rreg(:,7),'k.-','LineWidth',2,'Color',[0.6,0.6,0.6]);
plot([300 600 1800 3600],rspa(:,7),'k.--','LineWidth',2,'Color',[0.6,0.6,0.6]);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure 9) Hamming inh/exc errors vs T for sparse prior for N=50,100
clear
load('scan-0528-probe.mat')
T=zeros(size(result_50));
r_base=T;
r_spa=r_base;
r_dal=r_base;
for k=1:length(result_50)
  T(k)=result_50(k).T;
  c=dale(result_50(k).Weights_glm);
  r_base(k)=sum(c(1:40)<0)+sum(c(41:50)>0);
  c=dale(result_50(k).Weights_spa);
  r_spa(k)=sum(c(1:40)<0)+sum(c(41:50)>0);
  c=dale(result_50(k).Weights_dal);
  r_dal(k)=sum(c(1:40)<0)+sum(c(41:50)>0);
end
figure
plot(T,r_base/50,'k.-','LineWidth',2)
hold on
plot(T,r_spa/50,'k.--','LineWidth',2)
plot(T,r_dal/50,'k.:','LineWidth',2)
axis([0 3600 0 1])
xlabel('Data available, sec')
ylabel('Recovered variance')
legend({'regular GLM','sparse GLM','sparse & Dale GLM'},'Location','NorthEast')
title('Figure 9: inh/exc id vs T for N=50, 100')

%N=100 neurons
T=[];
r_base=[];
r_spa=[];
r_dal=[];
for k=300:300:2100
  sname=sprintf('scan-0528-result_%i_100.mat',k);
  load(sname);
  eval(sprintf('result_100=result_%i_100;',k));
  T(end+1)=result_100.T;
  c=dale(result_100.Weights_glm);
  r_base(end+1)=sum(c(1:80)<0)+sum(c(81:100)>0);
  c=dale(result_100.Weights_spa);
  r_spa(end+1)=sum(c(1:80)<0)+sum(c(81:100)>0);
  c=dale(result_100.Weights_dal);
  r_dal(end+1)=sum(c(1:80)<0)+sum(c(81:100)>0);
  clear(sprintf('result_%i_100',k));  
end
plot(T,r_base/100,'k.-','LineWidth',2,'Color',[0.3 0.3 0.3])
hold on
plot(T,r_spa/100,'k.--','LineWidth',2,'Color',[0.3 0.3 0.3])
plot(T,r_dal/100,'k.:','LineWidth',2,'Color',[0.3 0.3 0.3])
axis([0 3800 0 0.3])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure 10) Strong/weak solution
figure                          %Figure 10c) weak correlations
load fluor-0610-50-1-data.mat netSim
load fluor-0610-50-1-glm_base0proc1.mat
plot(netSim.weights(:),Weights(:),'.')
hold on
plot([-10 10],[-10 10],'k--','LineWidth',2,'Color',[0.5 0.5 0.5])
axis([-2 2 -2 2])
xlabel('Actual connection weights')
ylabel('Inferred connection weights')
title('Figure 10c: weak correlations solution')

figure                            %Figure 10d) strong correlations
load fluor-0706-50-9-result_bf.mat
plot(result_base.GT_W(:),result_base.Weights_glm(:),'.')
hold on
plot([-10 10],[-10 10],'k--','LineWidth',2,'Color',[0.5 0.5 0.5])
axis([-4 4 -4 4])
xlabel('Actual connection weights')
ylabel('Inferred connection weights')
title('Figure 10d: strong correlations solution')


%Figure X) variation in EPSP time constants, tau
figure
load fluor-0622-25-8-result_iid.mat
plot(result_iid.GT_W(:),result_iid.Weights_glm(:),'r*')
hold on
load fluor-0610-25-1-result_base.mat
plot(result_base.GT_W(:),result_base.Weights_glm(:),'b.')
hold on
plot([-10 10],[-10 10],'k--','LineWidth',2,'Color',[0.5 0.5 0.5])
axis([-2 1.5 -1.25 0.75])
xlabel('Actual connection weights')
ylabel('Inferred connection weights')
title('Figure X: 25% variability in {\bf \tau}')
legend({'Indep. approx.','Orig. spikes'},'Location','NorthWest')