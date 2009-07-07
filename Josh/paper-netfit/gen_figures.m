%Figure 1) MCMC vs IID weights inferrence, N=25 T=600s
clear
load('fluor-0610-25-3-glm_iid0proc1.mat')
W_iid=Weights;
load('fluor-0610-25-3-glm_neal2proc1.mat')
W_mcmc=Weights;
load('fluor-0610-25-3-result_base.mat')
W_GT=result_base.GT_W;

figure
plot(W_GT(:),W_iid(:),'*r')
hold on
plot(W_GT(:),W_mcmc(:),'.')
plot([-2 2],[-2 2],'k--','LineWidth',2)
axis([-2 1.5 -1.25 0.75])
xlabel('Actual connection weights')
ylabel('Inferred connection weights')
title('Figure 1: MCMC vs IID scatter plot')
legend({'IID approximation','Gibbs-MCMC solver'},'Location','NorthWest')

%Figure 2) MCMC vs IID weights inferrence, N=25 T=600s
clear
load('fluor-0610-25-4-glm_iid0proc1.mat')
W_iid=Weights;
load('fluor-0610-25-4-glm_base0proc1.mat')
W_base=Weights;
load('fluor-0610-25-4-result_base.mat')
W_GT=result_base.GT_W;

figure
plot(W_GT(:),W_base(:),'*r')
hold on
plot(W_GT(:),W_iid(:),'.')
plot([-2 2],[-2 2],'k--','LineWidth',2)
axis([-2 1.5 -1.25 0.75])
xlabel('Actual connection weights')
ylabel('Inferred connection weights')
title('Figure 2: Base vs IID scatter plot')
legend({'Baseline','IID approximation'},'Location','NorthWest')


%Figure 3) r^2 vs SNR for N=25, 50, 100
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
plot(gamma/1000,r,'k.-','LineWidth',2);

gamma=[1e3,5e3,10e3,20e3,40e3,80e3];
load('fluor-0610-50-1-data.mat','netSim')
W_GT=netSim.weights;
r=zeros(1,6);
for k=1:6
  try
    load(sprintf('fluor-0610-50-%i-glm_iid0proc1.mat',k))
    W_iid=Weights;
    c=corrcoef(W_GT(:),W_iid(:));
    r(k)=c(3)^2;
  catch
    r(k)=nan;
  end
end
figure(h)
hold on
plot(gamma/1000,r,'k.--','LineWidth',2);

xlabel('Photon budget, Kph/neuron/frame')
ylabel('{\bf\it r}^2')
legend({'N=25, T=600s','N=50, T=800s'},'Location','NorthWest')
title('Figure 3: r^2 vs. calcium imaging noise')

%Figure 5) r^2 vs T for sparse prior for N=50,100
clear
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
ylabel('{\bf\it r}^2')
legend({'GLM regular','GLM with sparse prior',...
  'GLM with sparse and dale prior'},'Location','SouthEast')
title('Figure 4: r^2 vs T for N=50')

%100 neurons
T=[];
r_base=[];
r_spa=[];
r_dal=[];
for k=300:300:3600
  sname=sprintf('scan-0528-result_%i_100.mat',k);
  try 
    load(sname);
  catch
    continue;
  end
  eval(sprintf('result_100=result_%i_100;',k));
  T(end+1)=result_100.T;
  c=corrcoef(result_100.GT_W(:),result_100.Weights_glm(:));
  r_base(end+1)=c(3)^2;
  c=corrcoef(result_100.GT_W(:),result_100.Weights_spa(:));
  r_spa(end+1)=c(3)^2;
  c=corrcoef(result_100.GT_W(:),result_100.Weights_dal(:));
  r_dal(end+1)=c(3)^2;
end
plot(T,r_base,'g.-','LineWidth',2)
hold on
plot(T,r_spa,'g.--','LineWidth',2)
plot(T,r_dal,'g.:','LineWidth',2)
axis([0 3600 0 1])


%Figure 6) histogram of weights
clear
load('fluor-0610-50-5-glm_iid0proc1.mat')
W_glm=Weights;
load('fluor-0610-50-5-spa_base0proc1.mat')
W_spa=Weights;
load('fluor-0610-50-1-data.mat','netSim')
W_GT=netSim.weights;

figure
[h,n]=hist(W_GT(W_GT>0),9);
fp=sum(n.*h)/sum(h);
plot(n,h/sum(h),'r','LineWidth',2);
hold on
[h,n]=hist(W_GT(W_GT<0),9);
fm=sum(n.*h)/sum(h);
plot(n,h/sum(h),'b','LineWidth',2);
[h,n]=hist(W_glm(W_GT>0),9);
fp=fp/(sum(n.*h)/sum(h));
plot(fp*n,h/sum(h),'r--','LineWidth',2);
hold on
[h,n]=hist(W_glm(W_GT<0),9);
fm=fm/(sum(n.*h)/sum(h));
plot(fm*n,h/sum(h),'b--','LineWidth',2);
[h,n]=hist(W_glm(W_GT==0),9);
fo=(fp*sum(W_GT(:)>0)+fm*sum(W_GT(:)<0))/sum(W_GT(:)~=0);
plot(fo*n,h/sum(h),'k--','LineWidth',2);
xlabel('Connection weights')
ylabel('Frequency')
title('Figure 5a: True and inferred distribution of connectivity weights')

figure
[h,n]=hist(W_GT(W_GT>0),9);
fp=sum(n.*h)/sum(h);
plot(n,h/sum(h),'r','LineWidth',2);
hold on
[h,n]=hist(W_GT(W_GT<0),9);
fm=sum(n.*h)/sum(h);
plot(n,h/sum(h),'b','LineWidth',2);
[h,n]=hist(W_spa(W_GT>0),9);
fp=fp/(sum(n.*h)/sum(h));
plot(fp*n,h/sum(h),'r--','LineWidth',2);
hold on
[h,n]=hist(W_spa(W_GT<0),9);
fm=fm/(sum(n.*h)/sum(h));
plot(fm*n,h/sum(h),'b--','LineWidth',2);
[h,n]=hist(W_spa(W_GT==0),9);
fo=(fp*sum(W_GT(:)>0)+fm*sum(W_GT(:)<0))/sum(W_GT(:)~=0);
plot(fo*n,h/sum(h),'k--','LineWidth',2);
xlabel('Connection weights')
ylabel('Frequency')
title('Figure 5b: True and inferred distribution of connectivity weights, sparse')


%Figure 7) r^2 vs N for diff T
clear
N=[10,25,50,75,100,150,200,250];
r=zeros(4,8);
for k=1:8
  clear result_300 result_600 result_1800 result_3600
  load(sprintf('probe-0528-result_%i.mat',N(k)));
  if(exist('result_300'))
    c=corrcoef(result_300.GT_W(:),result_300.Weights_glm(:));
    r(1,k)=c(3)^2;
  else
    r(1,k)=nan;
  end
  if(exist('result_600'))
    c=corrcoef(result_600.GT_W(:),result_600.Weights_glm(:));
    r(2,k)=c(3)^2;
  else
    r(2,k)=nan;
  end
  if(exist('result_1800'))
    c=corrcoef(result_1800.GT_W(:),result_1800.Weights_glm(:));
    r(3,k)=c(3)^2;
  else
    r(3,k)=nan;
  end
  if(exist('result_3600'))
    c=corrcoef(result_3600.GT_W(:),result_3600.Weights_glm(:));
    r(4,k)=c(3)^2;
  else
    r(4,k)=nan;
  end
end

h=figure;
hold on
plot([300 600 1800 3600],r);
axis([0 3600 0 1])
xlabel('Data available, sec')
ylabel('{\bf \it r}^2')
legend('N=10 neurons','N=25 neurons','N=50 neurons','N=75 neurons',...
  'N=100 neurons','N=150 neurons','N=200 neurons',...
  'Location','Southeast')
title('Figure 6a: r^2 vs. N')

%-------------------------------------------------------------
clear
N=[10,25,50,75,100,150,200,250];
r=zeros(4,8);
for k=1:8
  clear result_300 result_600 result_1800 result_3600
  load(sprintf('probe-0528-result_%i.mat',N(k)));
  if(exist('result_300'))
    c=corrcoef(result_300.GT_W(:),result_300.Weights_spa(:));
    r(1,k)=c(3)^2;
  else
    r(1,k)=nan;
  end
  if(exist('result_600'))
    c=corrcoef(result_600.GT_W(:),result_600.Weights_spa(:));
    r(2,k)=c(3)^2;
  else
    r(2,k)=nan;
  end
  if(exist('result_1800'))
    c=corrcoef(result_1800.GT_W(:),result_1800.Weights_spa(:));
    r(3,k)=c(3)^2;
  else
    r(3,k)=nan;
  end
  if(exist('result_3600'))
    c=corrcoef(result_3600.GT_W(:),result_3600.Weights_spa(:));
    r(4,k)=c(3)^2;
  else
    r(4,k)=nan;
  end
end

h=figure;
hold on
plot([300 600 1800 3600],r);
axis([0 max(N) 0 1])
axis([0 3600 0 1])
xlabel('Data available, sec')
ylabel('{\bf \it r}^2')
legend('N=10 neurons','N=25 neurons','N=50 neurons','N=75 neurons',...
  'N=100 neurons','N=200 neurons',...
  'Location','Southeast')
title('Figure 6b: r^2 vs. N sparse')


%Figure 8a) Strong/weak solutions
load fluor-0622-25-9-result_bf.mat
figure
plot(result_base.GT_W(:),result_base.Weights_glm(:),'.')
hold on
plot([-10 10],[-10 10],'k--','LineWidth',2)
axis([-10 10 -6.5 6.5])
xlabel('Actual connection weights')
ylabel('Inferred connection weights')
title('Figure 8a: strong coupling solution')


%Figure 8b) Strong/weak solutions
load fluor-0610-25-1-result_base.mat
figure
plot(result_base.GT_W(:),result_base.Weights_glm(:),'.')
hold on
plot([-10 10],[-10 10],'k--','LineWidth',2)
axis([-2 1.5 -1.25 0.75])
xlabel('Actual connection weights')
ylabel('Inferred connection weights')
title('Figure 8b: weak coupling solution')



%Figure 8c) variable tau
load fluor-0610-25-1-result_base.mat
figure
plot(result_base.GT_W(:),result_base.Weights_glm(:),'r*')
hold on
load fluor-0622-25-8-result_iid.mat
plot(result_iid.GT_W(:),result_iid.Weights_glm(:),'.')
hold on
plot([-10 10],[-10 10],'k--','LineWidth',2)
axis([-2 1.5 -1.25 0.75])
xlabel('Actual connection weights')
ylabel('Inferred connection weights')
title('Figure 8c: variable tau solution')
legend({'Baseline','25% variability in {\bf \tau}'},'Location','NorthWest')

% %Figure 9a) regular glm & sparse overlaid for N=50, T=800
% load fluor-0610-50-5-data.mat netSim
% Wt=full(netSim.weights);
% load fluor-0610-50-5-glm_iid0proc1.mat Weights
% figure
% plot(Wt(:),Weights(:),'r*')
% hold on
% load fluor-0610-50-5-spa_iid0proc1.mat Weights
% plot(Wt(:),Weights(:),'.')
% hold on
% plot([-10 10],[-10 10],'k--','LineWidth',2)
% axis([-2 1.5 -1.25 0.75])
% xlabel('Actual connection weights')
% ylabel('Inferred connection weights')
% title('Figure 9a: GLM regular vs sparse, IID')
% legend({'GLM regular solution','GLM with sparse prior'},'Location','NorthWest')

%Figure 9b) regular glm & sparse overlaid for N=50, T=800
load fluor-0610-50-5-data.mat netSim
Wt=full(netSim.weights);
load fluor-0610-50-5-glm_base0proc1.mat Weights
figure
plot(Wt(:),Weights(:),'r*')
hold on
load fluor-0610-50-5-spa_base0proc1.mat Weights
plot(Wt(:),Weights(:),'.')
hold on
plot([-10 10],[-10 10],'k--','LineWidth',2)
axis([-2 1.5 -1.25 0.75])
xlabel('Actual connection weights')
ylabel('Inferred connection weights')
title('Figure 9b: GLM regular vs sparse, baseline')
legend({'GLM regular solution','GLM with sparse prior'},'Location','NorthWest')
