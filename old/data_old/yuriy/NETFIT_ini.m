% this script generates stochastic network
clear, close, clc, fprintf('\nnetSim.SIMULATE\n')

netsim_name='netSim0315N50S5678mod';
short=0;                        %use 1 if adjusting params for firing rate

rndinit=1234;                   %initialization for rnd num generator                           
rand('state',rndinit);          %for repeatability
randn('state',rndinit);
                          
netSim.N=50;                    %number of cells
netSim.SP=0.1;                  %sparseness of connections
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

%DEFINED TEST
% netsim_name='netSim0315N50S5678short';
% netSim.T=1.2e5;                 %simulation duration, steps
% netSim.dt=0.005;                %simulation time-step
% netSim.FR=0.005;                %sampling time step
% netSim.hepsp=1;                 %EPSP profile
% netSim.hipsp=1;                 %IPSP profile
% netSim.hrefr=1;                 %REFR profile


%%GENERATE CONNECTIVITY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sparse matrix::exp-distr weights::matrix-order POST-PRE|(i,j)
rate=netSim.rate(2);                %mean firing rate
identie=false(netSim.N,1);          %I/E identity
Ninh=round(netSim.N*netSim.balance);%number of inhibitory neurons

netSim.weights=full(sprand(netSim.N,netSim.N,netSim.SP));
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


%%GENERATE POPULATION OF CELLS/calcium models
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
randtrn=inline('max(A(1),A(2)+randn*A(3))','A');%truncated Gaussian def
netSim.cells=[]; 
rate=[];
for k=1:netSim.N
  c=[];
  
  c.A=randtrn(netSim.A);
  c.C_0=randtrn(netSim.C_0)*c.A;  
  c.sigma=randtrn(netSim.sigma)*c.A;  
  c.tau=randtrn(netSim.tau);
  c.n=netSim.n;
  c.k_d=netSim.k_d;
  
  c.alpha=randtrn(netSim.alpha);
  c.beta=randtrn(netSim.beta);
  c.gamma=1/randtrn(netSim.gamma);
  c.zeta=randtrn(netSim.zeta)*c.gamma;
  
  c.omega=randtrn(netSim.omega);  
  if(identie(k))                                %admit two rates for subpops
    c.rate=randtrn(netSim.rate);
  else
    c.rate=randtrn(netSim.irate);
  end
  
  if(isempty(netSim.cells)) netSim.cells=c; else netSim.cells(k)=c; end
end

%%GENERATE FIRING PATTERNS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(short) netSim.T=round(5/netSim.dt); end    %if aiming, use 5 seconds
fprintf('Generating firing patterns...\n');
Np=netSim.N;
Nexc=sum(identie);
Ninh=sum(~identie);
n=cell(Np,1);                                 %spikes
C=zeros(Np,floor(netSim.T/netSim.K));         %calcium traces
H=zeros(Np,length(netSim.hepsp));             %spike histories
nbuf=zeros(Np,1);                             %buffer spikes
Cbuf=zeros(Np,1); Cbuf(:)=[netSim.cells.C_0]; %buffer calcium
for t=1:netSim.T                              %for all times
  Ji=zeros(Np,1);                             %injected currents:  
  Ji(identie)=sum(H(identie,:) .* repmat(netSim.hepsp,Nexc,1),2);%excitatory  
  Ji(~identie)=sum(H(~identie,:) .* repmat(netSim.hipsp,Ninh,1),2);%inhibitory
  Jj=sum(H .* repmat(netSim.hrefr,Np,1),2);   %refractory currents
  
  for k=1:Np
    %SPIKE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    f=netSim.invtransf(netSim.cells(k).rate); %spontaneous contribution
    
    g=-Jj(k) .* netSim.cells(k).omega;        %refractory contribution
    
    g=g + sum(netSim.weights(k,:)' .* Ji);    %PSP contribution    
    
    nbuf(k)=rand>exp(-netSim.transf(f+g)*netSim.dt);%Bernouli spike
    if(nbuf(k)>0) n{k}=[n{k},t]; end

    %CALCIUM
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    A     = netSim.cells(k).A;                %aliases
    sigma = netSim.cells(k).sigma;
    tau   = netSim.cells(k).tau;
    C_0   = netSim.cells(k).C_0;
    
    a     = netSim.dt/tau;  
    Cbuf(k)=max(0,(1-a)*Cbuf(k) + a*C_0 + A*nbuf(k) + sigma*sqrt(netSim.dt)*randn);
  end  
  H=[H(:,2:end),nbuf];                        %update spike-histories
    
  if(mod(t,netSim.K)==0) C(:,t/netSim.K)=Cbuf; end%record sample
    
  if(mod(t,round(netSim.T/25))==0) fprintf('.'); pause(0.01); end
end
fprintf('\n');

%tot ## of spikes subpop::exc
ntote=0; for k=find(identie)' ntote=ntote+length(n{k}); end
%tot ## of spikes subpop::inh
ntoti=0; for k=find(~identie)' ntoti=ntoti+length(n{k}); end
fprintf('Done, mean firing rate is %.3g/%.3g Hz\n',...
  ntote/sum(identie)/(netSim.T*netSim.dt),...
   ntoti/sum(~identie)/(netSim.T*netSim.dt));
if(short) return; end

    
%%GENERATE Fluorescence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Generating fluorescence...\n');
F=zeros(size(C));
for k=1:netSim.N
  alpha   = netSim.cells(k).alpha;            %aliases
  beta    = netSim.cells(k).beta;
  gamma   = netSim.cells(k).gamma;
  zeta    = netSim.cells(k).zeta;  
  P.n     = netSim.n;
  P.k_d   = netSim.k_d;
  
  S=netSim.hill(P,C(k,:));
  F(k,:)=max(0,alpha*S+beta+sqrt(gamma*S+zeta).*randn(size(S)));
end

fprintf('Saving data...\n');
save([netsim_name,'.mat'],'netSim','n','C','F');
