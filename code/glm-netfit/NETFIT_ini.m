%SCRIPT TO GENERATE CALCIUM IMAGING DATA FOR 
%A POPULATION OF COUPLED GLM NEURONS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INITIALIZATION
short=0;                            %use 1 if adjusting firing rate

if(exist('T')~=1)  T = 600; end     %recording time
if(exist('N')~=1)  N = 50; end      %# of neurons
if(exist('SP')~=1) SP = 0.1; end    %sparseness of matrix
if(exist('FR')~=1) FR = 66; end     %frame rate imaging
if(exist('G')~=1)  G = 4e4; end     %shot-noise multiple
if(exist('CM')~=1) CM = 0; end      %weak/strong coupling
if(exist('rndinit')~=1) rndinit = 3467; end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rand('state',rndinit);           %for repeatability
randn('state',rndinit);

netSim.rnd=rndinit;              %random seed
netSim.N=N;                      %number of cells
netSim.SP=SP;                    %sparseness of connections
netSim.scaleE=0.025;             %excitatiry abs scale
netSim.scaleI=0.1;               %inhibitory abs scale
netSim.balance=0.2;              %fraction of inhibitory neurons
netSim.dt=0.001;                 %simulation time-step
netSim.FR=1/FR;                  %sampling time step
netSim.T=ceil(T/netSim.dt);      %simulation duration, steps
netSim.transf='exp';             %nonlinearity, FIY

netSim.A=[20 40 12];             %Ca-spike intensity, min, mean & STD
netSim.sigma=[0.05 0.35 0.1];    %Ca-noise/spike ratio, min, mean & STD
netSim.tau=[0.15 0.25 0.05];     %Ca-decay constant, min, mean & STD
netSim.C_0=[0.05 0.30 0.1];      %Ca-baseline/spike ratio, min, mean & STD
netSim.n=1;                      %Hill-profile exponent
netSim.k_d=100;                  %Hill-profile scale

netSim.alpha=[0 1 0];            %fluorescence scaling, min, mean & STD
netSim.beta =[0 0 0];            %fluorescence offset, min mean & STD
netSim.gamma=G*[0.3 1 0.5];      %fluorescence generation, shot noise, min, mean & STD
netSim.zeta=[0 4 0];             %fluorescence generation, background noise
netSim.hill=inline('C.^(P.n)./(C.^(P.n)+P.k_d)','P','C');%nonlinearity

netSim.rate =[0 5 0];           %spontaneous ignition rate, min, mean & STD
netSim.irate=[0 5 0];           %spontaneous ignition rate for inh neurons, same
netSim.omega=[0 100 0];         %refraction strength, min, mean & STD

netSim.tauexc=10*[0.5 1 0.];  %excitatory time scale, decay constant
netSim.tauinh=20*[0.5 1 0.];  %inhibitory time scale, decay constant
netSim.tauref=10*[0.5 1 0.];  %refractory time scale, decay constant

netSim.K=round(netSim.FR/netSim.dt);%sampling/sim time-steps ratio


%%%%%%%%% ADJUSTMENTS ----
if(N>100) netSim.scaleI=0.2; end%stronger inhibition is necessary for larger N
if(N>400) netSim.scaleI=0.25; end

if(CM==1)                       %strong coupling subpopulation
  netSim.SP=0.2;
  netSim.fractCM=0.2;
  netSim.scaleES=15.0;
  netSim.scaleIS=30.0;
  netSim.rate =[0 1 0];         %lower base rate
  netSim.irate=[0 1 0]; 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CONNECTIVITY MATRIX -- sparse::exp-distr::order POST-PRE|(i,j)
dT=0.01;                            %EPSP typical time-scale
rate=netSim.rate(2);                %mean firing rate
signs=zeros(netSim.N,1);            %neuron sign, I/E
Ninh=round(netSim.N*netSim.balance);%number of inhibitory neurons

netSim.weights=sprand(netSim.N,netSim.N,netSim.SP);

%excitatory neurons -- log-rates
Zmin=1e-6;                          %truncation for log
idx=1:netSim.N-Ninh; signs(idx)=1;
W=netSim.weights(:,idx); ind=find(W>0);
z=netSim.scaleE*exprnd(1,length(ind),1);
if(CM==1)
  idxCM=find(rand(size(ind))<netSim.fractCM);
  z(idxCM)=netSim.scaleES*exprnd(1,length(idxCM),1);
  z=log(1+z/(rate*dT));
else
  z=log(-log(max(Zmin,exp(-rate*dT)-z))/(rate*dT));
end
W(ind)=z; netSim.weights(:,idx)=W;

%inhibitory neurons -- log-rates
Zmin=1e-6;                        %truncation for log
idx=netSim.N-Ninh+1:netSim.N; signs(idx)=-1;
W=netSim.weights(:,idx); ind=find(W>0);
z=netSim.scaleI*exprnd(1,length(ind),1);
if(CM==1)
  idxCM=find(rand(size(ind))<netSim.fractCM);
  z(idxCM)=netSim.scaleIS*exprnd(1,length(idxCM),1);  
  z=-log(1+z/(rate*dT));
else
  z=-log(-log(max(Zmin,exp(-rate*dT)-z))/(rate*dT));
end
W(ind)=z; netSim.weights(:,idx)=W;

%remove self-coupling, refractory is Omega, not W
for k=1:netSim.N netSim.weights(k,k)=0; end
fprintf('W-scale: %g;%g\n',...
  mean(abs(full(netSim.weights(netSim.weights>0)))),...
   mean(abs(full(netSim.weights(netSim.weights<0)))));

netSim.signs=signs;               %take note of signs


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GENERATE POPULATION OF CELLS/calcium models
randtrn=inline('max(A(1),A(2)+randn*A(3))','A');%truncated Gaussian def
netSim.cells=[];
rate=[];
for k=1:netSim.N
  c=[];

  c.A=randtrn(netSim.A);                %calcium jump
  c.C_0=randtrn(netSim.C_0)*c.A;        %calcium base
  c.sigma=randtrn(netSim.sigma)*c.A;    %calcium noise
  c.tau=randtrn(netSim.tau);            %calcium decay
  c.n=netSim.n;                         %fluorescence Hill-exp
  c.k_d=netSim.k_d;                     %fluorescence Hill-K_d

  c.alpha=randtrn(netSim.alpha);        %fluorescence alpha
  c.beta=randtrn(netSim.beta);          %fluorescence beta
  c.gamma=1/randtrn(netSim.gamma);      %fluorescence gamma
  c.zeta=randtrn(netSim.zeta)*c.gamma;  %fluorescence zeta
  
  c.sigma_j=0;                          %intrinsic noise, if any
  c.sigma_h=0;                          %noise in currents, if any

  c.omega=randtrn(netSim.omega);        %refractory Omega
  if(signs(k)>0)                        %admit two rates for I/E
    c.rate=randtrn(netSim.rate);        %base firing rate
    rate=[rate,c.rate];
  else
    c.rate=randtrn(netSim.irate);
    rate=[rate,c.rate];    
  end
  c.b0=log(c.rate);                     %base firing rate, log-units
  
  %FEEDBACK PROFILES/should be same length/
  t=(120:-1:0);                         %profile length
  tauexc=randtrn(netSim.tauexc);
  tauinh=randtrn(netSim.tauinh);
  tauref=randtrn(netSim.tauref);
  if(signs(k)>0)
    h=exp(-t/tauexc)-exp(-t/1);
    c.htau=tauexc;
    c.h=h(:)/max(h);                    %EPSP profile
  else  
    h=exp(-t/tauinh)-exp(-t/1);
    c.htau=tauinh;
    c.h=h(:)/max(h);                     %IPSP profile
  end
  h=exp(-t/tauref);
  c.rtau=tauref;
  c.hrefr=h(:)/max(h);                  %REFR profile  

  if(isempty(netSim.cells)) netSim.cells=c; else netSim.cells(k)=c; end
end


%FIRING PATTERNS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(short) netSim.T=round(5/netSim.dt); end    %if aiming, use 5 seconds
fprintf('Generating firing patterns...\n');
Np=netSim.N;
Nexc=sum(signs>0);
Ninh=sum(signs<0);
n=cell(Np,1);                                 %spikes
H=false(Np,length(netSim.cells(1).h));        %contianer spike-history
nbuf=false(Np,1);                             %buffer spike-history
for t=1:netSim.T                              %for all times
  Ji=zeros(Np,1);                             %injected currents:
  Ji=sum(H.*[netSim.cells.h]',2); 
  Jj=sum(H.*[netSim.cells.hrefr]',2);         %refractory currents

  Ji=Ji+[netSim.cells.sigma_h]'.*randn(Np,1); %currents noise
  
  %SPIKE
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  f=[netSim.cells.b0]';                       %spontaneous contribution
  
  g=-Jj.*[netSim.cells.omega]';               %refractory contribution
  
  g=g+netSim.weights*Ji;                      %PSP contribution
  
  g=g+[netSim.cells.sigma_j]'.*randn(Np,1);   %intrinsic noise
  
  nbuf=rand(Np,1)>exp(-exp(f+g)*netSim.dt);   %Bernouli spike

  if(~isempty(find(nbuf,1))) for k=find(nbuf)' n{k}=[n{k},t]; end; end

  H=[H(:,2:end),nbuf];                        %shift spike-history container

  if(mod(t,round(netSim.T/25))==0) fprintf('.'); pause(0.01); end
end
fprintf('\n');
%tot ## of spikes subpop::exc
ntote=0; for k=find(signs>0)' ntote=ntote+length(n{k}); end
%tot ## of spikes subpop::inh
ntoti=0; for k=find(signs<0)' ntoti=ntoti+length(n{k}); end
fprintf('Base rate %.3g mean firing rates are %.3g/%.3g Hz\n',mean(rate),...
  ntote/sum(signs>0)/(netSim.T*netSim.dt),...
  ntoti/sum(signs<0)/(netSim.T*netSim.dt));
if(short) return; end


%CALCIUM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Generating calcium...\n');
normC=zeros(Np,1);
maxC=65535;                                   %uint16 conv factor
C=cell(Np,1);                                 %calcium traces
A     = [netSim.cells.A]';                    %aliases
sigma = [netSim.cells.sigma]';
tau   = [netSim.cells.tau]';
C_0   = [netSim.cells.C_0]';
a     = netSim.dt./tau;

Cbuf=C_0;
nbuf=false(Np,netSim.T);
for k=1:Np nbuf(k,n{k})=1; end
Ctmp=zeros(Np,ceil(netSim.T/netSim.K));
for t=1:netSim.T
  Cbuf=(1-a).*Cbuf + a.*C_0 + A.*nbuf(:,t) + sigma.*sqrt(netSim.dt).*randn(Np,1);
  Cbuf=max(0,Cbuf);                           %can't get smaller than 0
  
  if(mod(t,netSim.K)==0) Ctmp(:,t/netSim.K)=Cbuf; end%record calcium sample
end
normC=max(Ctmp,[],2);
for k=1:Np C{k}=uint16(round(Ctmp(k,:)/normC(k)*maxC)); end


%FLUORESCENCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Generating fluorescence...\n');
F=cell(Np,1);
for k=1:netSim.N
  alpha   = netSim.cells(k).alpha;            %aliases
  beta    = netSim.cells(k).beta;
  gamma   = netSim.cells(k).gamma;
  zeta    = netSim.cells(k).zeta;
  P.n     = netSim.n;
  P.k_d   = netSim.k_d;

  S=netSim.hill(P,normC(k)*double(C{k})/maxC);
  F{k}=max(0,alpha*S+beta+sqrt(gamma*S+zeta).*randn(size(S)));
end

fprintf('Saving data...\n');
n_GT=n;
save(fname,'netSim','n_GT','maxC','normC','C','F');

%clean the mess
clear t h dt rate signs Ninh Zmin idx z W k randtrn c Np Nexc H
clear nbuf Cbuf Ctmp A sigma tau C_0 a  Ji Jj f g ntote ntoti n
clear alpha beta gamma zeta P S
clear SP FR B T G varTau CM rndinit
