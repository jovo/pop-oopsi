function run_netfit(N,cmode,runid,id_proc,N_proc,rndinit)
%LAUNCH script for netfit package-wrapper NETFIT_main
if(nargin<1) N=25; end    %neurons number
if(nargin<2) cmode=1; end %selection mode
if(nargin<3) runid=5; end %selection of other run par (N,T, etc)
if(nargin<4) id_proc=1; end
if(nargin<5) N_proc=1; end;
if(nargin<6) rndinit=3711; end%random seed

fprintf('In run_netfit %i\n',cmode);
% cd /hmt/sardine/hpc/scratch/stats/users/ym2289/glm-netfit

%%%%%%%%%%  SIMULATION NAME  -- unique to the project
netsim_name='data/fluor-0706';
netsim_name=[netsim_name,sprintf('-%i-%i',N,runid)];
sync_name=[netsim_name,sprintf('-%i-sync',cmode)];

%INITIALIZE PARAMETERS
% N=10;               %number of neurons
T=250;              %trace length

lambda=0;           %sparse L1-prior weight
flambda=0;          %find best weight based on actual weights?
fdale=0;            %dale prior enacted [0/1]?

FR=66;              %imaging frame rate
SP=max(0.1,2/N);    %connections sparseness

spkM=1;             %spike-train samples from spk-sampler, for GLM
tmin=1;             %min coupling time-depth, >1 
tmax=1;             %max coupling time-depth

Tp=tmax-tmin+1;     %couplings temporal depth

setH=1;             %h-dimensions in PF
setFQ=1;            %set skip-freq in PF (oversample)

flgSparse=1;        %compute sparse solution?
flgDale=1;          %compute dale solution?

if(N<10)            %ASSIGN TIME
  T=600;
elseif(N<=25)           
  T=600;
elseif(N<=50)
  T=600;
elseif(N<=100)
  T=600;
else
  T=600;
end
switch(runid)     %RUNID SPECIFIC ADJ
  case 1
    gamma=1e3;
  case 2
    gamma=5e3;
  case 3
    gamma=10e3;
  case 4
    gamma=20e3;
  case 5
    gamma=40e3;
  case 6
    gamma=80e3;  
  case 7
    FR=33;
    setFQ=6;
  case 8
    varTau=0.25;
  case 9
    CM=1;
  case 50
    N=50;
    T=600;
  otherwise
    N=10;
    T=250;
end
switch(cmode)       %CHOOSEINFERENCE ALGORITHM
  case 1
    mode='base';    %spike sampler base/iid/neal/indep-gibbs
    spkM=1;         %adjust glm-sample size
  case 2
    mode='iid';
    spkM=25;        %adjust glm-sample size
  case 3
    mode='neal';
    spkM=10;        %adjust glm-sample size
  case 4
    mode='gibbs';   %gibbs using indep-PF inside
    spkM=25;        %THIS IS A HACK NOT RECOMMEND TO USE
  otherwise
    error('Unknown main sampler type');    
end

D=ceil(N/N_proc);
for i=1:N_proc nrange{i}=(1+(i-1)*D):min(N,i*D); end

NETFIT_main         %THIS ACTUALLY DOES THINGS
