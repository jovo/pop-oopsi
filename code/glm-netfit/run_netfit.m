function run_netfit(N,cmode,runid,id_proc,N_proc,rndinit)
%LAUNCH script for netfit package-wrapper NETFIT_main
if(nargin<1) N=25; end    %neurons number
if(nargin<2) cmode=1; end %selection mode
if(nargin<3) runid=5; end %selection of other run par (N,T, etc)
if(nargin<4) id_proc=1; end
if(nargin<5) N_proc=1; end;
if(nargin<6) rndinit=3711; end%random seed

% cd /hmt/sardine/hpc/scratch/stats/users/ym2289/glm-netfit
fprintf('In run_netfit %i\n',cmode);

%%%%%%%%%%  SIMULATION NAME  -- unique to the project
netsim_name='data/fluor-0717';
netsim_name=[netsim_name,sprintf('-%i-%i',N,runid)];
sync_name  =[netsim_name,sprintf('-%i-sync',cmode)];

%INITIALIZE PARAMETERS
T=600;              %trace length

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

switch(runid)       %RUNID SPECIFIC ADJ
  case 1            % different noise amounts
    G=1e3;
  case 2
    G=5e3;
  case 3
    G=10e3;
  case 4
    G=20e3;
  case 5
    G=40e3;
  case 6
    G=80e3;  
  case 9            % strong coupling model
    CM=1;
  case 11            % different noise amounts
    G=1e3; FR=33;
  case 12
    G=5e3; FR=33;
  case 13
    G=10e3; FR=33;
  case 14
    G=20e3; FR=33;
  case 15
    G=40e3; FR=33;
  case 16
    G=80e3; FR=33;  
  case 21            % different noise amounts
    G=1e3; FR=15;
  case 22
    G=5e3; FR=15;
  case 23
    G=10e3; FR=15;
  case 24
    G=20e3; FR=15;
  case 25
    G=40e3; FR=15;
  case 26
    G=80e3; FR=15;
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
