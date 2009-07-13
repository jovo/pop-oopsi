function run_datafit(fname,cmode,runid,id_proc,N_proc,rndinit)
%LAUNCH script for netfit package-wrapper NETFIT_main
if(nargin<2) cmode=2; end % sampler mode
if(nargin<3) runid=1; end % job id
if(nargin<4) id_proc=1; end% proc-id
if(nargin<5) N_proc=1; end;% proc ##
if(nargin<6) rndinit=3711; end% random seed for reproducibility

load (fname);       %LOAD DATA
for k=1:length(F) F{k}=im2double(F{k}); end

fprintf('In run_datafit %s\n',fname);
% cd /hmt/sardine/hpc/scratch/stats/users/ym2289/glm-netfit

%%%%%%%%%%  SIMULATION NAME  -- unique to the project
N=length(F);        %number of neurons
T=length(F{1});     %trace length

netsim_name='data/s1m1-0713';
netsim_name=[netsim_name,sprintf('-%i-%i',N,runid)];
sync_name=[netsim_name,sprintf('-%i-sync',cmode)];

%INITIALIZE PARAMETERS
lambda=20;          %sparse L1-prior weight
flambda=0;          %find best weight based on actual weights?

if(exist('frame_rate')==1) FR=frame_rate; else FR=30; end%imaging frame rate
if(exist('freq')==1) setFQ=freq; else setFQ=1; end %set skip-freq in PF (oversample)

spkM=1;             %spike-train samples from spk-sampler, for GLM
tmin=1;             %min coupling time-depth, >1 
tmax=1;             %max coupling time-depth

Tp=tmax-tmin+1;     %couplings temporal depth

setH=1;             %h-dimensions in PF

flgSparse=1;        %compute sparse solution?
flgDale=1;          %compute dale solution?

netSim.K=1;         %some dummies for compatibility with NETFIT_main
netSim.FR=1/FR;
netSim.dt=1/FR;
netSim.weights=zeros(N,N);
n_GT=cell(N,1);


% CUSTOM INITIALS FOR particle filter parameters
cP.rate     = 0;                    % expected rate
cP.k        = [cP.rate,zeros(1,0)]';% linear filter
cP.tau      = 0.5;                  % calcium decay time constant (sec)
cP.A        = 150;                  % jump size ($\mu$M)
cP.sig      = 25;                   % standard deviation of noise (\mu M)
cP.C_0      = 25;                   % baseline [Ca++]
cP.n        = 1;                    % hill equation exponent
cP.k_d      = 100^cP.n;              % hill coefficient
cP.alpha    = 1;                    % F_max
cP.beta     = 0;                    % F_min
cP.gamma    = 1e-4;                 % scaled variance
cP.zeta     = 4*cP.gamma;            % constant variance

cP.C_init   = cP.C_0;                % initial [Ca++]
cP.tau_c    = cP.tau;
cP.sigma_c  = cP.sig;
cP.a        = 1/FR/cP.tau_c;
cP.lam      = T/(cP.rate*cP.A)/FR;% expected jump size ber time bin

if setH>0                                        % if there are spike history terms
    cP.omega = -1;                               % weight
    cP.tau_h = 0.01;                             % time constant
    cP.tau_h=1/FR/(1-exp(-1/FR/cP.tau_h));       % correction -- large dt
    cP.sigma_h = 0.01;                           % stan dev of noise
end



flgData=1;          %tell NETFIT_main that this is running with data (no sim)


switch(runid)     %RUNID SPECIFIC ADJ
  case 1
  otherwise
end
switch(cmode)       %CHOOSE INFERENCE ALGORITHM
  case 1
    error('Base sampler is impossible wiht real data'); 
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
